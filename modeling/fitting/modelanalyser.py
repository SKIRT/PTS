#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelanalyser Contains the FitModelAnalyser class, used for analysing the goodness
#  of the radiative transfer model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import tables, time
from ...core.basics.table import SmartTable
from .component import load_fitting_run

# -----------------------------------------------------------------

class FluxDifferencesTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FluxDifferencesTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Instrument", str, None, "Instrument")
        self.add_column_info("Band", str, None, "Band")
        self.add_column_info("Flux difference", float, None, "Flux difference")
        self.add_column_info("Relative difference", float, None, "Relative flux difference")
        self.add_column_info("Chi squared term", float, None, "Chi squared term")

    # -----------------------------------------------------------------

    def add_entry(self, instrument, band, difference, relative_difference, chi_squared_term):

        """
        This function ...
        :param instrument:
        :param band:
        :param difference:
        :param relative_difference:
        :param chi_squared_term:
        :return:
        """

        self.add_row([instrument, band, difference, relative_difference, chi_squared_term])

# -----------------------------------------------------------------

class FitModelAnalyser(FittingComponent):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(FitModelAnalyser, self).__init__(config, interactive)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The fitting run
        self.fitting_run = None

        # The name of the generation
        self.generation_name = None

        # The flux calculator
        self.flux_calculator = None

        # The weights given to each band for the calculation of the chi squared
        self.weights = None

        # The flux differences table
        self.differences = None

        # The calculated chi squared value
        self.chi_squared = None

        # The chi squared table
        self.chi_squared_table = None

    # -----------------------------------------------------------------

    @classmethod
    def for_simulation(cls, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Create the instance
        analyser = cls()

        # Set the modeling path as the working path for this class
        analyser.config.path = simulation.analysis.modeling_path

        # Set the task
        analyser.simulation = simulation

        # Return the instance
        return analyser

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 4. Calculate the differences
        self.calculate_differences()

        # 5. Calculate the chi squared for this model
        self.calculate_chi_squared()

        # 6. Load the chi squared table
        self.load_chi_squared_table()

        # 7. Update the status of the generation if necessary
        self.update_generation()

        # 8. Write
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the fit model analyser ...")

        # Set the attributes to default values
        self.simulation = None
        self.generation_name = None
        self.flux_calculator = None
        self.differences = None
        self.chi_squared = None
        self.chi_squared_table = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FitModelAnalyser, self).setup(**kwargs)

        # Get a reference to the flux calculator
        simulationanalyser = kwargs.pop("simulation_analyser")
        flux_calculator = simulationanalyser.basic_analyser.flux_calculator

        # Make a local reference to the flux calculator
        if flux_calculator is None:
            raise RuntimeError("No ObservedFluxCalculator found; the calculate_observed_fluxes flag must be enabled on "
                               "each simulation that is part of the radiative transfer modeling")
        self.flux_calculator = flux_calculator

        # Determine the generation directory
        generation_path = fs.directory_of(self.simulation.base_path)

        # Set the name of the generation
        self.generation_name = fs.name(fs.directory_of(self.simulation.base_path))

        # Set the name of the fitting run
        fitting_run_name = fs.name(fs.directory_of(fs.directory_of(generation_path)))

        print(generation_path)
        print(self.generation_name)
        print(fitting_run_name)
        print(self.config.path)
        exit()

        # Load the fitting run
        self.fitting_run = load_fitting_run(self.config.path, fitting_run_name)

        # Load the weights table
        #self.weights = tables.from_file(self.weights_table_path, fix_floats=True) # For some reason, the weights are parsed as strings instead of floats (but not from the command line!!??)
        self.weights = tables.from_file(self.fitting_run.weights_table_path)

        # Initialize the differences table
        self.differences = FluxDifferencesTable()

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences between the observed and simulated SED ...")

        # In the flux-density tables derived from the simulation (created by the ObservedFluxCalculator object),
        # search the one corresponding to the "earth" instrument
        mock_sed_name = self.object_name + "_earth"
        if mock_sed_name not in self.flux_calculator.mock_seds: raise RuntimeError("Could not find a mock observation SED for the 'earth' instrument")

        # Get the mock SED
        mock_sed = self.flux_calculator.mock_seds[mock_sed_name]

        # Loop over the entries in the fluxdensity table (SED) derived from the simulation
        for i in range(len(mock_sed)):

            # Get instrument, band and flux density
            instrument = mock_sed["Instrument"][i]
            band = mock_sed["Band"][i]
            fluxdensity = mock_sed["Photometry"][i]

            # Find the corresponding flux in the SED derived from observation
            observed_fluxdensity = self.observed_sed.photometry_for_band(instrument, band, unit="Jy").value

            # Find the corresponding flux error in the SED derived from observation
            observed_fluxdensity_error = self.observed_sed.error_for_band(instrument, band, unit="Jy").average.to("Jy").value

            # If no match with (instrument, band) is found in the observed SED
            if observed_fluxdensity is None:
                log.warning("The observed flux density could not be found for the " + instrument + " " + band + " band")
                continue

            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # Find the index of the current band in the weights table
            index = tables.find_index(self.weights, key=[instrument, band], column_name=["Instrument", "Band"])

            if index is None: continue # Skip this band if a weight is not found

            # Get the weight
            weight = self.weights["Weight"][index] # apparently, this is a string, so parsing the table went wrong ...
            weight = float(weight)

            # Calculate the chi squared term
            chi_squared_term = weight * difference ** 2 / observed_fluxdensity_error ** 2

            # Add an entry to the differences table
            self.differences.add_entry(instrument, band, difference, relative_difference, chi_squared_term)

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the chi squared value for this model ...")

        # Calculate the degrees of freedom
        dof = len(self.observed_sed) - 3. - 1.  # number of data points - number of fitted parameters - 1

        # The (reduced) chi squared value is the sum of all the terms (for each band),
        # divided by the number of degrees of freedom
        self.chi_squared = np.sum(self.differences["Chi squared term"]) / dof

        # Debugging
        log.debug("Found a (reduced) chi squared value of " + str(self.chi_squared))

    # -----------------------------------------------------------------

    def load_chi_squared_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the chi squared table ...")

        # Open the table
        self.chi_squared_table = self.fitting_run.chi_squared_table_for_generation(self.generation_name)

    # -----------------------------------------------------------------

    def update_generation(self):

        """
        This function ...
        :return:
        """

        # Find the index in the table for this generation
        index = tables.find_index(self.fitting_run.generations_table, self.generation_name, "Generation name")

        # Get the number of simulations for this generation
        nsimulations = self.fitting_run.generations_table["Number of simulations"][index]

        # Get the number of entries in the chi squared table
        nfinished_simulations = len(self.chi_squared_table)

        # If this is the last simulation
        if nsimulations == nfinished_simulations + 1:

            # Update the generations table
            self.fitting_run.generations_table.set_finishing_time(self.generation_name, time.timestamp())
            self.fitting_run.generations_table.save()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the flux differences
        self.write_differences()

        # Write the chi-squared value
        self.write_chi_squared()

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with the flux density differences for the current model ...")

        # Determine the path to the differences table
        path = fs.join(self.simulation.analysis.misc.path, "differences.dat")

        # Save the differences table
        self.differences.saveto(path)

    # -----------------------------------------------------------------

    def write_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding the chi squared value for the current model to the chi squared data file ...")

        # Add entry
        self.chi_squared_table.add_entry(self.simulation.name, self.chi_squared)

        # Save the table
        self.chi_squared_table.save()

# -----------------------------------------------------------------
