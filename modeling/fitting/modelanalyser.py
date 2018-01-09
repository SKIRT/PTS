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
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import tables, time
from ...core.basics.table import SmartTable
from ...core.filter.filter import parse_filter_from_instrument_and_band
from .tables import WeightsTable
from ...core.simulation.remote import get_simulation_id, get_simulation_for_host
from ...core.launch.analyser import SimulationAnalyser

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

    def instruments(self):

        """
        Thisf ucntion ...
        :return:
        """

        return list(self["Instrument"])

    # -----------------------------------------------------------------

    def bands(self):

        """
        This function ...
        :return:
        """

        return list(self["Band"])

    # -----------------------------------------------------------------

    def filters(self, sort=False):

        """
        This function ...
        :param sort:
        :return:
        """

        # Get the filter objects
        filters = [parse_filter_from_instrument_and_band(instrument, band) for instrument, band in zip(self.instruments(), self.bands())]

        # Sort
        if sort: return list(sorted(filters, key=lambda fltr: fltr.wavelength.to("micron").value))
        else: return filters

    # -----------------------------------------------------------------

    def wavelengths(self, sort=False, add_unit=True, unit="micron"):

        """
        This function ...
        :param sort:
        :param add_unit:
        :param unit:
        :return:
        """

        # Return
        if add_unit: return [fltr.wavelength.to(unit) for fltr in self.filters(sort=sort)]
        else: return [fltr.wavelength.to(unit).value for fltr in self.filters(sort=sort)]

    # -----------------------------------------------------------------

    def index_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return tables.find_index(self, [fltr.instrument, fltr.band], ["Instrument", "Band"])

    # -----------------------------------------------------------------

    def filter_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        instrument = self["Instrument"][index]
        band = self["Band"][index]
        return parse_filter_from_instrument_and_band(instrument, band)

    # -----------------------------------------------------------------

    def wavelength_for_index(self, index):

        """
        This function ....
        :param index:
        :return:
        """

        return self.filter_for_index(index).wavelength

    # -----------------------------------------------------------------

    def difference_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get the index
        index = self.index_for_filter(fltr)

        # Return the difference
        return self["Flux difference"][index]

    # -----------------------------------------------------------------

    def relative_difference_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get the index
        index = self.index_for_filter(fltr)

        # Return the difference
        return self["Relative difference"][index]

    # -----------------------------------------------------------------

    def chi_squared_term_for_filter(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        # Get the index
        index = self.index_for_filter(fltr)

        # Return the chi squared term
        return self["Chi squared term"] if not self["Chi squared term"].mask[index] else None

    # -----------------------------------------------------------------

    def add_entry(self, instrument, band, difference, relative_difference, chi_squared_term=None, sort=False):

        """
        This function ...
        :param instrument:
        :param band:
        :param difference:
        :param relative_difference:
        :param chi_squared_term:
        :param sort:
        :return:
        """

        self.add_row([instrument, band, difference, relative_difference, chi_squared_term])

        # Sort the table
        if sort: self.sort()

    # -----------------------------------------------------------------

    def sort(self):

        """
        This function ...
        :return:
        """

        # Add a column with the wavelengths
        self["Wavelength"] = self.wavelengths(add_unit=False)
        self["Wavelength"].unit = "micron"

        # Sort the table on the wavelength column
        super(FluxDifferencesTable, self).sort("Wavelength")

        # Delete the wavelength column
        self.remove_column("Wavelength")

    # -----------------------------------------------------------------

    def add_from_filter_and_fluxes(self, fltr, flux, reference_flux, sort=False):

        """
        This function ...
        :param fltr:
        :param flux:
        :param reference_flux:
        :param sort:
        :return:
        """

        # Calculate the difference
        difference = flux - reference_flux
        relative_difference = difference / reference_flux

        # Add
        self.add_entry(fltr.instrument, fltr.band, difference, relative_difference, sort=sort)

# -----------------------------------------------------------------

class SEDFitModelAnalyser(FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SEDFitModelAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The fitting run
        self.fitting_run = None

        # The weights given to each band for the calculation of the chi squared
        self.weights = None

        # The flux differences table
        self.differences = None

        # The calculated chi squared value
        self.chi_squared = None

        # The chi squared table
        self.chi_squared_table = None

        # The mock observed SEDs
        self.mock_seds = None

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

        # 2. Calculate the differences
        self.calculate_differences()

        # 3. Calculate the chi squared for this model
        self.calculate_chi_squared()

        # 4. Load the chi squared table
        self.load_chi_squared_table()

        # 5. Update the status of the generation if necessary
        self.update_generation()

        # 6. Write
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
        self.differences = None
        self.chi_squared = None
        self.chi_squared_table = None

    # -----------------------------------------------------------------

    @property
    def has_simulation(self):

        """
        This function ...
        :return:
        """

        return self.simulation is not None

    # -----------------------------------------------------------------

    @property
    def simulation_base_path(self):

        """
        This function ...
        :return:
        """

        return self.simulation.base_path

    # -----------------------------------------------------------------

    @property
    def generation_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.simulation_base_path)

    # -----------------------------------------------------------------

    @property
    def generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.generation_path)

    # -----------------------------------------------------------------

    @property
    def fitting_run_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.generations_path)

    # -----------------------------------------------------------------

    @property
    def generation_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.generation_path)

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.fitting_run_path)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SEDFitModelAnalyser, self).setup(**kwargs)

        # Get the simulation
        if not self.has_simulation: self.get_simulation(**kwargs)

        # Get the mock seds
        self.get_mock_seds(**kwargs)

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.fitting_run_name)

        # Load the weights table
        self.weights = WeightsTable.from_file(self.fitting_run.weights_table_path)

        # Initialize the differences table
        self.differences = FluxDifferencesTable()

    # -----------------------------------------------------------------

    def get_simulation(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        if "simulation" in kwargs: self.simulation = kwargs.pop("simulation")
        else: self.load_simulation()

    # -----------------------------------------------------------------

    @property
    def simulation_id(self):

        """
        This function ...
        :return:
        """

        # Determine simulation ID
        if self.config.name is not None:
            if self.config.id is not None: raise ValueError("Cannot specifiy both name and simulation ID")
            return get_simulation_id(self.config.remote, self.config.name)
        else: return self.config.id

    # -----------------------------------------------------------------

    def load_simulation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading simulation with ID '" + self.simulation_id + "' from remote host '" + self.config.remote + "' ...")

        # Load
        self.simulation = get_simulation_for_host(self.config.remote, self.simulation_id)

    # -----------------------------------------------------------------

    def get_mock_seds(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get a reference to the flux calculator
        simulationanalyser = self.get_analyser(**kwargs)

        # Get mock observed fluxes
        if simulationanalyser.basic_analyser.mock_seds_from_images is not None: self.mock_seds = simulationanalyser.basic_analyser.mock_seds_from_images
        elif simulationanalyser.basic_analyser.mock_seds is not None: self.mock_seds = simulationanalyser.basic_analyser.mock_seds
        else: raise ValueError("No mock observed fluxes could be obtained from the simulation analyser")

    # -----------------------------------------------------------------

    def get_analyser(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Load from kwargs
        if "simulation_analyser" in kwargs: return kwargs.pop("simulation_analyser")
        else: return self.load_analyser()

    # -----------------------------------------------------------------

    def load_analyser(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the simulation analyser ...")

        # Create simulation analyser
        analyser = SimulationAnalyser()
        analyser.config.extra = False

        # Run the analyser on the simulation
        analyser.run(simulation=self.simulation)

        # Return the analyser
        return analyser

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences between the observed and simulated SED ...")

        earth_instrument_name = "earth"
        # In the flux-density tables derived from the simulation (created by the ObservedFluxCalculator object),
        # search the one corresponding to the "earth" instrument
        #mock_sed_name = self.object_name + "_earth"
        #if mock_sed_name not in self.flux_calculator.mock_seds: raise RuntimeError("Could not find a mock observation SED for the 'earth' instrument")
        if earth_instrument_name not in self.mock_seds: raise RuntimeError("Could not find a mock observed SED for the 'earth' instrument")

        # Get the mock SED
        #mock_sed = self.flux_calculator.mock_seds[mock_sed_name]
        mock_sed = self.mock_seds[earth_instrument_name]

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

    @property
    def ndifferences(self):

        """
        This function ...
        :return:
        """

        return len(self.differences)

    # -----------------------------------------------------------------

    @property
    def nfree_parameters(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.nfree_parameters

    # -----------------------------------------------------------------

    @property
    def ndof(self):

        """
        This function ...
        :return:
        """

        return self.ndifferences - self.nfree_parameters - 1 # number of data points - number of fitted parameters - 1

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the chi squared value for this model ...")

        # The (reduced) chi squared value is the sum of all the terms (for each band),
        # divided by the number of degrees of freedom
        self.chi_squared = np.sum(self.differences["Chi squared term"]) / self.ndof

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

class ImageResidualsTable(SmartTable):

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
        super(ImageResidualsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Instrument", str, None, "Instrument")
        self.add_column_info("Band", str, None, "Band")
        self.add_column_info("Mean relative difference", float, None, "mean relative residual")
        self.add_column_info("Median relative difference", float, None, "median relative residual")
        self.add_column_info("Standard deviation of relative difference", float, None, "standard deviation of relative residual")
        self.add_column_info("Chi squared term", float, None, "Chi squared term")

    # -----------------------------------------------------------------

    def add_entry(self, instrument, band, mean, median, standard_deviation, chi_squared_term):

        """
        This function ...
        :param instrument:
        :param band:
        :param mean:
        :param median:
        :param standard_deviation:
        :param chi_squared_term
        :return:
        """

        self.add_row([instrument, band, mean, median, standard_deviation, chi_squared_term])

# -----------------------------------------------------------------

class ImagesFitModelAnalyser(FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImagesFitModelAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The fitting run
        self.fitting_run = None

        # The observed image maker
        self.image_maker = None

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

        # 2. Calculate the residual maps
        self.calculate_residuals()

        # 3. Calculate the chi squared for this model
        self.calculate_chi_squared()

        # 4. Load the chi squared table
        self.load_chi_squared_table()

        # 5. Update the status of the generation if necessary
        self.update_generation()

        # 6. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ImagesFitModelAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the residual maps of the observed and simulated images ...")

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def load_chi_squared_table(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def update_generation(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
