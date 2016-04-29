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
from ...core.tools import filesystem, tables
from ..core.sed import ObservedSED

# -----------------------------------------------------------------

class FitModelAnalyser(FittingComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FitModelAnalyser, self).__init__(config)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The timeline extractor
        self.te = None

        # The log file of the simulation
        self.log_file = None

        # The ski file corresponding to the simulation
        self.ski = None

        # The flux calculator
        self.flux_calculator = None

        # The weights given to each band for the calculation of the chi squared
        self.weights = None

        # The observed fluxes
        self.fluxes = None

        # The flux differences table
        self.differences = None

        # The calculated chi squared value
        self.chi_squared = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline_extractor, flux_calculator):

        """
        This function ...
        :param simulation:
        :param timeline_extractor:
        :param flux_calculator:
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation, timeline_extractor, flux_calculator)

        # 2. Load the log file of the simulation
        self.load_log_file()

        # 3. Load the observed SED
        self.load_observed_sed()

        # 4. Calculate the differences
        self.calculate_differences()

        # 5. Calculate the chi squared for this model
        self.calculate_chi_squared()

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
        self.te = None
        self.log_file = None
        self.ski = None
        self.flux_calculator = None
        self.fluxes = None
        self.differences = None
        self.chi_squared = None

    # -----------------------------------------------------------------

    def setup(self, simulation, timeline_extractor, flux_calculator):

        """
        This function ...
        :param simulation:
        :param timeline_extractor:
        :param flux_calculator:
        :return:
        """

        # Call the setup function of the base class
        super(FitModelAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

        # Make a reference to the timeline extractor
        self.te = timeline_extractor

        # Make a local reference to the flux calculator
        if flux_calculator is None:
            raise RuntimeError("No ObservedFluxCalculator found; the calculate_observed_fluxes flag must be enabled on "
                               "each simulation that is part of the radiative transfer modeling")
        self.flux_calculator = flux_calculator

        # Load the weights table
        #self.weights = tables.from_file(self.weights_table_path, fix_floats=True) # For some reason, the weights are parsed as strings instead of floats (but not from the command line!!??)
        self.weights = tables.from_file(self.weights_table_path, format="ascii.ecsv")

        # Initialize the differences table
        names = ["Instrument", "Band", "Flux difference", "Relative difference", "Chi squared term"]
        data = [[], [], [], [], []]
        dtypes = ["S5", "S7", "float64", "float64", "float64"]
        self.differences = tables.new(data, names, dtypes=dtypes)

        # Load the ski file
        self.ski = self.simulation.parameters()

    # -----------------------------------------------------------------

    def load_log_file(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the log file produced by the simulation ...")

        # Get the log file produced by the simulation
        self.log_file = self.simulation.log_file

    # -----------------------------------------------------------------

    def load_observed_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed SED ...")

        # Determine the path to the fluxes table
        fluxes_path = filesystem.join(self.phot_path, "fluxes.dat")

        # Load the observed SED
        self.fluxes = ObservedSED.from_file(fluxes_path)

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
        table_name = self.galaxy_name + "_earth"
        if table_name not in self.flux_calculator.tables: raise RuntimeError("Could not find a flux-density table for the 'earth' instrument")

        # Get the table
        table = self.flux_calculator.tables[table_name]

        # Loop over the entries in the fluxdensity table (SED) derived from the simulation
        for i in range(len(table)):

            #observatory = table["Observatory"][i]
            instrument = table["Instrument"][i]
            band = table["Band"][i]
            wavelength = table["Wavelength"][i]
            fluxdensity = table["Flux"][i]

            # Find the corresponding flux in the SED derived from observation
            observed_fluxdensity = self.fluxes.flux_for_band(instrument, band, unit="Jy").value

            # Find the corresponding flux error in the SED derived from observation
            observed_fluxdensity_error = self.fluxes.error_for_band(instrument, band, unit="Jy").average.to("Jy").value

            # If no match with (instrument, band) is found in the observed SED
            if observed_fluxdensity is None:
                log.warning("The observed flux density could not be found for the " + instrument + " " + band + " band")
                continue

            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # Find the index of the current band in the weights table
            index = tables.find_index(self.weights, key=[instrument, band], column_name=["Instrument", "Band"])

            # Get the weight
            weight = self.weights["Weight"][index] # apparently, this is a string, so parsing the table went wrong ...
            weight = float(weight)

            # Calculate the chi squared term
            chi_squared_term = weight * difference ** 2 / observed_fluxdensity_error ** 2

            # Add an entry to the differences table
            self.differences.add_row([instrument, band, difference, relative_difference, chi_squared_term])

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the chi squared value for this model ...")

        # Calculate the degrees of freedom
        dof = len(self.fluxes.table) - 3. - 1.  # number of data points - number of fitted parameters - 1

        # The (reduced) chi squared value is the sum of all the terms (for each band),
        # divided by the number of degrees of freedom
        self.chi_squared = np.sum(self.differences["Chi squared term"]) / dof

        # Debugging
        log.debug("Found a chi squared value of " + str(self.chi_squared))

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the runtime
        self.write_runtime()

        # Write the flux differences
        self.write_differences()

        # Write the chi-squared value
        self.write_chi_squared()

    # -----------------------------------------------------------------

    def write_runtime(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the runtime of this simulation ...")

        # Get the name of the host on which the simulation was run
        #host = self.log_file.host
        host_id = self.simulation.host_id
        cluster_name = self.simulation.cluster_name
        if cluster_name is None: cluster_name = "--"

        # Get the parallelization object from the simulation
        parallelization = self.simulation.parallelization

        # Get the paralleliation properties
        cores = parallelization.cores
        hyperthreads = parallelization.threads_per_core
        processes = parallelization.processes

        # Get the runtime in seconds
        runtime = self.log_file.total_runtime

        # Get the number of photon packages
        packages = self.ski.packages()

        # Open the timing file in 'append' mode
        timing_file = open(self.timing_table_path, 'a')

        # Initialize a list to contain the values of the row
        row = []

        # Columns:
        # "Submission time"
        # "Host id"
        # "Cluster name"
        # "Cores"
        # "Hyperthreads per core"
        # "Processes"
        # "Packages"
        # "Total runtime"
        # "Serial runtime"
        # "Parallel runtime"
        # "
        row.append(self.simulation.name)
        row.append(host_id)
        row.append(cluster_name)
        row.append(str(cores))
        row.append(str(hyperthreads))
        row.append(str(processes))
        row.append(str(packages))
        row.append(str(runtime))
        row.append(str())

        # Add the row to the runtime file
        timing_file.write(" ".join(row) + "\n")

        # Close the file
        timing_file.close()

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with the flux-density differences for the current model ...")

        # Determine the path to the differences table
        path = filesystem.join(self.fit_res_path, self.simulation.name, "differences.dat")

        # Save the differences table
        tables.write(self.differences, path)

    # -----------------------------------------------------------------

    def write_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding the chi squared value for the current model to the chi squared data file ...")

        # Open the chi squared file (in 'append' mode)
        resultfile = open(self.chi_squared_table_path, 'a')

        # Add a line to the chi squared file containing the simulation name and the chi squared value
        resultfile.write(self.simulation.name + " " + str(self.chi_squared) + "\n")

        # Close the output file
        resultfile.close()

# -----------------------------------------------------------------
