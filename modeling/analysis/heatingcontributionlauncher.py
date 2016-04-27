#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.launch Contains the BestModelLauncher class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import filesystem
from ...core.simulation.skifile import SkiFile
from ...core.launch.batchlauncher import BatchLauncher
from ...core.tools.logging import log
from ...core.launch.options import AnalysisOptions
from ...core.simulation.arguments import SkirtArguments

# -----------------------------------------------------------------

class HeatingContributionLauncher(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(HeatingContributionLauncher, self).__init__(config)

        # -- Attributes --

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The different contributing components
        self.contributions = ["tot", "old", "young", "ionizing"]

        # The path to the directory with the best model parameters
        self.best_path = None

        # The ski file corresponding to the best model
        self.ski = None

        # The ski files for the different contributions
        self.ski_files = {}

        # The paths to the analysis/heating/tot, analysis/heating/old, analysis/heating/young and analysis/heating/ionizing directory
        self.simulation_paths = {}
        self.output_paths = {}

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new BestModelLauncher instance
        launcher = cls()

        # Set the modeling path
        launcher.config.path = arguments.path

        # Return the new instance
        return launcher

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the ski file describing the best model
        self.load_ski()

        # 5. Adjust the ski files for the different contributors
        self.adjust_ski_files()

        # 6. Writing
        self.write()

        # 7. Launch the simulations
        self.launch()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(HeatingContributionLauncher, self).setup()

        # The path to the directory with the best model parameters
        self.best_path = filesystem.join(self.fit_path, "best")

        # Set the paths to the different simulation directories and corresponding output directories
        for contribution in self.contributions:

            # Set the simulation path
            simulation_path = filesystem.join(self.analysis_heating_path, contribution)

            # Create the directory if it is not present
            if not filesystem.is_directory(simulation_path): filesystem.create_directory(simulation_path)

            # Set the path to the output directory
            output_path = filesystem.join(simulation_path, "out")

            # Add the paths to the appropriate dictionaries
            self.simulation_paths[contribution] = simulation_path
            self.output_paths[contribution] = output_path

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file for the best fitting model ...")

        # Determine the path to the best model ski file
        path = filesystem.join(self.best_path, self.galaxy_name + ".ski")

        # Load the ski file
        self.ski = SkiFile(path)

    # -----------------------------------------------------------------

    def adjust_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Parameters of the wavelength grid
        min_wavelength = self.config.wavelengths.min * Unit(self.config.wavelengths.unit)
        max_wavelength = self.config.wavelengths.max * Unit(self.config.wavelengths.unit)
        points = self.config.wavelengths.points

        # Set the logarithmic wavelength grid
        self.ski.set_log_wavelength_grid(min_wavelength, max_wavelength, points, write=True)

        # Set the number of photon packages
        self.ski.setpackages(self.config.packages)

        # Enable all writing options for analysis
        self.ski.enable_all_writing_options()

        # Loop over the different contributions, create seperate ski file instance
        for contribution in self.contributions:

            # Create a copy of the ski file instance
            ski = copy.deepcopy(self.ski)

            if contribution == "tot": self.ski_files[contribution] = ski
            elif contribution == "old":

                # Remove all stellar components except for the old stellar bulge and disk
                ski.remove_stellar_components_except(["Old"])

                # Add the ski file instance to the dictionary
                self.ski_files[contribution] = ski

            elif contribution == "young":

                # Remove all stellar components except for the young stellar component
                ski.remove_stellar_components_except("")

                # Add the ski file instance to the dictionary
                self.ski_files[contribution] = ski

            elif contribution == "ionizing":

                # Remove all stellar components except for the ionizing stellar component
                ski.remove_stellar_components_except("")

                # Add the ski file instance to the dictionary
                self.ski_files[contribution] = ski

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file
        self.write_ski_files()

    # -----------------------------------------------------------------

    def write_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Loop over the contributions
        for contribution in self.ski_files:

            # Determine the path to the ski file
            path = filesystem.join(self.simulation_paths[contribution], self.galaxy_name + ".ski")

            # Save the ski file
            self.ski_files[contribution].save(path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Get the names of the filters for which we have photometry
        filter_names = self.get_filter_names()

        # Scheduling options
        scheduling_options = None

        # Analysis options
        analysis_options = AnalysisOptions()

        # Set options for extraction
        analysis_options.extraction.path = self.analysis_extr_path
        analysis_options.extraction.progress = True
        analysis_options.extraction.timeline = True

        # Set options for plotting
        analysis_options.plotting.path = self.analysis_plot_path
        analysis_options.plotting.progress = True
        analysis_options.plotting.timeline = True
        analysis_options.plotting.seds = True
        analysis_options.plotting.grids = True
        analysis_options.plotting.reference_sed = filesystem.join(self.phot_path, "fluxes.dat")

        # Set miscellaneous options
        analysis_options.misc.path = self.analysis_misc_path
        analysis_options.misc.rgb = True
        analysis_options.misc.wave = True
        analysis_options.misc.fluxes = True
        analysis_options.misc.images = True
        analysis_options.misc.observation_filters = filter_names

        # Create the SKIRT arguments object
        arguments = SkirtArguments()

        # Set the arguments
        arguments.ski_pattern = self.analysis_ski_path
        arguments.single = True
        arguments.input_path = self.analysis_in_path
        arguments.output_path = self.analysis_out_path
        arguments.logging.verbose = True

        # Run the simulation
        simulation = self.remote.run(arguments, scheduling_options=scheduling_options, analysis_options=analysis_options)

# -----------------------------------------------------------------
