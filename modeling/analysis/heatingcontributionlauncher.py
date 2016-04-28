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
from .dustheatinganalysiscomponent import DustHeatingAnalysisComponent
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ...core.launch.batchlauncher import BatchLauncher
from ...core.tools.logging import log
from ...core.simulation.arguments import SkirtArguments

# -----------------------------------------------------------------

class HeatingContributionLauncher(DustHeatingAnalysisComponent):
    
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

        # The path to the directory with the best model parameters
        self.best_path = None

        # The ski file corresponding to the best model
        self.ski = None

        # The ski files for the different contributions
        self.ski_files = {}

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
        self.best_path = fs.join(self.fit_path, "best")

        # Set options for the BatchLauncher
        self.launcher.config.shared_input = True  # The input directories for the different simulations are shared
        #self.launcher.config.group_simulations = True  # group multiple simulations into a single job
        self.launcher.config.remotes = self.config.remotes  # the remote hosts on which to run the simulations

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file for the best fitting model ...")

        # Determine the path to the best model ski file
        path = fs.join(self.best_path, self.galaxy_name + ".ski")

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
        points = self.config.wavelengths.npoints

        # Set the logarithmic wavelength grid
        self.ski.set_log_wavelength_grid(min_wavelength, max_wavelength, points, write=True)

        # Set the number of photon packages
        self.ski.setpackages(self.config.packages)

        # Enable all writing options for analysis
        self.ski.enable_all_writing_options()

        # Loop over the different contributions, create seperate ski file instance
        for contribution in self.contributions:

            # Debugging
            log.debug("Adjusting ski file for the contribution of the " + contribution + " stellar population ...")

            # Create a copy of the ski file instance
            ski = copy.deepcopy(self.ski)

            # Remove other stellar components, except for the contribution of the total stellar population
            if contribution != "total": ski.remove_stellar_components_except(self.component_names[contribution])

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

        # Copy the input maps (if necessary)
        self.copy_maps()

        # Write the ski file
        self.write_ski_files()

    # -----------------------------------------------------------------

    def copy_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Copying the input maps ...")

        # Determine the paths to the input maps in the fit/in directory
        fit_in_path = fs.join(self.fit_path, "in")
        old_path = fs.join(fit_in_path, "old_stars.fits")
        young_path = fs.join(fit_in_path, "young_stars.fits")
        ionizing_path = fs.join(fit_in_path, "ionizing_stars.fits")
        dust_path = fs.join(fit_in_path, "dust.fits")

        # Copy the files to the analysis/in directory (if necessary)
        if not fs.has_file(self.analysis_in_path, fs.name(old_path)): fs.copy_file(old_path, self.analysis_in_path)
        if not fs.has_file(self.analysis_in_path, fs.name(young_path)): fs.copy_file(young_path, self.analysis_in_path)
        if not fs.has_file(self.analysis_in_path, fs.name(ionizing_path)): fs.copy_file(ionizing_path, self.analysis_in_path)
        if not fs.has_file(self.analysis_in_path, fs.name(dust_path)): fs.copy_file(dust_path, self.analysis_in_path)

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
            path = self.ski_paths[contribution]

            # Debugging
            log.debug("Writing the ski file for the " + contribution + " stellar population to '" + path + "' ...")

            # Save the ski file
            self.ski_files[contribution].saveto(path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Determine the path to the analysis/heating/scripts path (store batch scripts there for manual inspection)
        scripts_path = fs.join(self.analysis_heating_path, "scripts")
        if not fs.is_directory(scripts_path): fs.create_directory(scripts_path)
        for host_id in self.launcher.host_ids:
            script_dir_path = fs.join(scripts_path, host_id)
            if not fs.is_directory(script_dir_path): fs.create_directory(script_dir_path)
            self.launcher.set_script_path(host_id, script_dir_path)

        # Loop over the contributions
        for contribution in self.ski_paths:

            # Determine a name for this simulation
            simulation_name = self.galaxy_name + "_heating_" + contribution

            # Get the ski path for this simulation
            ski_path = self.ski_paths[contribution]

            # Get the local output path for the simulation
            output_path = self.output_paths[contribution]

            # Create the SKIRT arguments object
            arguments = create_arguments(ski_path, self.analysis_in_path, output_path)

            # Debugging
            log.debug("Adding the simulation of the contribution of the " + contribution + " stellar population to the queue ...")

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(arguments, simulation_name)

            # Set scheduling options (for the different remote hosts with a scheduling system)
            #for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, simulation_name, self.scheduling_options[host_id])

        # Run the launcher, schedules the simulations
        simulations = self.launcher.run()

        # Loop over the scheduled simulations (if something has to be set)
        #for simulation in simulations: pass

# -----------------------------------------------------------------

def create_arguments(ski_path, input_path, output_path):

    """
    This function ...
    :param ski_path:
    :param input_path:
    :param output_path:
    :return:
    """

    # Create a new SkirtArguments object
    arguments = SkirtArguments()

    # The ski file pattern
    arguments.ski_pattern = ski_path
    arguments.recursive = False
    arguments.relative = False

    # Input and output
    arguments.input_path = input_path
    arguments.output_path = output_path

    # Parallelization settings
    arguments.parallel.threads = None
    arguments.parallel.processes = None

    # Return the SKIRT arguments object
    return arguments

# -----------------------------------------------------------------