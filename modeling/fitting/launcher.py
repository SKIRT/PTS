#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.launcher Contains the abstract FittingModelLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod
from collections import defaultdict

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools import time, tables
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.launch.options import LoggingOptions
from ...core.basics.filter import Filter
from ...core.simulation.skifile import SkiFile
from ...core.launch.batchlauncher import BatchLauncher
from ...core.tools.logging import log
from ...core.launch.options import AnalysisOptions
from ...magic.misc.kernels import AnianoKernels

# -----------------------------------------------------------------

class FittingModelLauncher(FittingComponent):

    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingModelLauncher, self).__init__(config)

        # -- Attributes --

        # The parameter ranges
        self.ranges = dict()

        # The names of the filters for which we have photometry
        self.filter_names = []

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The ski file
        self.ski = None

        # The parameter combinations
        self.parameters = defaultdict(list)

        # The table with the parameter values for each simulation
        self.table = None

        # A dictionary with the scheduling options for the different remote hosts
        self.scheduling_options = dict()

        # The animations
        self.scatter_animation = None
        self.fuv_young_animation = None
        self.fuv_ionizing_animation = None
        self.dust_mass_animation = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary input
        self.load_input()

        # 3. Set the parameter ranges
        self.set_parameter_ranges()

        # 4. Initialize the animations, if requested
        if self.config.visualise: self.initialize_animations()

        # 5. Set the combinations of parameter values
        self.set_parameters()

        # 6. Set the parallelization schemes for the different remote hosts
        self.set_parallelization()

        # 7. Estimate the runtimes for the different remote hosts
        self.estimate_runtimes()

        # 8. Launch the simulations for different parameter values
        self.launch()

        # 9. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingModelLauncher, self).setup()

        # Get the names of the filters for which we have photometry
        self.filter_names = self.get_observed_filter_names()

        # Set options for the BatchLauncher: basic options
        self.launcher.config.shared_input = True  # The input directories for the different simulations are shared
        self.launcher.config.group_simulations = True  # group multiple simulations into a single job (because a very large number of simulations will be scheduled)
        self.launcher.config.remotes = self.config.remotes  # the remote hosts on which to run the simulations
        self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file

        # Set options for the BatchLauncher: simulation analysis options
        self.launcher.config.analysis.extraction.path = self.fit_res_path
        self.launcher.config.analysis.misc.path = self.fit_res_path # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        self.launcher.config.analysis.plotting.path = self.fit_plot_path # The base directory where all of the simulations will have a seperate directory with the plotting analysis output
        self.launcher.config.analysis.extraction.timeline = True # extract the simulation timeline
        self.launcher.config.analysis.plotting.seds = True  # Plot the output SEDs
        self.launcher.config.analysis.plotting.reference_sed = self.observed_sed_path # the path to the reference SED (for plotting the simulated SED against the reference points)
        #self.launcher.config.analysis.plotting.reference_sed = self.observed_sed_dustpedia_path # the path to the DustPedia SED
        #self.launcher.config.analysis.misc.fluxes = True  # Calculate observed fluxes
        #self.launcher.config.analysis.misc.images = True  # Make observed images
        self.launcher.config.analysis.misc.observation_filters = self.filter_names  # The filters for which to create the observations
        self.launcher.config.analysis.plotting.format = "png" # plot in PNG format so that an animation can be made from the fit SEDs

        # Set remote for the 'extra' simulations
        self.launcher.config.extra_remote = "nancy"

    # -----------------------------------------------------------------

    @abstractmethod
    def load_input(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    def set_parameter_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # The given ranges are relative to the best or initial value
        if self.config.relative:

            # TODO: get parameter values of best fitting model

            # Get the current values in the ski file prepared by InputInitializer
            young_luminosity_guess, young_filter = self.ski.get_stellar_component_luminosity("Young stars")
            ionizing_luminosity_guess, ionizing_filter = self.ski.get_stellar_component_luminosity("Ionizing stars")
            dust_mass_guess = self.ski.get_dust_component_mass(0)

            # Set the range of the FUV luminosity of the young stellar population
            min_value = self.config.young[0] * young_luminosity_guess
            max_value = self.config.young[1] * young_luminosity_guess
            self.ranges["FUV young"] = (min_value, max_value)

            # Set the range of the FUV luminosity of the ionizing stellar population
            min_value = self.config.ionizing[0] * ionizing_luminosity_guess
            max_value = self.config.ionizing[1] * ionizing_luminosity_guess
            self.ranges["FUV ionizing"] = (min_value, max_value)

            # Set the range of the dust mass
            min_value = self.config.dust[0] * dust_mass_guess
            max_value = self.config.dust[1] * dust_mass_guess
            self.ranges["Dust mass"] = (min_value, max_value)

        else:

            # Set the ranges directly from the command-line arguments
            self.ranges["FUV young"] = self.config.young
            self.ranges["FUV ionizing"] = self.config.ionizing
            self.ranges["Dust mass"] = self.config.dust

    # -----------------------------------------------------------------

    @abstractmethod
    def initialize_animations(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    @abstractmethod
    def set_parameters(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    @abstractmethod
    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    @abstractmethod
    def estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    def load_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the parameter table ...")

        # Load the parameter table
        self.table = tables.from_file(self.parameter_table_path, format="ascii.ecsv", fix_string_length=("Simulation name", 24))

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file ...")

        # Open the ski file (created by InputInitializer)
        self.ski = SkiFile(self.fit_ski_path)

    # -----------------------------------------------------------------

    @property
    def number_of_models(self):

        """
        This function ...
        :return:
        """

        return len(self.parameters["FUV young"])

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Set the paths to the directories to contain the launch scripts (job scripts) for the different remote hosts
        for host_id in self.launcher.host_ids:
            script_dir_path = fs.join(self.fit_scripts_path, host_id)
            if not fs.is_directory(script_dir_path): fs.create_directory(script_dir_path)
            self.launcher.set_script_path(host_id, script_dir_path)

        # Set the paths to the screen output directories (for debugging) for remotes without a scheduling system for jobs
        for host_id in self.launcher.no_scheduler_host_ids: self.launcher.enable_screen_output(host_id)

        # Create a FUV filter object
        fuv = Filter.from_string("FUV")

        # Loop over the different parameter combinations
        for i in range(self.number_of_models):

            # Get the parameter values
            young_luminosity = self.parameters["FUV young"][i]
            ionizing_luminosity = self.parameters["FUV ionizing"][i]
            dust_mass = self.parameters["Dust mass"][i]

            # Create a unique name for this combination of parameter values
            simulation_name = time.unique_name()

            # Change the parameter values in the ski file
            self.ski.set_stellar_component_luminosity("Young stars", young_luminosity, fuv.centerwavelength() * Unit("micron"))
            self.ski.set_stellar_component_luminosity("Ionizing stars", ionizing_luminosity, fuv.centerwavelength() * Unit("micron"))
            self.ski.set_dust_component_mass(0, dust_mass)

            # Determine the directory for this simulation
            simulation_path = fs.join(self.fit_out_path, simulation_name)

            # Create the simulation directory
            fs.create_directory(simulation_path)

            # Create an 'out' directory within the simulation directory
            output_path = fs.join(simulation_path, "out")
            fs.create_directory(output_path)

            # Put the ski file with adjusted parameters into the simulation directory
            ski_path = fs.join(simulation_path, self.galaxy_name + ".ski")
            self.ski.saveto(ski_path)

            # Create the SKIRT simulation definition
            definition = SingleSimulationDefinition(ski_path, self.fit_in_path, output_path)

            # Debugging
            log.debug("Adding a simulation to the queue with:")
            log.debug(" - ski path: " + definition.ski_path)
            log.debug(" - output path: " + definition.output_path)

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name)

            # Set scheduling options (for the different remote hosts with a scheduling system)
            for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, simulation_name, self.scheduling_options[host_id])

            # Add an entry to the parameter table
            self.table.add_row([simulation_name, young_luminosity, ionizing_luminosity, dust_mass])

        # Add simulations to the extra queue to calculate the contribution of the various stellar components
        self.add_contribution_simulations()

        # Add simulation to the extra queue to create simulated images
        self.add_images_simulation()

        # Run the launcher, schedules the simulations
        simulations = self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in simulations:

            # Skip simulations for which analysis options are not defined (the contribution simulations and the images simulation)
            if simulation.analysis is None or simulation.name == "images": continue

            # Add the path to the modeling directory to the simulation object
            simulation.analysis.modeling_path = self.config.path

            # Save the simulation object
            simulation.save()

    # -----------------------------------------------------------------

    def add_contribution_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding simulations to check the contribution of the various stellar components ...")

        # Loop over the contributions
        contributions = ["old", "young", "ionizing"]
        for contribution in contributions:

            # Simulation name
            simulation_name = contribution

            # Set the ski file path
            ski_path = fs.join(self.fit_best_path, contribution, self.galaxy_name + ".ski")

            # Create the SKIRT simulation definition
            definition = SingleSimulationDefinition(ski_path, self.fit_in_path, fs.join(self.fit_best_path, contribution))

            # Create the AnalysisOptions instance
            #analysis_options = AnalysisOptions()

            # Add the arguments object
            self.launcher.add_to_extra_queue(definition, simulation_name, share_input=True)

            # Set scheduling options if necessary
            #for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, simulation_name, self.scheduling_options[host_id])

    # -----------------------------------------------------------------

    def add_images_simulation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding simulation to generate simulated images ...")

        # Set the ski file path
        ski_path = fs.join(self.fit_best_path, "images", self.galaxy_name + ".ski")

        # Set the input and output path
        input_path = self.fit_in_path
        output_path = fs.join(self.fit_best_path, "images")

        # Create the SKIRT simulation definition
        definition = SingleSimulationDefinition(ski_path, input_path, output_path)

        # Create the logging options instance
        logging = LoggingOptions(verbose=True, memory=True)

        # Create the analysis options object
        analysis_options = self.create_analysis_options_images_simulation(output_path)

        # Add the arguments object
        self.launcher.add_to_extra_queue(definition, "images", logging, analysis_options, share_input=True)

        # Set scheduling options if necessary
        #for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, "images", self.scheduling_options[host_id])

    # -----------------------------------------------------------------

    def create_analysis_options_images_simulation(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Create the AnalysisOptions instance
        analysis_options = AnalysisOptions()

        # Set extraction options
        analysis_options.extraction.path = output_path
        analysis_options.extraction.progress = True
        analysis_options.extraction.timeline = True
        analysis_options.extraction.memory = True

        # Set plotting options
        analysis_options.plotting.path = output_path
        analysis_options.plotting.progress = True
        analysis_options.plotting.timeline = True
        analysis_options.plotting.memory = True
        analysis_options.plotting.seds = True
        analysis_options.plotting.reference_sed = self.observed_sed_path

        # Set the paths to the for each image (except for the SPIRE images)
        kernel_paths = dict()
        aniano = AnianoKernels()
        pacs_red_psf_path = aniano.get_psf_path("PACS_160")
        for filter_name in self.filter_names:
            if "SPIRE" in filter_name: continue
            kernel_paths[filter_name] = pacs_red_psf_path

        # Set misc options
        analysis_options.misc.path = output_path
        analysis_options.misc.images = True
        analysis_options.misc.observation_filters = self.filter_names
        analysis_options.misc.make_images_remote = "nancy"
        analysis_options.misc.images_wcs = self.reference_path
        analysis_options.misc.images_unit = "MJy/sr"
        analysis_options.misc.images_kernels = kernel_paths

        # Return the analysis options
        return analysis_options

    # -----------------------------------------------------------------

    @abstractmethod
    def write(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    def write_parameter_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameter table ...")

        # Set the units of the parameter table
        self.table["FUV young"].unit = "Lsun_FUV"
        self.table["FUV ionizing"].unit = "Lsun_FUV"
        self.table["Dust mass"].unit = "Msun"

        # Write the parameter table
        tables.write(self.table, self.parameter_table_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def write_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the animations ...")

        # Save the animation of the parameter values as points in a 3D plot
        path = fs.join(self.visualisation_path, time.unique_name("fittingmodellauncher") + ".gif")
        self.scatter_animation.save(path)

        # Save the animation of the distribution of values for the FUV luminosity of the young stars
        path = fs.join(self.visualisation_path, time.unique_name("fittingmodellauncher_fuvyoung") + ".gif")
        if self.fuv_young_animation is not None: self.fuv_young_animation.save(path)

        # Save the animation of the distribution of values for the FUV luminosity of the ionizing stars
        path = fs.join(self.visualisation_path, time.unique_name("fittingmodellauncher_fuvionizing") + ".gif")
        if self.fuv_ionizing_animation is not None: self.fuv_ionizing_animation.save(path)

        # Save the animation of the distribution of values for the dust mass
        path = fs.join(self.visualisation_path, time.unique_name("fittingmodellauncher_dustmass") + ".gif")
        if self.dust_mass_animation is not None: self.dust_mass_animation.save(path)

# -----------------------------------------------------------------
