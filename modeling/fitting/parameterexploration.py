#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.parameterexploration Contains the ParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools import time, tables
from ...core.simulation.arguments import SkirtArguments
from ...core.basics.filter import Filter
from ...core.simulation.skifile import SkiFile
from ...core.launch.batchlauncher import BatchLauncher
from ...core.tools.logging import log
from ...core.launch.parallelization import Parallelization
from ...core.launch.options import AnalysisOptions

# -----------------------------------------------------------------

class ParameterExplorer(FittingComponent):
    
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
        super(ParameterExplorer, self).__init__(config)

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

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ParameterExplorer, self).setup()

        # Get the names of the filters for which we have photometry
        filter_names = []
        fluxes_table_path = fs.join(self.phot_path, "fluxes.dat")
        fluxes_table = tables.from_file(fluxes_table_path, format="ascii.ecsv")
        # Loop over the entries in the fluxes table, get the filter
        for entry in fluxes_table:
            # Get the filter
            filter_id = entry["Instrument"] + "." + entry["Band"]
            filter_names.append(filter_id)

        # Set options for the BatchLauncher: basic options
        self.launcher.config.shared_input = True  # The input directories for the different simulations are shared
        self.launcher.config.group_simulations = True  # group multiple simulations into a single job (because a very large number of simulations will be scheduled)
        self.launcher.config.remotes = self.config.remotes  # the remote hosts on which to run the simulations

        # Set options for the BatchLauncher: simulation analysis options
        self.launcher.config.analysis.extraction.path = self.fit_res_path
        self.launcher.config.analysis.misc.path = self.fit_res_path # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        self.launcher.config.analysis.plotting.path = self.fit_plot_path # The base directory where all of the simulations will have a seperate directory with the plotting analysis output
        self.launcher.config.analysis.extraction.timeline = True # extract the simulation timeline
        self.launcher.config.analysis.plotting.seds = True  # Plot the output SEDs
        self.launcher.config.analysis.plotting.reference_sed = fs.join(self.phot_path, "fluxes.dat") # the path to the reference SED (for plotting the simulated SED against the reference points)
        self.launcher.config.analysis.misc.fluxes = True  # Calculate observed fluxes
        #self.launcher.config.analysis.misc.images = True  # Make observed images
        self.launcher.config.analysis.misc.observation_filters = filter_names  # The filters for which to create the observations
        self.launcher.config.analysis.plotting.format = "png" # plot in PNG format so that an animation can be made from the fit SEDs
        self.launcher.config.analysis.timing_table_path = self.timing_table_path # The path to the timing table file
        self.launcher.config.analysis.memory_table_path = self.memory_table_path # The path to the memory table file

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

    def set_parallelization(self):

        """
        This function sets the parallelization scheme for those remote hosts used by the batch launcher that use
        a scheduling system (the parallelization for the other hosts is left up to the batch launcher and will be
        based on the current load of the correponding system).
        :return:
        """

        # Loop over the IDs of the hosts used by the batch launcher that use a scheduling system
        for host in self.launcher.scheduler_hosts:

            # Get the number of cores per node for this host
            cores_per_node = host.clusters[host.cluster_name].cores

            # Determine the number of cores corresponding to 4 full nodes
            cores = cores_per_node * 4

            # Use 1 core for each process (assume there is enough memory)
            processes = cores

            # Determine the number of threads per core
            if host.use_hyperthreading: threads_per_core = host.clusters[host.cluster_name].threads_per_core
            else: threads_per_core = 1

            # Create a Parallelization instance
            parallelization = Parallelization(cores, threads_per_core, processes)

            # Set the parallelization for this host
            self.launcher.set_parallelization_for_host(host.id, parallelization)

    # -----------------------------------------------------------------

    @property
    def number_of_models(self):

        """
        This function ...
        :return:
        """

        return len(self.parameters["FUV young"])

    # -----------------------------------------------------------------

    def simulate(self):

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
            self.ski.set_stellar_component_luminosity("Young stars", young_luminosity, fuv)
            self.ski.set_stellar_component_luminosity("Ionizing stars", ionizing_luminosity, fuv)
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

            # Create the SKIRT arguments object
            arguments = create_arguments(ski_path, self.fit_in_path, output_path)

            # Debugging
            log.debug("Adding a simulation to the queue with:")
            log.debug(" - ski path: " + arguments.ski_pattern)
            log.debug(" - output path: " + arguments.output_path)

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(arguments, simulation_name)

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

            # Create the SkirtArguments instance
            arguments = SkirtArguments()

            # Set the ski file path
            ski_path = fs.join(self.fit_best_path, contribution, self.galaxy_name + ".ski")
            arguments.ski_pattern = ski_path
            arguments.single = True
            arguments.recursive = False
            arguments.relative = False
            arguments.parallel.threads = None
            arguments.parallel.processes = None
            arguments.input_path = self.fit_in_path
            arguments.output_path = fs.join(self.fit_best_path, contribution)

            # Create the AnalysisOptions instance
            analysis_options = AnalysisOptions()


            # Add the arguments object
            self.launcher.add_to_extra_queue(arguments, analysis_options, simulation_name)

            # Set scheduling options if necessary
            for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, simulation_name, self.scheduling_options[host_id])

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

        # Create the SkirtArguments instance
        arguments = SkirtArguments.single(ski_path, input_path, output_path, verbose=True, memory=True)

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
        analysis_options.plotting.reference_sed = fs.join(self.phot_path, "fluxes.dat")

        # Set misc options
        analysis_options.misc.path = output_path
        analysis_options.misc.images = True
        analysis_options.misc.observation_filters =
        analysis_options.misc.make_images_remote =
        analysis_options.misc.images_wcs =
        analysis_options.misc.images_unit =
        analysis_options.misc.images_kernels =

        # Add the arguments object
        self.launcher.add_to_extra_queue(arguments, analysis_options, "images")

        # Set scheduling options if necessary
        for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, "images", self.scheduling_options[host_id])

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Create and write a table with the parameter values for each simulation
        self.write_parameter_table()

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
