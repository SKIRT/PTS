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
import numpy as np
from collections import defaultdict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem, time, tables
from ...core.simulation.arguments import SkirtArguments
from ...core.basics.filter import Filter
from ...core.simulation.skifile import SkiFile
from ...core.launch.batchlauncher import BatchLauncher
from ...core.tools.logging import log
from ...core.launch.parallelization import Parallelization

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
        self.parameters = []

        # The table with the parameter values for each simulation
        self.table = None

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
        fluxes_table_path = filesystem.join(self.phot_path, "fluxes.dat")
        #fluxes_table = tables.from_file(fluxes_table_path, format="ascii.ecsv")
        fluxes_table = tables.from_file(fluxes_table_path)
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
        self.launcher.config.analysis.plotting.reference_sed = filesystem.join(self.phot_path, "fluxes.dat") # the path to the reference SED (for plotting the simulated SED against the reference points)
        self.launcher.config.analysis.misc.fluxes = True  # Calculate observed fluxes
        self.launcher.config.analysis.misc.images = True  # Make observed images
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

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Set the paths to the directories to contain the launch scripts (job scripts) for the different remote hosts
        for host_id in self.launcher.host_ids:
            script_dir_path = filesystem.join(self.fit_scripts_path, host_id)
            if not filesystem.is_directory(script_dir_path): filesystem.create_directory(script_dir_path)
            self.launcher.set_script_path(host_id, script_dir_path)

        # Create a FUV filter object
        fuv = Filter.from_string("FUV")

        # Loop over the different parameter combinations
        for young_luminosity, ionizing_luminosity, dust_mass in self.parameters:

            # Create a unique name for this combination of parameter values
            simulation_name = time.unique_name()

            # Change the parameter values in the ski file
            self.ski.set_stellar_component_luminosity("Young stars", young_luminosity, fuv)
            self.ski.set_stellar_component_luminosity("Ionizing stars", ionizing_luminosity, fuv)
            self.ski.set_dust_component_mass(0, dust_mass)

            # Determine the directory for this simulation
            simulation_path = filesystem.join(self.fit_out_path, simulation_name)

            # Create the simulation directory
            filesystem.create_directory(simulation_path)

            # Create an 'out' directory within the simulation directory
            output_path = filesystem.join(simulation_path, "out")
            filesystem.create_directory(output_path)

            # Put the ski file with adjusted parameters into the simulation directory
            ski_path = filesystem.join(simulation_path, self.galaxy_name + ".ski")
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

        # Run the launcher, schedules the simulations
        simulations = self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in simulations:

            # Add the path to the modeling directory to the simulation object
            simulation.analysis.modeling_path = self.config.path

            # Save the simulation object
            simulation.save()

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
