#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.remotelauncher Contains the SkirtRemoteLauncher class, which can be used to launch
#  SKIRT simulations on a remote machine in a convenient way.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import the relevant PTS classes and modules
from .analyser import SimulationAnalyser
from ..simulation.remote import SkirtRemote
from ..simulation.arguments import SkirtArguments
from ..basics.configurable import Configurable
from ..test.resources import ResourceEstimator
from ..tools import filesystem
from ..tools.logging import log
from .options import SchedulingOptions

# -----------------------------------------------------------------

class SkirtRemoteLauncher(Configurable):

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
        super(SkirtRemoteLauncher, self).__init__(config, "core")

        # -- Attributes --

        # Create the SKIRT remote execution context
        self.remote = SkirtRemote()

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # Initialize a list to contain the retrieved finished simulations
        self.simulations = []

        # Set the paths to None initialy
        self.base_path = None
        self.input_path = None
        self.output_path = None
        self.extr_path = None
        self.plot_path = None
        self.misc_path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create and a new SkirtLauncher instance
        launcher = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Remote host and cluster (if applicable)
        launcher.config.remote = arguments.remote
        launcher.config.cluster = arguments.cluster

        # Walltime
        launcher.config.walltime = arguments.walltime

        # Ski file
        launcher.config.arguments.ski_pattern = arguments.filepath
        launcher.config.arguments.recursive = False
        launcher.config.arguments.relative = False

        # Simulation logging
        launcher.config.arguments.logging.brief = arguments.brief
        launcher.config.arguments.logging.verbose = arguments.verbose
        launcher.config.arguments.logging.memory = arguments.memory
        launcher.config.arguments.logging.allocation = arguments.allocation

        # Parallelization
        if arguments.parallel is not None:
            launcher.config.arguments.parallel.processes = arguments.parallel[0]
            launcher.config.arguments.parallel.threads = arguments.parallel[1]
        else:
            launcher.config.arguments.parallel.processes = None
            launcher.config.arguments.parallel.threads = None

        # Other simulation arguments
        launcher.config.arguments.emulate = arguments.emulate
        launcher.config.arguments.single = True  # For now, we only allow single simulations

        # Extraction
        launcher.config.analysis.extraction.progress = arguments.extractprogress
        launcher.config.analysis.extraction.timeline = arguments.extracttimeline
        launcher.config.analysis.extraction.memory = arguments.extractmemory

        # Plotting
        launcher.config.analysis.plotting.seds = arguments.plotseds
        launcher.config.analysis.plotting.grids = arguments.plotgrids
        launcher.config.analysis.plotting.progress = arguments.plotprogress
        launcher.config.analysis.plotting.timeline = arguments.plottimeline
        launcher.config.analysis.plotting.memory = arguments.plotmemory
        launcher.config.analysis.plotting.reference_sed = arguments.refsed

        # Miscellaneous
        launcher.config.analysis.misc.rgb = arguments.makergb
        launcher.config.analysis.misc.wave = arguments.makewave
        launcher.config.analysis.misc.fluxes = arguments.fluxes
        launcher.config.analysis.misc.images = arguments.images
        launcher.config.analysis.misc.observation_filters = arguments.filters
        launcher.config.analysis.misc.observation_instruments = arguments.instruments
        launcher.config.analysis.misc.images_wcs = arguments.wcs
        launcher.config.analysis.misc.images_unit = arguments.unit

        # Keep remote input and output
        launcher.config.keep = arguments.keep

        # Retrieve types
        launcher.config.retrieve_types = arguments.retrieve

        # Return the new launcher
        return launcher

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Set the parallelization scheme
        if not self.has_parallelization: self.set_parallelization()
        else: self.check_parallelization()

        # 3. Run the simulation
        self.simulate()

        # 4. Retrieve the simulations that are finished
        self.retrieve()

        # 5. Analyse the output of the retrieved simulations
        self.analyse()

    # -----------------------------------------------------------------

    @property
    def has_parallelization(self):

        """
        This function ...
        :return:
        """

        # Check whether the number of processes and the number of threads are both defined
        return self.config.arguments.parallel.processes is not None and self.config.arguments.parallel.threads is not None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtRemoteLauncher, self).setup()

        # Setup the remote execution context
        self.remote.setup(self.config.remote, self.config.cluster)

        # Set the paths
        self.base_path = filesystem.directory_of(self.config.arguments.ski_pattern) if "/" in self.config.arguments.ski_pattern else filesystem.cwd()
        self.input_path = filesystem.join(self.base_path, "in")
        self.output_path = filesystem.join(self.base_path, "out")
        self.extr_path = filesystem.join(self.base_path, "extr")
        self.plot_path = filesystem.join(self.base_path, "plot")
        self.misc_path = filesystem.join(self.base_path, "misc")

        # Check if an input directory exists
        if not filesystem.is_directory(self.input_path): self.input_path = None

        # Set the paths for the simulation
        self.config.arguments.input_path = self.input_path
        self.config.arguments.output_path = self.output_path

        # Create the output directory if necessary
        if not filesystem.is_directory(self.output_path): filesystem.create_directory(self.output_path, recursive=True)

        # Create the extraction directory if necessary
        if self.config.extraction.progress or self.config.extraction.timeline or self.config.extraction.memory:
            if not filesystem.is_directory(self.extr_path): filesystem.create_directory(self.extr_path, recursive=True)

        # Create the plotting directory if necessary
        if self.config.plotting.seds or self.config.plotting.grids or self.config.plotting.progress \
            or self.config.plotting.timeline or self.config.plotting.memory:
            if not filesystem.is_directory(self.plot_path): filesystem.create_directory(self.plot_path, recursive=True)

        # Create the 'misc' directory if necessary
        if self.config.misc.fluxes or self.config.misc.images:
            if not filesystem.is_directory(self.misc_path): filesystem.create_directory(self.misc_path, recursive=True)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        #log.info("free cores: " + str(self.remote.free_cores))
        #log.info("free memory: " + str(self.remote.free_memory))
        #log.info("free space: " + str(self.remote.free_space))
        #log.info("cores: " + str(self.remote.cores))
        #log.info("cpu load: " + str(self.remote.cpu_load))
        #log.info("memory load: " + str(self.remote.memory_load))

        # Inform the user
        log.info("Determining the parallelization scheme by estimating the memory requirements...")

        # Calculate the amount of required memory for this simulation
        estimator = ResourceEstimator()
        estimator.run(self.config.arguments.ski_pattern)

        # Calculate the maximum number of processes based on the memory requirements
        processes = int(self.remote.free_memory / estimator.memory)

        # If there is too little free memory for the simulation, the number of processes will be smaller than one
        if processes < 1:

            # Exit with an error
            log.error("Not enough memory available to run this simulation")
            exit()

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        threads = int(self.remote.free_cores / processes)

        # If hyperthreading should be used for the remote host, we can even use more threads
        if self.remote.use_hyperthreading: threads *= self.remote.threads_per_core

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:

            processes = int(self.remote.free_cores)
            threads = 1

        # Set the parallelization options
        self.config.arguments.parallel.processes = processes
        self.config.arguments.parallel.threads = threads

    # -----------------------------------------------------------------

    def check_parallelization(self):

        """
        This function checks whether the parallelization scheme that is asked by the user is possible given the
        number of cores and hyperthreads per core on the remote host.
        Returns:
        """

        # If the remote host uses a scheduling system, check whether the parallelization options are possible
        # based on the cluster properties defined in the configuration
        if self.remote.scheduler:

            # Determine the total number of hardware threads that can be used on the remote cluster
            hardware_threads_per_node = self.remote.cores
            if self.remote.use_hyperthreading: hardware_threads_per_node *= self.remote.threads_per_core

            # Raise an error if the number of requested threads per process exceeds the number of hardware threads
            # per node
            if self.config.arguments.parallel.threads > hardware_threads_per_node:
                raise RuntimeError("The number of requested threads per process exceeds the number of allowed threads per node")

            # Determine the number of processes per node (this same calculation is also done in JobScript)
            # self.remote.cores = cores per node
            processes_per_node = self.remote.cores // self.config.arguments.parallel.threads

            # Determine the amount of requested nodes based on the total number of processes and the number of processes per node
            requested_nodes = math.ceil(self.config.arguments.parallel.processes / processes_per_node)

            # Raise an error if the number of requested nodes exceeds the number of nodes of the system
            if requested_nodes > self.remote.nodes: raise RuntimeError("The required number of computing nodes for"
                                                                       "the requested number of processes and threads "
                                                                       "exceeds the existing number of nodes")

        # No scheduling system
        else:

            # Determine the total number of requested threads
            requested_threads = self.config.arguments.parallel.processes * self.config.arguments.parallel.threads

            # Determine the total number of hardware threads that can be used on the remote host
            hardware_threads = self.remote.cores
            if self.remote.use_hyperthreading: hardware_threads *= self.remote.threads_per_core

            # If the number of requested threads is greater than the allowed number of hardware threads, raise
            # an error
            if requested_threads > hardware_threads: raise RuntimeError("The requested number of processes and threads "
                                                                        "exceeds the total number of hardware threads")

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the simulation...")

        # Add the walltime to the scheduling options
        if self.config.walltime is not None:
            scheduling_options = SchedulingOptions()
            scheduling_options.walltime = self.config.walltime
        else: scheduling_options = None

        # Run the simulation
        arguments = SkirtArguments(self.config.arguments)
        simulation = self.remote.run(arguments, scheduling_options=scheduling_options)

        # Set the analysis options for the simulation
        self.set_analysis_options(simulation)

        # Save the simulation object
        simulation.save()

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving finished simulations...")

        # Get a list of the simulations that have been succesfully retrieved
        self.simulations = self.remote.retrieve()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the output of retrieved simulations...")

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation)

            # Clear the analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def set_analysis_options(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Set the analysis options from the configuration settings
        simulation.set_analysis_options(self.config.analysis)

        # Remove remote files
        simulation.remove_remote_input = not self.config.keep
        simulation.remove_remote_output = not self.config.keep
        simulation.remove_remote_simulation_directory = not self.config.keep

        # Retrieval
        simulation.retrieve_types = self.config.retrieve_types

# -----------------------------------------------------------------
