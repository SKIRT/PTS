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
import os
import json

# Import the relevant PTS classes and modules
from .analyser import SimulationAnalyser
from ..simulation.remote import SkirtRemote
from ..simulation.arguments import SkirtArguments
from ..basics.configurable import Configurable
from ..test.resources import ResourceEstimator

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

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :return:
        """

        # Create and a new SkirtLauncher instance
        launcher = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Logging
        if arguments.debug: launcher.config.logging.level = "DEBUG"

        # Remote host and cluster (if applicable)
        launcher.config.remote = arguments.remote
        launcher.config.cluster = arguments.cluster

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
        launcher.config.extraction.progress = arguments.extractprogress
        launcher.config.extraction.timeline = arguments.extracttimeline
        launcher.config.extraction.memory = arguments.extractmemory

        # Plotting
        launcher.config.plotting.seds = arguments.plotseds
        launcher.config.plotting.grids = arguments.plotgrids
        launcher.config.plotting.progress = arguments.plotprogress
        launcher.config.plotting.timeline = arguments.plottimeline
        launcher.config.plotting.memory = arguments.plotmemory

        # Advanced options
        launcher.config.advanced.rgb = arguments.makergb
        launcher.config.advanced.wavemovie = arguments.makewave

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
        self.base_path = os.path.dirname(self.config.arguments.ski_pattern) if "/" in self.config.arguments.ski_pattern else os.getcwd()
        self.input_path = os.path.join(self.base_path, "in")
        self.output_path = os.path.join(self.base_path, "out")
        self.extr_path = os.path.join(self.base_path, "extr")
        self.plot_path = os.path.join(self.base_path, "plot")

        # Check if an input directory exists
        if not os.path.isdir(self.input_path): self.input_path = None

        # Set the paths for the simulation
        self.config.arguments.input_path = self.input_path
        self.config.arguments.output_path = self.output_path

        # Create the output directory if necessary
        if not os.path.isdir(self.output_path): os.makedirs(self.output_path)

        # Create the extraction directory if necessary
        if self.config.extraction.progress or self.config.extraction.timeline or self.config.extraction.memory:
            if not os.path.isdir(self.extr_path): os.makedirs(self.extr_path)

        # Create the plotting directory if necessary
        if self.config.plotting.seds or self.config.plotting.grids or self.config.plotting.progress \
            or self.config.plotting.timeline or self.config.plotting.memory:
            if not os.path.isdir(self.plot_path): os.makedirs(self.plot_path)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        #self.log.info("free cores: " + str(self.remote.free_cores))
        #self.log.info("free memory: " + str(self.remote.free_memory))
        #self.log.info("free space: " + str(self.remote.free_space))
        #self.log.info("cores: " + str(self.remote.cores))
        #self.log.info("cpu load: " + str(self.remote.cpu_load))
        #self.log.info("memory load: " + str(self.remote.memory_load))

        # Inform the user
        self.log.info("Determining the parallelization scheme by estimating the memory requirements...")

        # Calculate the amount of required memory for this simulation
        estimator = ResourceEstimator()
        estimator.run(self.config.arguments.ski_pattern)

        # Calculate the maximum number of processes based on the memory requirements
        processes = int(self.remote.free_memory / estimator.memory)

        # If there is too little free memory for the simulation, the number of processes will be smaller than one
        if processes < 1:

            # Exit with an error
            self.log.error("Not enough memory available to run this simulation")
            exit()

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        threads = int(self.remote.free_cores / processes)

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:

            processes = int(self.remote.free_cores)
            threads = 1

        # Set the parallelization options
        self.config.arguments.parallel.processes = processes
        self.config.arguments.parallel.threads = threads

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Performing the simulation...")

        # Run the simulation
        arguments = SkirtArguments(self.config.arguments)
        simulation = self.remote.run(arguments)

        # Add additional information to the simulation object
        self.add_analysis_info(simulation)

        # Save the simulation object
        simulation.save()

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Get a list of the simulations that have been succesfully retrieved
        self.simulations = self.remote.retrieve()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation)

            # Clear the analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def add_analysis_info(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Extraction
        simulation.extraction.progress = self.config.extraction.progress
        simulation.extraction.timeline = self.config.extraction.timeline
        simulation.extraction.memory = self.config.extraction.memory
        simulation.extraction_path = self.extr_path

        # Plotting
        simulation.plotting.progress = self.config.plotting.progress
        simulation.plotting.timeline = self.config.plotting.timeline
        simulation.plotting.memory = self.config.plotting.memory
        simulation.plotting.seds = self.config.plotting.seds
        simulation.plotting.grids = self.config.plotting.grids
        simulation.plotting_path = self.plot_path

        # Advanced
        simulation.advanced.rgb = self.config.advanced.rgb
        simulation.advanced.wave = self.config.advanced.wavemovie

        # Remove remote files
        simulation.remove_remote_input = not self.config.keep
        simulation.remove_remote_output = not self.config.keep
        simulation.remove_remote_simulation_directory = not self.config.keep

        # Retrieval
        simulation.retrieve_types = self.config.retrieve_types

# -----------------------------------------------------------------
