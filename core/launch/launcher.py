#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module can be used to launch SKIRT/FitSKIRT simulations in a convenient way
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from .analyser import SimulationAnalyser
from ..simulation import SkirtExec
from ..simulation import SkirtArguments
from ..basics import Configurable
from ..test import ResourceEstimator
from ..tools import monitoring

# -----------------------------------------------------------------

class SkirtLauncher(Configurable):

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
        super(SkirtLauncher, self).__init__(config)

        ## Attributes

        # Create the SKIRT execution context
        self.skirt = SkirtExec()

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # Set the simulation instance to None initially
        self.simulation = None

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

        # Create a new SkirtLauncher instance
        launcher = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Logging
        if arguments.debug: launcher.config.logging.level = "DEBUG"

        # Ski file
        launcher.config.arguments.ski_pattern = arguments.filepath

        # Simulation logging
        launcher.config.arguments.logging.brief = arguments.brief
        launcher.config.arguments.logging.verbose = arguments.verbose
        launcher.config.arguments.logging.memory = arguments.memory
        launcher.config.arguments.logging.allocation = arguments.allocation

        # Parallelization
        if arguments.parallel is not None:
            launcher.config.arguments.parallel.processes = arguments.parallel[0]
            launcher.config.arguments.parallel.threads = arguments.parallel[1]

        # Other simulation arguments
        launcher.config.arguments.emulate = arguments.emulate
        launcher.config.arguments.single = True  # For now, we only allow single simulations

        # Extraction
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

        # 4. Analyse the simulation output
        self.analyse()

    # -----------------------------------------------------------------

    @property
    def has_parallelization(self):

        """
        This function ...
        :return:
        """

        # Check whether the number of processes and the number of threads are both defined
        return self.config.parameters.parallel.processes is not None and self.config.parameters.parallel.threads is not None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtLauncher, self).setup()

        # Set the paths
        self.base_path = os.path.dirname(self.config.parameters.ski_pattern) if "/" in self.config.parameters.ski_pattern else os.getcwd()
        self.input_path = os.path.join(self.base_path, "in")
        self.output_path = os.path.join(self.base_path, "out")
        self.extr_path = os.path.join(self.base_path, "extr")
        self.plot_path = os.path.join(self.base_path, "plot")

        # Check if an input directory exists
        if not os.path.isdir(self.input_path): self.input_path = None

        # Set the paths for the simulation
        self.config.parameters.input_path = self.input_path
        self.config.parameters.output_path = self.output_path

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

        # Inform the user
        self.log.info("Determining the parallelization scheme by estimating the memory requirements...")

        # Calculate the amount of required memory for this simulation
        estimator = ResourceEstimator()
        estimator.run(self.config.arguments.ski_pattern)

        # Calculate the maximum number of processes based on the memory requirements
        processes = int(monitoring.free_memory() / estimator.memory)

        # If there is too little free memory for the simulation, the number of processes will be smaller than one
        if processes < 1:

            # Exit with an error
            self.log.error("Not enough memory available to run this simulation")
            exit()

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        threads = int(monitoring.free_cpus() / processes)

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:

            processes = int(monitoring.free_cpus())
            threads = 1

        # Set the parallelization options
        self.config.simulation.parallel.processes = processes
        self.config.simulation.parallel.threads = threads

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
        self.simulation = self.skirt.run(arguments, silent=True)

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Analysing the simulation output...")

        # Set simulation analysis options
        self.simulation.extract_progress = self.config.extraction.progress
        self.simulation.extract_timeline = self.config.extraction.timeline
        self.simulation.extract_memory = self.config.extraction.memory
        self.simulation.plot_seds = self.config.plotting.seds
        self.simulation.plot_grids = self.config.plotting.grids
        self.simulation.plot_progress = self.config.plotting.progress
        self.simulation.plot_timeline = self.config.plotting.timeline
        self.simulation.plot_memory = self.config.plotting.memory
        self.simulation.make_rgb = self.config.advanced.rgb
        self.simulation.make_wave = self.config.advanced.wavemovie

        self.simulation.extraction_path = self.config.extr_path
        self.simulation.plot_path = self.config.plot_path

        # Run the analyser on the simulation
        self.analyser.run(self.simulation)

# -----------------------------------------------------------------
