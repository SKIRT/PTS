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
from ..simulation import SkirtRemote
from ..simulation import SkirtParameters
from ..basics import Configurable
from ..test import ResourceEstimator
from ..extract import ProgressExtractor, TimeLineExtractor, MemoryExtractor
from ..plot import ProgressPlotter, TimeLinePlotter, MemoryPlotter
from ..tools import monitoring

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
        super(SkirtRemoteLauncher, self).__init__(config)

        ## Attributes

        # Create the SKIRT remote execution context
        self.remote = SkirtRemote()

        # Initialize a list to contain the retreived finished simulations
        self.simulations = []

        # Set the paths to None initialy
        self.base_path = None
        self.input_path = None
        self.output_path = None
        self.extr_path = None
        self.plot_path = None

        # Tables
        #self.progress = None
        #self.timeline = None
        #self.memory = None

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

        # Logging (no options here yet)
        # ...

        # Ski file
        launcher.config.parameters.ski_pattern = arguments.filepath

        # Simulation logging
        launcher.config.parameters.logging.brief = arguments.brief
        launcher.config.parameters.logging.verbose = arguments.verbose
        launcher.config.parameters.logging.memory = arguments.memory
        launcher.config.parameters.logging.allocation = arguments.allocation

        # Other simulation parameters
        launcher.config.parameters.emulate = arguments.emulate
        launcher.config.parameters.single = True  # For now, we only allow single simulations

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
        self.set_parallelization()

        # 3. Run the simulation
        self.simulate()

        # 4. Retrieve the simulations that are finished
        self.retreive()

        # 5. Extract information from the simulation's log files
        self.extract()

        # 6. Make plots based on the simulation output
        self.plot()

        # 7. Advanced output
        self.advanced()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtRemoteLauncher, self).setup()

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
        estimator.run(self.config.parameters.ski_pattern)

        # Calculate the maximum number of processes based on the memory requirements
        processes = int(monitoring.memory() / estimator.memory)

        # If there is too little free memory for the simulation, the number of processes will be smaller than one
        if processes < 1:

            # Exit with an error
            self.log.error("Not enough memory available to run this simulation")
            exit()

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        threads = int(monitoring.cpu() / processes)

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:

            processes = int(monitoring.cpu())
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
        parameters = SkirtParameters(self.config.parameters)
        self.remote.run(parameters)

    # -----------------------------------------------------------------

    def retreive(self):

        """
        This function ...
        :return:
        """

        # Get a list of the simulations that have been succesfully retreived
        self.simulations = self.remote.retreive()

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Loop over the different retreived simulations
        for simulation in self.simulations:

            # Extract the progress information
            if self.config.extraction.progress: self.extract_progress()

            # Extract the timeline information
            if self.config.extraction.timeline: self.extract_timeline()

            # Extract the memory information
            if self.config.extraction.memory: self.extract_memory()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Loop over the different retreived simulations
        for simulation in self.simulations:

            # If requested, plot the SED's
            if self.config.plotting.seds: self.plot_seds()

            # If requested, make plots of the dust grid
            if self.config.plotting.grids: self.plot_grids()

            # If requested, plot the simulation progress as a function of time
            if self.config.plotting.progress: self.plot_progress()

            # If requested, plot a timeline of the different simulation phases
            if self.config.plotting.timeline: self.plot_timeline()

            # If requested, plot the memory usage as a function of time
            if self.config.plotting.memory: self.plot_memory()

    # -----------------------------------------------------------------

    def advanced(self):

        """
        This function ...
        :return:
        """

        # Loop over the different retreived simulations
        for simulation in self.simulations:

            # If requested, make RGB images of the output FITS files
            if self.config.advanced.rgb: self.make_rgb()

            # If requested, make wave movies from the ouput FITS files
            if self.config.advanced.wavemovie: self.make_wave()

    # -----------------------------------------------------------------

    def extract_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the progress information...")

        # Determine the path to the progress file
        path = os.path.join(self.extr_path, "progress.dat")

        # Create and run a ProgressExtractor object
        extractor = ProgressExtractor()
        extractor.run(self.simulation, path)

        # Set the table
        #self.progress = extractor.table

    # -----------------------------------------------------------------

    def extract_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the timeline information...")

        # Determine the path to the timeline file
        path = os.path.join(self.extr_path, "timeline.dat")

        # Create and run a TimeLineExtractor object
        extractor = TimeLineExtractor()
        extractor.run(self.simulation, path)

        # Set the table
        #self.timeline = extractor.table

    # -----------------------------------------------------------------

    def extract_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the memory information...")

        # Determine the path to the memory file
        path = os.path.join(self.extr_path, "memory.dat")

        # Create and run a MemoryExtractor object
        extractor = MemoryExtractor()
        extractor.run(self.simulation, path)

        # Set the table
        self.memory = extractor.table

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting SEDs...")

    # -----------------------------------------------------------------

    def plot_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting grids...")

    # -----------------------------------------------------------------

    def plot_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the progress information...")

        # Determine the path to the progress plot file
        path = os.path.join(self.plot_path, "progress.pdf")

        # Create and run a ProgressPlotter object
        plotter = ProgressPlotter()
        #plotter.run(self.progress, path)

    # -----------------------------------------------------------------

    def plot_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the timeline...")

        # Determine the path to the timeline plot file
        path = os.path.join(self.plot_path, "timeline.pdf")

        # Create and run a TimeLinePlotter object
        #plotter = TimeLinePlotter()
        #plotter.run(self.timeline, path)

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the memory information...")

        # Determine the path to the memory plot file
        path = os.path.join(self.plot_path, "memory.pdf")

        # Create and run a MemoryPlotter object
        plotter = MemoryPlotter()
        plotter.run(self.memory, path)

    # -----------------------------------------------------------------

    def make_rgb(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Making RGB images...")

    # -----------------------------------------------------------------

    def make_wave(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Making wave movies...")

# -----------------------------------------------------------------
