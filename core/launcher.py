#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.launcher This module can be used to launch SKIRT/FitSKIRT simulations in a convenient way

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import astronomical modules
from astropy import log
import astropy.logger

# Import the relevant PTS modules
from .skirtexec import SkirtExec, FitSkirtExec
from ..core.parameters import SkirtParameters
from ..performance.resources import ResourceEstimator
from ..tools import configuration
from ..extract.progress import ProgressExtractor
from ..extract.timeline import TimeLineExtractor
from ..extract.memory import MemoryExtractor
from ..tools import monitoring
from ..plotting.progress import ProgressPlotter
from ..plotting.timeline import TimeLinePlotter
from ..plotting.memory import MemoryPlotter

# -----------------------------------------------------------------

class SkirtLauncher(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        ## Configuration

        self.config = configuration.set("skirtlauncher", config)

        ## Logging

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ## Attributes

        # Create the SKIRT execution context
        self.skirt = SkirtExec()

        # Set the simulation instance to None initially
        self.simulation = None

        # Set the paths to None initialy
        self.base_path = None
        self.input_path = None
        self.output_path = None
        self.extr_path = None
        self.plot_path = None

        # Tables
        self.progress = None
        self.timeline = None
        self.memory = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :return:
        """

        # Load the default configuration for this class
        config = configuration.set("skirtlauncher")

        ## Adjust the configuration settings according to the command-line arguments

        # Logging (no options here yet)
        # ...

        # Ski file
        config.parameters.ski_pattern = arguments.filepath

        # Simulation logging
        config.parameters.logging.brief = arguments.brief
        config.parameters.logging.verbose = arguments.verbose
        config.parameters.logging.memory = arguments.memory
        config.parameters.logging.allocation = arguments.allocation

        # Other simulation parameters
        config.parameters.emulate = arguments.emulate
        config.parameters.single = True  # For now, we only allow single simulations

        # Extraction
        config.extraction.memory = arguments.extractmemory

        # Plotting
        config.plotting.seds = arguments.plotseds
        config.plotting.grids = arguments.plotgrids
        config.plotting.progress = arguments.plotprogress
        config.plotting.timeline = arguments.plottimeline
        config.plotting.memory = arguments.plotmemory

        # Advanced options
        config.advanced.rgb = arguments.makergb
        config.advanced.wavemovie = arguments.makewave

        # Create and return a SkirtLauncher instance
        return cls(config)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # Setup
        self.setup()

        # Set the parallelization scheme
        self.set_parallelization()

        ## Simulation

        # Run the simulation
        self.simulate()

        ## Extraction

        # Extract the progress information
        if self.config.extraction.progress: self.extract_progress()

        # Extract the timeline information
        if self.config.extraction.timeline: self.extract_timeline()

        # Extract the memory information
        if self.config.extraction.memory: self.extract_memory()

        ## Plotting

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

        ## Advanced

        # If requested, make RGB images of the output FITS files
        if self.config.advanced.rgb: self.make_rgb()

        # If requested, make wave movies from the ouput FITS files
        if self.config.advanced.wavemovie: self.make_wave()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

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
        log.info("Determining the parallelization scheme by estimating the memory requirements...")

        # Calculate the amount of required memory for this simulation
        estimator = ResourceEstimator()
        estimator.run(self.config.parameters.ski_pattern)

        # Calculate the maximum number of processes based on the memory requirements
        processes = int(monitoring.memory() / estimator.memory)

        # If there is too little free memory for the simulation, the number of processes will be smaller than one
        if processes < 1:

            # Exit with an error
            log.error("Not enough memory available to run this simulation")
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
        log.info("Performing the simulation...")

        # Run the simulation
        parameters = SkirtParameters(self.config.parameters)
        self.simulation = self.skirt.run(parameters, silent=True)

    # -----------------------------------------------------------------

    def extract_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the progress information...")

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
        log.info("Extracting the timeline information...")

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
        log.info("Extracting the memory information...")

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
        log.info("Plotting SEDs...")

    # -----------------------------------------------------------------

    def plot_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting grids...")

    # -----------------------------------------------------------------

    def plot_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the progress information...")

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
        log.info("Plotting the timeline...")

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
        log.info("Plotting the memory information...")

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
        log.info("Making RGB images...")

    # -----------------------------------------------------------------

    def make_wave(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making wave movies...")

# -----------------------------------------------------------------

class FitSkirtLauncher(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        ## Configuration

        self.config = configuration.set("fitskirtlauncher", config)

        ## Logging

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ## Attributes

        # Create the FitSKIRT execution context
        self.fitskirt = FitSkirtExec()

        # Set the simulation instance to None initially
        self.simulation = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Load the default configuration for this class
        config = configuration.set("fitskirtlauncher")

        ## Adjust the configuration settings according to the command-line arguments
        # ...

        # Create and return a FitSkirtLauncher instance
        return cls(config)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # Run the simulation
        self.simulate()

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        simulation = self.fitskirt.run()

# -----------------------------------------------------------------
