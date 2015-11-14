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
import numpy as np
import inspect
import psutil
import multiprocessing

# Import astronomical modules
from astropy import log
import astropy.logger

# Import the relevant PTS modules
from pts.skirtexec import SkirtExec, FitSkirtExec, SkirtMemoryExec
from pts import configuration

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

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "skirtlauncher.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

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

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :return:
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "skirtlauncher.cfg")

        # Open the default configuration
        config = configuration.open(default_config)

        ## Adjust the configuration settings according to the command-line arguments

        # Logging (no options here yet)

        # Simulation input and output
        config.simulation.skifile_path = arguments.filepath
        config.simulation.input_path = arguments.input_path
        config.simulation.output_path = arguments.output_path

        # Simulation logging
        config.simulation.logging.brief = arguments.brief
        config.simulation.logging.verbose = arguments.verbose
        config.simulation.logging.memory = arguments.logmemory

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

        # Set the parallelization scheme
        self.set_scheme()

        # Run the simulation
        self.simulate()

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

        # If requested, make RGB images of the output FITS files
        if self.config.advanced.rgb: self.make_rgb()

        # If requested, make wave movies from the ouput FITS files
        if self.config.advanced.wavemovie: self.make_wave()

    # -----------------------------------------------------------------

    def set_scheme(self):

        """
        This function ...
        :return:
        """

        # Get the total number of processors on this system
        total = multiprocessing.cpu_count()

        # Get the load of the different processors
        load = np.array(psutil.cpu_percent(percpu=True))/100.0

        # Calculate the the number of full processors (where 2 processors with loads x and y contribute as a processor with load x+y)
        full = np.sum(load)

        # Get the number of free processors
        free = total - full

        # Get the currently available virtual memory (in gigabytes)
        memory = psutil.virtual_memory().available / 1e9

        # Inform the user
        #log.debug("The number of currently available processors on this system is " + str(free))
        #log.debug("The amount of currently available memory on this system is " + str(memory) + " gigabytes")

        # Calculate the amount of required memory for this simulation (in gigabytes)
        skirtmemory = SkirtMemoryExec()
        simulation = skirtmemory.run()
        required = simulation.memory()

        # Inform the user
        #log.info("The estimated memory requirement for this simulation is " + str(required) + " gigabytes")

        # Calculate the maximum number of MPI processes
        processes = int(memory/required)

        if processes < 1:

            log.error("Not enough memory available to run this simulation")
            exit()

        # Calculate the number of threads per process
        threads = int(free / processes)

        # Set the parallelization options
        self.config.simulation.parallel.processes = processes
        self.config.simulation.parallel.threads = threads

        # Inform the user
        #log.debug("The number of processes that will be used is " + str(processes))
        #log.debug("The number of threads per process that will be used is " + str(threads))

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Run the simulation
        self.simulation = self.skirt.run(self.config.simulation)

    # -----------------------------------------------------------------

    def plot_seds(self):

        pass

    # -----------------------------------------------------------------

    def plot_grids(self):

        pass

    # -----------------------------------------------------------------

    def plot_progress(self):

        pass

    # -----------------------------------------------------------------

    def plot_timeline(self):

        pass

    # -----------------------------------------------------------------

    def plot_memory(self):

        pass

    # -----------------------------------------------------------------

    def make_rgb(self):

        pass

    # -----------------------------------------------------------------

    def make_wave(self):

        pass

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

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "fitskirtlauncher.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

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

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "fitskirtlauncher.cfg")

        # Open the default configuration
        config = configuration.open(default_config)

        # Adjust the configuration settings according to the command-line arguments


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

class SkirtMemoryLauncher(object):

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

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "skirtmemorylauncher.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        ## Logging

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ## Attributes

        # Create the SkirtMemory execution context
        self.skirtmemory = SkirtMemoryExec()

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

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "skirtmemorylauncher.cfg")

        # Open the default configuration
        config = configuration.open(default_config)

        ## Adjust the configuration settings according to the command-line arguments

        # Logging (no options here yet)

        # Simulation input and output
        config.simulation.skifile_path = arguments.filepath
        config.simulation.input_path = arguments.input_path
        config.simulation.output_path = arguments.output_path

        # Plotting
        config.plotting.memory = arguments.plotmemory

        # Create and return a SkirtMemoryLauncher instance
        return cls(config)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # Run the simulation
        self.simulate()

        # If requested, plot the memory usage as a function of time
        if self.config.plotting.memory: self.plot_memory()

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Run the simulation
        self.simulation = self.skirtmemory.run(self.config.simulation)

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------