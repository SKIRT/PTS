#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module can be used to analyse SKIRT simulations
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..extract.progress import ProgressExtractor
from ..extract.timeline import TimeLineExtractor
from ..extract.memory import MemoryExtractor
from ..plot.progress import ProgressPlotter
from ..plot.timeline import TimeLinePlotter
from ..plot.memory import MemoryPlotter

# -----------------------------------------------------------------

class SimulationAnalyser(Configurable):

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
        super(SimulationAnalyser, self).__init__(config)

        ## Attributes

        # Set the simulation object to None initially
        self.simulation = None

        # Tables
        self.progress = None
        self.timeline = None
        self.memory = None

    # -----------------------------------------------------------------

    def run(self, simulation):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation)

        # 1. Extract information from the simulation's log files
        self.extract()

        # 2. Make plots based on the simulation output
        self.plot()

        # 3. Advanced output
        self.advanced()

    # -----------------------------------------------------------------

    def setup(self, simulation):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SimulationAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the simulation to None
        self.simulation = None

        # Set the tables to None
        self.progress = None
        self.timeline = None
        self.memory = None

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Extract the progress information
        if self.simulation.extract_progress: self.extract_progress()

        # Extract the timeline information
        if self.simulation.extract_timeline: self.extract_timeline()

        # Extract the memory information
        if self.simulation.extract_memory: self.extract_memory()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # If requested, plot the SED's
        if self.simulation.plot_seds: self.plot_seds()

        # If requested, make plots of the dust grid
        if self.simulation.plot_grids: self.plot_grids()

        # If requested, plot the simulation progress as a function of time
        if self.simulation.plot_progress: self.plot_progress()

        # If requested, plot a timeline of the different simulation phases
        if self.simulation.plot_timeline: self.plot_timeline()

        # If requested, plot the memory usage as a function of time
        if self.simulation.plot_memory: self.plot_memory()

    # -----------------------------------------------------------------

    def advanced(self):

        """
        This function ...
        :return:
        """

        # If requested, make RGB images of the output FITS files
        if self.simulation.make_rgb: self.make_rgb()

        # If requested, make wave movies from the ouput FITS files
        if self.simulation.make_wave: self.make_wave()

    # -----------------------------------------------------------------

    def extract_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the progress information...")

        # Determine the path to the progress file
        path = os.path.join(self.simulation.extraction_path, "progress.dat")

        # Create and run a ProgressExtractor object
        extractor = ProgressExtractor()
        extractor.run(self.simulation, path)

        # Set the table
        self.progress = extractor.table

    # -----------------------------------------------------------------

    def extract_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the timeline information...")

        # Determine the path to the timeline file
        path = os.path.join(self.simulation.extraction_path, "timeline.dat")

        # Create and run a TimeLineExtractor object
        extractor = TimeLineExtractor()
        extractor.run(self.simulation, path)

        # Set the table
        self.timeline = extractor.table

    # -----------------------------------------------------------------

    def extract_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the memory information...")

        # Determine the path to the memory file
        path = os.path.join(self.simulation.extraction_path, "memory.dat")

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
        path = os.path.join(self.simulation.plot_path, "progress.pdf")

        # Create and run a ProgressPlotter object
        plotter = ProgressPlotter()
        plotter.run(self.progress, path)

    # -----------------------------------------------------------------

    def plot_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the timeline...")

        # Determine the path to the timeline plot file
        path = os.path.join(self.simulation.plot_path, "timeline.pdf")

        # Create and run a TimeLinePlotter object
        plotter = TimeLinePlotter()
        plotter.run(self.timeline, path)

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the memory information...")

        # Determine the path to the memory plot file
        path = os.path.join(self.simulation.plot_path, "memory.pdf")

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
