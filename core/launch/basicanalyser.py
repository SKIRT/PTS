#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.analyser Contains the BasicAnalyser class, used for analysing simulation output.

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
from ..plot.seds import plotseds
from ..plot.grids import plotgrids
from ..plot.rgbimages import makergbimages
from ..plot.wavemovie import makewavemovie
from ..tools.logging import log

# -----------------------------------------------------------------

class BasicAnalyser(Configurable):

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
        super(BasicAnalyser, self).__init__(config, "core")

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The extractors
        self.progress_extractor = ProgressExtractor()
        self.timeline_extractor = TimeLineExtractor()
        self.memory_extractor = MemoryExtractor()

    # -----------------------------------------------------------------

    def run(self, simulation):

        """
        This function ...
        :param simulation
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation)

        # 2. Extract information from the simulation's log files
        self.extract()

        # 3. Make plots based on the simulation output
        self.plot()

        # 4. Advanced output
        self.advanced()

    # -----------------------------------------------------------------

    def setup(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Call the setup function of the base class
        super(BasicAnalyser, self).setup()

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

        # Clear the extractors
        self.progress_extractor.clear()
        self.timeline_extractor.clear()
        self.memory_extractor.clear()

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
        log.info("Extracting the progress information ...")

        # Determine the path to the progress file
        path = os.path.join(self.simulation.extraction_path, "progress.dat")

        # Run the progress extractor
        self.progress_extractor.run(self.simulation, path)

    # -----------------------------------------------------------------

    def extract_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the timeline information ...")

        # Determine the path to the timeline file
        path = os.path.join(self.simulation.extraction_path, "timeline.dat")

        # Run the timeline extractor
        self.timeline_extractor.run(self.simulation, path)

    # -----------------------------------------------------------------

    def extract_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the memory information ...")

        # Determine the path to the memory file
        path = os.path.join(self.simulation.extraction_path, "memory.dat")

        # Run the memory extractor
        self.memory_extractor.run(self.simulation, path)

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs ...")

        # Plot the SEDs for the simulation
        plotseds(self.simulation, output_path=self.simulation.plotting_path)

    # -----------------------------------------------------------------

    def plot_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting grids ...")

        # Plot the dust grid for the simulation
        plotgrids(self.simulation, output_path=self.simulation.plotting_path)

    # -----------------------------------------------------------------

    def plot_progress(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the progress information ...")

        # Create and run a ProgressPlotter object
        plotter = ProgressPlotter()
        plotter.run(self.progress_extractor.table, self.simulation.plotting_path)

    # -----------------------------------------------------------------

    def plot_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the timeline ...")

        # Create and run a TimeLinePlotter object
        plotter = TimeLinePlotter()
        plotter.run(self.timeline_extractor.table, self.simulation.plotting_path)

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the memory information ...")

        # Create and run a MemoryPlotter object
        plotter = MemoryPlotter()
        plotter.run(self.memory_extractor.table, self.simulation.plotting_path)

    # -----------------------------------------------------------------

    def make_rgb(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making RGB images ...")

        # Make RGB images from the output images
        makergbimages(self.simulation, output_path=self.simulation.plotting_path)

    # -----------------------------------------------------------------

    def make_wave(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making wave movies ...")

        # Make wave movies from the output images
        makewavemovie(self.simulation, output_path=self.simulation.plotting_path)

# -----------------------------------------------------------------
