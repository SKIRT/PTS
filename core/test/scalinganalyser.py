#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.scalinganalyser Contains the ScalingAnalyser class, used for analysing the scaling results
#  of a SKIRT simulation that is part of a scaling test.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..extract.scaling import ScalingExtractor
from ..plot.scaling import ScalingPlotter
from ..tools.logging import log

# -----------------------------------------------------------------

class ScalingAnalyser(Configurable):

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
        super(ScalingAnalyser, self).__init__(config, "core")

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The timeline and memory usage tables
        self.timeline = None
        self.memory = None

        # Set the scaling table to None initially
        self.scaling = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline, memory, plot=True):

        """
        This function ...
        :return:
        :param simulation:
        :param timeline:
        :param memory:
        :param plot:
        """

        # 1. Call the setup function
        self.setup(simulation, timeline, memory)

        # 2. Extract scaling information
        self.extract()

        # 3. Make the scaling plots
        if plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, simulation, timeline, memory):

        """
        This function ...
        :param simulation:
        :param timeline:
        :param memory:
        :return:
        """

        # Call the setup function of the base class
        super(ScalingAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

        # Make local references to the timeline and memory extractors
        self.timeline = timeline
        self.memory = memory

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the output table to None
        self.scaling = None

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the scaling information ...")

        # Create a ScalingExtractor object
        extractor = ScalingExtractor()

        # Run the scaling extractor
        self.scaling = extractor.run(self.simulation, self.timeline, self.memory)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the scaling information ...")

        # Create a ScalingPlotter object
        plotter = ScalingPlotter()

        # Run the scaling plotter
        plotter.run(self.scaling, self.simulation.analysis.scaling_plot_path)

# -----------------------------------------------------------------
