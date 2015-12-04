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
from ..extract.scaling import ScalingExtractor
from ..plot.scaling import ScalingPlotter

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
        super(ScalingAnalyser, self).__init__(config)

        ## Attributes

        # Set the simulation object to None initially
        self.simulation = None

        # Set the input timeline and memory tables to None initially
        self.timeline = None
        self.memory = None

        # Set the output scaling table to None initially
        self.scaling = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline, memory):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation, timeline, memory)

        # 1. Extract scaling information
        self.extract()

        # 2. Make the scaling plots
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, simulation, timeline, memory):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ScalingAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

        # Make local references to the input tables
        self.timeline = timeline
        self.memory = memory

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the input tables to None
        self.timeline = None
        self.memory = None

        # Set the output table to None
        self.scaling = None

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the scaling information...")

        # Create and run a ScalingExtractor object
        extractor = ScalingExtractor()
        extractor.run(self.timeline, self.memory, self.simulation.scaling_file_path)

        # Set the table
        self.scaling = extractor.table

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the scaling information...")

        # Determine the path to the plot directory
        #path = os.path.join(self.plot_path, ...)

        # Create and run a ScalingPlotter object
        plotter = ScalingPlotter()
        plotter.run(self.scaling, path)

# -----------------------------------------------------------------
