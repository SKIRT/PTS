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

        # Set the timeline and memory extractors to None initially
        self.te = None
        self.me = None

        # Set the scaling table to None initially
        self.scaling = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline_extractor, memory_extractor, plot=True):

        """
        This function ...
        :return:
        :param simulation:
        :param timeline_extractor:
        :param memory_extractor:
        """

        # 1. Call the setup function
        self.setup(simulation, timeline_extractor, memory_extractor)

        # 2. Extract scaling information
        self.extract()

        # 3. Make the scaling plots
        if plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, simulation, timeline_extractor, memory_extractor):

        """
        This function ...
        :param simulation:
        :param timeline_extractor:
        :param memory_extractor:
        :return:
        """

        # Call the setup function of the base class
        super(ScalingAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

        # Make local references to the timeline and memory extractors
        self.te = timeline_extractor
        self.me = memory_extractor

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
        log.info("Extracting the scaling information...")

        # Create and run a ScalingExtractor object
        extractor = ScalingExtractor()
        extractor.run(self.simulation, self.te, self.me)

        # Set the table
        self.scaling = extractor.table

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the scaling information...")

        # Create and run a ScalingPlotter object
        plotter = ScalingPlotter()
        plotter.run(self.scaling, self.simulation.scaling_plot_path)

# -----------------------------------------------------------------
