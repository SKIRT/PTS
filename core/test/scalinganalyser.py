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

        # Set the ski file path to None initially
        self.ski_path = None

        # Tables
        self.scaling = None

    # -----------------------------------------------------------------

    def run(self, ski_path):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(ski_path)

        # 1. Extract scaling information from the simulations' log files
        self.extract()

        # 2. Make the scaling plots
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, ski_path):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ScalingAnalyser, self).setup()

        # Make a local reference to the ski file path
        self.ski_path = ski_path

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the ski file path to None
        self.ski_path = None

        # Set the scaling table to None
        self.table = None

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Extracting the scaling information...")

        # Determine the path to the scaling data file
        path = os.path.join(self.extraction_path, "scaling.dat")

        # Create and run a ScalingExtractor object
        extractor = ScalingExtractor()
        extractor.run(self.ski_path, path)

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
