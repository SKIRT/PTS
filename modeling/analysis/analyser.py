#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.bestmodelanalyser Contains the BestModelAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools.logging import log

# -----------------------------------------------------------------

class BestModelAnalyser(AnalysisComponent):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisComponent, self).__init__(config, interactive)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The log file of the simulation
        self.log_file = None

        # The ski file corresponding to the simulation
        self.ski = None

    # -----------------------------------------------------------------

    def run(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation)

        # 2. Load the log file of the simulation
        self.load_log_file()

        # 6. Write
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the best model analyser ...")

        # Set the attributes to default values
        self.simulation = None
        self.log_file = None
        self.ski = None

    # -----------------------------------------------------------------

    def setup(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Call the setup function of the base class
        super(BestModelAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

        # Load the ski file
        self.ski = self.simulation.parameters()

    # -----------------------------------------------------------------

    def load_log_file(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulation log file ...")

        # Get the log file produced by the simulation
        self.log_file = self.simulation.log_file

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
