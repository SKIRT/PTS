#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelanalyser Contains the ModelAnalyser class, used for analysing the goodness
#  of the radiative transfer model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools.logging import log

# -----------------------------------------------------------------

class ModelAnalyser(ModelingComponent):

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
        super(ModelAnalyser, self).__init__(config)

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

    # -----------------------------------------------------------------

    def run(self, simulation):

        """
        This function ...
        :param simulation
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation)

    # -----------------------------------------------------------------

    def setup(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Call the setup function of the base class
        super(ModelAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the model analyser ...")

        # Set the simulation to None
        self.simulation = None

# -----------------------------------------------------------------
