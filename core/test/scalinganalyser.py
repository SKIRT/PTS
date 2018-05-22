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
from ..basics.log import log

# -----------------------------------------------------------------

class ScalingAnalyser(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ScalingAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The timeline and memory usage tables
        self.timeline = None
        self.memory = None

        # Set the scaling table to None initially
        self.scaling = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        :param kwargs:
        """

        # 2. Extract scaling information
        self.extract()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ScalingAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = kwargs.pop("simulation")

        # Make local references to the timeline and memory extractors
        self.timeline = kwargs.pop("timeline")
        self.memory = kwargs.pop("memory")

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
