#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.emission Contains the EmissionDustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent

# -----------------------------------------------------------------

class EmissionDustMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(EmissionDustMapMaker, self).__init__()

        # -- Attributes --

        # The dust map
        self.map = None

    # -----------------------------------------------------------------

    @classmethod
    def requirements(cls, config):

        """
        This function ...
        :param config:
        :return:
        """

        config = cls.get_config(config)
        return [str(config.filter)]

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Make the map
        self.make_map()

        # 3. Make everything positive
        self.make_positive()

        # 3. Normalize the map
        self.normalize()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(EmissionDustMapMaker, self).setup()

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Make the map
        self.map = self.dataset.get_frame_for_filter(self.config.filter)

    # -----------------------------------------------------------------

    def make_positive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Replacing NaNs and negative pixels by zeros ...")

        # Set negatives and NaNs to zero
        self.map.replace_nans(0.0)
        self.map.replace_negatives(0.0)

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the dust map ...")

        self.map.normalize()
        self.map.unit = None

# -----------------------------------------------------------------
