#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.oldstars.total Contains the TotalOldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.basics.configurable import Configurable
from ...core.list import FrameList

# -----------------------------------------------------------------

def make_map(frame):

    """
    This function ...
    :param frame: 
    :return: 
    """

    pass

# -----------------------------------------------------------------

class TotalOldStellarMapMaker(Configurable):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(TotalOldStellarMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The frames
        self.frames = None

        # THe maps
        self.maps = dict()

        # The origins
        self.origins = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 4. Make the map of old stars
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TotalOldStellarMapMaker, self).setup(**kwargs)

        # Get the frames
        self.frames = kwargs.pop("frames")

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This fucntion ...
        :return: 
        """

        return self.frames.filters

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of old stars ...")

        # Loop over the filters
        for fltr in self.filters:

            # Set name
            name = str(fltr)

            # Set map
            total = self.frames[fltr]

            # Normalize
            total.normalize()

            # Add
            self.maps[name] = total

            # Set origin
            self.origins[name] = [fltr]

# -----------------------------------------------------------------
