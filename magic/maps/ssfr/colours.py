#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.ssfr.colours Contains the ColoursSSFRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.tools import filesystem as fs
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

ssfr_colours = ["FUV-H", "FUV-i", "FUV-r", "FUV-g", "FUV-B"]

# -----------------------------------------------------------------

def make_map():

    """
    This function ...
    :param:
    :return: 
    """

    # Create the sSFR map maker
    maker = ColoursSSFRMapsMaker()

    colours = dict()

    maker.run(colours=colours)

    # Get the maps
    maps = maker.maps

# -----------------------------------------------------------------

class ColoursSSFRMapsMaker(Configurable):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ColoursSSFRMapsMaker, self).__init__(*args, **kwargs)

        # The colour maps
        self.colours = dict()

        # The sSFR maps
        self.maps = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 3. Make maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ColoursSSFRMapsMaker, self).setup(**kwargs)

        # Get the colours
        self.colours = kwargs.pop("colours")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the sSFR maps ...")

        # Loop over the colour maps
        for colour in self.colours:

            # Get the map
            colour_map = self.colours[colour]

            # Set as sSFR map
            self.maps[colour] = colour_map

# -----------------------------------------------------------------
