#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.ssfr.ssfr Contains the SSFRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent
from ....magic.tools.colours import get_filters_for_colour
from ....core.tools import filesystem as fs

# -----------------------------------------------------------------

ssfr_colours = ["FUV-H", "FUV-i", "FUV-r", "FUV-g", "FUV-B"]

# -----------------------------------------------------------------

class ColoursSSFRMapMaker(MapsComponent):

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
        super(ColoursSSFRMapMaker, self).__init__(config, interactive)

        # -- Attributes --

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

        # 2. Load the colour maps
        self.load_colours()

        # 3. Make maps
        self.make_maps()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ColoursSSFRMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def available_colours(self):

        """
        This function ...
        :return: 
        """

        colours = []

        # Loop over the colours
        for colour in self.config.colours:

            # Get the two filters
            fltr_a, fltr_b = get_filters_for_colour(colour)

            # has_colour_map_for_filters
            if not self.has_colour_map_for_filters(fltr_a, fltr_b): continue

            # Add
            colours.append(colour)

        # Return colours
        return colours

    # -----------------------------------------------------------------

    def load_colours(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the colour maps ...")

        # Loop over the available colours
        for colour in self.available_colours:

            # Debugging
            log.debug("Loading the '" + colour + "' colour map ...")

            # Get the two filters
            fltr_a, fltr_b = get_filters_for_colour(colour)

            # Get the colour map
            colour_map = self.get_colour_map_for_filters(fltr_a, fltr_b)

            # Add the colour map
            self.colours[colour] = colour_map

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

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the colour maps
        self.write_maps()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the sSFR maps ...")

        # Loop over the maps
        for colour in self.maps:

            # Determine path
            path = fs.join(self.maps_ssfr_path, colour + ".fits")

            # Save
            self.maps[colour].saveto(path)

# -----------------------------------------------------------------
