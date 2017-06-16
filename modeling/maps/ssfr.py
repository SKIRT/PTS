#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.ssfr Contains the SSFRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ...magic.maps.ssfr.colours import ColoursSSFRMapsMaker, ssfr_colours

# -----------------------------------------------------------------

class SSFRMapMaker(MapsComponent):

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
        super(SSFRMapMaker, self).__init__(*args, **kwargs)

        # The colour maps
        self.colours = dict()

        # THe maps
        self.maps = dict()

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_ssfr_path

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Load the colour maps
        self.load_colours()

        # 2. Make SSFR maps based on colours
        self.make_ssfr_colours()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SSFRMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the colour maps ...")

        # Loop over the possible colours for tracing sSFR
        for colour in ssfr_colours:

            # Get colour map and the name
            colour_map, colour_name = self.get_colour_map_and_name(colour)

            # Check if found
            if colour_map is None:
                log.warning("Could not find a '" + colour + "' colour map")
                continue

            # Add the colour map
            self.colours[colour_name] = colour_map

    # -----------------------------------------------------------------

    def make_ssfr_colours(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Making sSFR maps based on colours ...")

        # Create the map maker
        maker = ColoursSSFRMapsMaker()

        # Run the maker
        maker.run(colours=self.colours, colours_origins=self.get_colours_origins(), maps=self.current_maps)

        # Get the maps
        self.maps = maker.maps

        # Get the origins
        self.origins = maker.origins

    # -----------------------------------------------------------------

    def write(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the origins
        self.write_origins()

# -----------------------------------------------------------------
