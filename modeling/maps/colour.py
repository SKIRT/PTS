#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.colour Contains the ColourMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import standard modules
from .component import MapsComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...magic.maps.colour.colour import ColourMapsMaker, colour_strings
from ...magic.tools.colours import get_filters_for_colour

# -----------------------------------------------------------------

class ColourMapMaker(MapsComponent):

    """
    This class ...
    """
        
    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(ColourMapMaker, self).__init__(config, interactive)

        # The frames
        self.frames = dict()

        # The maps
        self.maps = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the data
        self.load_data()

        # 3. Make the maps
        self.make_maps()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup fucntion of the base class
        super(ColourMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def available_colours(self):

        """
        This function ...
        :return: 
        """

        colours = []

        # Loop over the colours
        for colour in colour_strings:

            # Get the two filters
            for fltr in get_filters_for_colour(colour):

                # If either one of the two images is not available, we can not calculate the colour
                if not self.dataset.has_frame_for_filter(fltr): break

            # Break not encountered
            else: colours.append(colour)

        # Return colours
        return colours

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        # Loop over the colours
        for colour in self.available_colours:

            # Debugging
            log.debug("Loading frames for the '" + colour + "' colour ...")

            # Get the two filters, load frames
            for fltr in get_filters_for_colour(colour):

                # Already loaded
                if fltr in self.frames: continue

                # Debugging
                log.debug("Loading the '" + str(fltr) + "' frame ...")

                # Load
                frame = self.dataset.get_frame_for_filter(fltr)
                self.frames[fltr] = frame

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the colour maps ...")

        # Create the map maker
        maker = ColourMapsMaker()

        # Set the colours
        maker.config.colours = self.available_colours

        # Run the map maker
        maker.run()

        # Set the maps
        self.maps = maker.maps

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

