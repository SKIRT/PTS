#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.colour.colour Contains the ColourMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.utils import lazyproperty

# Import standard modules
from ..component import MapsComponent
from ....core.tools.logging import log
from ....core.filter.filter import parse_filter
from ....magic.core.frame import Frame
from ....core.tools import filesystem as fs

# -----------------------------------------------------------------

# This list is not exclusive
colour_strings = ["FUV-NUV", "FUV-H", "FUV-u", "FUV-g", "FUV-r", "FUV-i", "FUV-z", "Pacs 70-Pacs 100", "Pacs 100-Pacs 160",
                  "Pacs 160-SPIRE 250", "SPIRE 250-SPIRE 350", "SPIRE 350-SPIRE 500"]

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
        self.maps = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Setup
        self.setup(**kwargs)

        # Load the data
        self.load_data()

        # Make the maps
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
        for colour in self.config.colours:

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

        # Loop over the colours
        for colour in self.available_colours:

            # Get the two filters
            fltr_a, fltr_b = get_filters_for_colour(colour)

            # Rebin the frames to the same pixelgrid
            frame_a, frame_b = self.rebin_to_highest_pixelscale(self.frames[fltr_a], self.frames[fltr_b])

            # Convert the frames to the same unit
            frame_a, frame_b = self.convert_to_same_unit(frame_a, frame_b)

            # Create the map
            colour_map = make_colour_map(frame_a, frame_b)

            # Add the map
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
        log.info("Writing the colour maps ...")

        # Loop over the colours
        for colour in self.available_colours:

            # Determine the path
            path = fs.join(self.maps_colours_path, colour + ".fits")

            # Save the map
            self.maps[colour].saveto(path)

# -----------------------------------------------------------------

def get_filters_for_colour(colour):

    """
    This function ...
    :param colour:
    :return: 
    """

    str_a, str_b = colour.split("-")
    fltr_a, fltr_b = parse_filter(str_a), parse_filter(str_b)
    return fltr_a, fltr_b

# -----------------------------------------------------------------

def make_colour_map(frame_a, frame_b):

    """
    This function ...
    :param frame_a:
    :param frame_b:
    :return:
    """

    return Frame(-2.5 * np.log10(frame_a / frame_b), wcs=frame_a.wcs)

# -----------------------------------------------------------------

