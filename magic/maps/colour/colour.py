#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.colour.colour Contains the ColourMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from ....core.basics.configurable import Configurable
from ....core.basics.log import log
from ...tools.colours import make_colour_map, get_filters_for_colour
from ...core.list import FrameList
from ...tools import colours

# -----------------------------------------------------------------

uv_colour_strings = ["FUV__NUV"]
fuv_optical_colour_strings = ["FUV__H", "FUV__u", "FUV__g", "FUV__r", "FUV__i", "FUV__z"]
nuv_optical_colour_strings = ["NUV__u", "NUV__g", "NUV__r", "NUV__i", "NUV__z"]
infrared_colour_strings = ["Pacs_70__Pacs_100", "Pacs_100__Pacs_160", "Pacs_160__SPIRE_250", "SPIRE_250__SPIRE_350", "SPIRE_350__SPIRE_500"]

# -----------------------------------------------------------------

colour_strings = uv_colour_strings + fuv_optical_colour_strings + nuv_optical_colour_strings + infrared_colour_strings

# -----------------------------------------------------------------

def make_map(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return: 
    """

    # Create the colour map maker
    maker = ColourMapsMaker()

    # Create frame list
    frames = FrameList(*args) # indexed on filter

    # Check length is 2
    if len(frames) != 2: raise ValueError("Need 2 frames")

    # Determine the colour
    colour = colours.get_colour_name_for_filters(frames[0].filter, frames[1].filter)

    # Run the map maker
    maker.run(frames=frames, colours=[colour])

    # Get the maps
    return maker.single_map

# -----------------------------------------------------------------

class ColourMapsMaker(Configurable):

    """
    This class ...
    """
        
    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ColourMapsMaker, self).__init__(*args, **kwargs)

        # The frames
        self.frames = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

        # The colours
        self.colours = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Make the maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup fucntion of the base class
        super(ColourMapsMaker, self).setup(**kwargs)

        # Get input
        self.frames = kwargs.pop("frames")
        self.colours = kwargs.pop("colours")

        # Get maps that have already been created
        if "maps" in kwargs: self.maps = kwargs.pop("maps")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the colour maps ...")

        # Loop over the colours
        for colour in self.colours:

            # Get the two filters
            fltr_a, fltr_b = get_filters_for_colour(colour)

            # Add the origins
            self.origins[colour] = [fltr_a, fltr_b]

            # Add the methods
            self.methods[colour] = []

            # Check whether a colour map is already present
            if colour in self.maps:
                log.success("The " + colour + " colour map is already created: not creating it again")
                continue

            # Create frame list
            frames = FrameList(self.frames[fltr_a], self.frames[fltr_b])

            # Convolve
            frames.convolve_to_highest_fwhm()

            # Rebin the frames to the same pixelgrid
            frames.rebin_to_highest_pixelscale()

            # Convert the frames to the same unit
            frames.convert_to_same_unit(unit="Jy")

            # Get the frames
            frame_a, frame_b = frames[0], frames[1]

            # Create the map
            colour_map = make_colour_map(frame_a, frame_b)

            # Add the map
            self.maps[colour] = colour_map

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This function ...
        :return: 
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

# -----------------------------------------------------------------
