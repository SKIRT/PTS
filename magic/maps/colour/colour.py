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
from ....core.tools.logging import log
from ...tools.colours import make_colour_map, get_filters_for_colour
from ...core.list import FrameList

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

    # Run the map maker
    maker.run(frames=frames)

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

        # The colours
        self.colours = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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

            # Add the origins
            self.origins[colour] = [fltr_a, fltr_b]

    # -----------------------------------------------------------------

    def single_map(self):

        """
        This function ...
        :return: 
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

# -----------------------------------------------------------------
