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
from ...core.dataset import DataSet

# -----------------------------------------------------------------

# This list is not exclusive
colour_strings = ["FUV-NUV", "FUV-H", "FUV-u", "FUV-g", "FUV-r", "FUV-i", "FUV-z", "Pacs 70-Pacs 100", "Pacs 100-Pacs 160",
                  "Pacs 160-SPIRE 250", "SPIRE 250-SPIRE 350", "SPIRE 350-SPIRE 500"]

# -----------------------------------------------------------------

def make_map(modeling_path):

    """
    This function ...
    :return: 
    """

    # Create the colour map maker
    maker = ColourMapsMaker()

    maker.config.check_database = False
    maker.config.colours = [""]
    maker.config.write = False

    maker.config.path = modeling_path

    frames = dict() # indexed on filter

    maker.run(frames=frames)

    # Get the maps
    maps = maker.maps

# -----------------------------------------------------------------

class ColourMapsMaker(Configurable):

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
        super(ColourMapsMaker, self).__init__(config, interactive)

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

        if "frames" in kwargs: self.frames = kwargs.pop("frames")
        elif "dataset" in kwargs: self.load_data(kwargs.pop("dataset"))
        elif self.config.dataset is not None: self.load_data(DataSet.from_file(self.config.dataset))

    # -----------------------------------------------------------------

    def load_data(self, dataset):

        """
        This function ...
        :param dataset:
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        # Loop over the colours
        for colour in self.config.colours:

            # Debugging
            log.debug("Loading frames for the '" + colour + "' colour ...")

            # Get the two filters, load frames
            for fltr in get_filters_for_colour(colour):

                # Already loaded
                if fltr in self.frames: continue

                # Debugging
                log.debug("Loading the '" + str(fltr) + "' frame ...")

                # Load
                frame = dataset.get_frame_for_filter(fltr)
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
        for colour in self.config.colours:

            # Get the two filters
            fltr_a, fltr_b = get_filters_for_colour(colour)

            # Rebin the frames to the same pixelgrid
            frame_a, frame_b = self.rebin_to_highest_pixelscale(self.frames[fltr_a], self.frames[fltr_b])

            # Convert the frames to the same unit
            frame_a, frame_b = self.convert_to_same_unit(frame_a, frame_b, unit="Jy")

            # Create the map
            colour_map = make_colour_map(frame_a, frame_b)

            # Add the map
            self.maps[colour] = colour_map

# -----------------------------------------------------------------
