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
from ...magic.maps.colour.colour import ColourMapsMaker
from ...magic.tools.colours import get_filters_for_colour
from ...magic.core.list import FrameList

# -----------------------------------------------------------------

colour_strings = ["FUV-NUV", "FUV-H", "FUV-u", "FUV-g", "FUV-r", "FUV-i", "FUV-z", "Pacs 70-Pacs 100", "Pacs 100-Pacs 160",
                  "Pacs 160-SPIRE 250", "SPIRE 250-SPIRE 350", "SPIRE 350-SPIRE 500"]

# -----------------------------------------------------------------

class ColourMapMaker(MapsComponent):

    """
    This class ...
    """
        
    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ColourMapMaker, self).__init__(*args, **kwargs)

        # The frames
        self.frames = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_colours_path

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

        # 4. Write
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

        # Load the current maps
        self.current_maps = self.get_current_maps()

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

        # Create frame list
        self.frames = FrameList()

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
                frame = self.get_frame_for_filter(fltr)
                self.frames.append(frame, fltr)

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

        # Run the map maker
        maker.run(colours=self.available_colours, frames=self.frames, maps=self.current_maps)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

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

        # Write origins
        self.write_origins()

# -----------------------------------------------------------------
