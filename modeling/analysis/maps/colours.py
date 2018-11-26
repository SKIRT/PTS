#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.colours Contains the ColourMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.colour.colour import ColourMapsMaker, colour_strings
from ....magic.tools.colours import get_filters_for_colour
from ....core.tools.utils import lazyproperty
from ....magic.core.list import FrameList

# -----------------------------------------------------------------

class ColourMapsAnalyser(MapsAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ColourMapsAnalyser, self).__init__(*args, **kwargs)

        # The frames
        self.frames = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # Load the data
        self.load_data()

        # Make the maps
        self.make_maps()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ColourMapsAnalyser, self).setup(**kwargs)

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
                if not self.simulated_dataset.has_frame_for_filter(fltr): break

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
        log.info("Making colour maps ...")

        # Current maps
        if self.config.remake: current = dict()
        else: current = self.current_maps

        # Create the map maker
        maker = ColourMapsMaker()

        # Run the map maker
        maker.run(colours=self.available_colours, frames=self.frames, maps=current)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):
        return self.colours_path

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

        # Write the methods
        self.write_methods()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps()

# -----------------------------------------------------------------
