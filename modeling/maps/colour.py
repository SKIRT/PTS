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

# Import standard modules
from .component import MapsComponent
from ...core.basics.log import log
from ...magic.maps.colour.colour import ColourMapsMaker, colour_strings
from ...magic.tools.colours import get_filters_for_colour
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

methods = None

# -----------------------------------------------------------------

class ColoursMapMaker(MapsComponent):

    """
    This class ...
    """
        
    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ColoursMapMaker, self).__init__(*args, **kwargs)

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

        # 5. Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup fucntion of the base class
        super(ColoursMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def all_colour_filters(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the colours
        for colour in colour_strings: filters.extend(get_filters_for_colour(colour))

        # Return
        return filters

    # -----------------------------------------------------------------

    @lazyproperty
    def available_colours(self):

        """
        This function ...
        :return: 
        """

        colours = []

        # Get names
        names = self.dataset.get_names_for_filters(self.all_colour_filters, as_dict=True)

        # Loop over the colours
        for colour in colour_strings:

            # Get the two filters
            for fltr in get_filters_for_colour(colour):

                # If either one of the two images is not available, we can not calculate the colour
                #if not self.dataset.has_frame_for_filter(fltr): break
                if fltr not in names: break

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

        # List of filters for which we need to load the data
        filters = []

        # Loop over the colours
        for colour in self.available_colours:

            # Debugging
            log.debug("Loading frames for the '" + colour + "' colour ...")

            # Get the two filters, load frames
            for fltr in get_filters_for_colour(colour):

                # Already loaded
                #if fltr in self.frames: continue
                if fltr in filters: continue

                # Add the filter
                filters.append(fltr)

        # Get frames
        self.frames = self.get_frames_for_filters(filters, framelist=True)

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the colour maps ...")

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

    @property
    def scales(self):

        """
        This function ...
        :return:
        """

        scales = dict()

        from ...magic.tools.colours import is_fir_colour, is_fir_or_submm_colour

        # Loop over the colours
        for colour_name in self.maps:

            # Check whether it is a FIR colour
            if is_fir_colour(colour_name) or is_fir_or_submm_colour(colour_name):
                #around_zero = True
                scale = "linear"
            else:
                #around_zero = False
                scale = "log"

            # Set the scale
            scales[colour_name] = scale

        # Return the scales
        return scales

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps(scales=self.scales, share_limits=False)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

# -----------------------------------------------------------------
