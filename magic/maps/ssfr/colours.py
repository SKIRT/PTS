#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.ssfr.colours Contains the ColoursSSFRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from copy import copy

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

ssfr_colours = ["FUV-H", "FUV-i", "FUV-r", "FUV-g", "FUV-B", "NUV-i", "NUV-r", "NUV-g", "NUV-B"]

# -----------------------------------------------------------------

def make_map(**kwargs):

    """
    This function ...
    :param:
    :return: 
    """

    # Get smoothing factor
    smoothing_factor = kwargs.pop("smoothing_factor", None)

    # Create the sSFR map maker
    maker = ColoursSSFRMapsMaker()

    # Set smoothing factor
    if smoothing_factor is not None:
        maker.config.smooth = True
        maker.config.smoothing_factor = smoothing_factor

    # Make
    maker.run(colours=kwargs)

    # Return the single colour map based on the input
    return maker.single_map

# -----------------------------------------------------------------

class ColoursSSFRMapsMaker(Configurable):

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
        super(ColoursSSFRMapsMaker, self).__init__(*args, **kwargs)

        # Input
        self.colours = dict()
        self.colours_origins = None
        self.colours_methods = None

        # The sSFR maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 3. Make maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ColoursSSFRMapsMaker, self).setup(**kwargs)

        # Get the colours
        self.colours = kwargs.pop("colours")

        # Get origins
        self.colours_origins = kwargs.pop("colours_origins", None)

        # Get methods
        self.colours_methods = kwargs.pop("colours_methods", None)

        # Get method name
        self.method_name = kwargs.pop("method_name", None)
        if self.has_methods and self.method_name is None: raise ValueError("Method name has to be specified when methods are given")

        # Get maps that have already been created
        if "maps" in kwargs: self.maps = kwargs.pop("maps")

    # -----------------------------------------------------------------

    @property
    def has_origins(self):
        return self.colours_origins is not None

    # -----------------------------------------------------------------

    @property
    def has_methods(self):
        return self.colours_methods is not None

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

            # Set origin
            # Set the origins
            if self.has_origins:

                # Set the origins
                origins = copy(self.colours_origins[colour])
                self.origins[colour] = origins

            # Set methods
            if self.has_methods:

                methods = copy(self.colours_methods[colour])
                methods.append(self.method_name)
                self.methods[colour] = methods

            # Check whether a colour map is already present
            if colour in self.maps:
                log.success("The " + colour + " sSFR map is already created: not creating it again")
                continue

            # Get the map
            colour_map = self.colours[colour]

            # Replace infinities by NaN
            colour_map.replace_infs(float("nan"))

            # Smooth?
            if self.config.smooth: colour_map.smooth(self.config.smoothing_factor, preserve_nans=False)

            # Set as sSFR map
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
