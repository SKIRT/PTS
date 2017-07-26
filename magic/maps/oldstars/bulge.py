#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.oldstars.bulge Contains the BulgeOldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.basics.configurable import Configurable
from ....magic.core.list import FrameList
from ....core.tools.stringify import tostr

# -----------------------------------------------------------------

def make_map(bulge):

    """
    This function ...
    :param bulge: 
    :return: 
    """

    # Create maker
    maker = BulgeOldStellarMapMaker()

    # Set input
    bulges = FrameList(bulge)

    # Make
    maker.run(bulges=bulges)

    # Return the bulge frame
    return maker.single_map

# -----------------------------------------------------------------

def make_maps(**args):

    """
    This function ...
    :param args: 
    :return: 
    """

    # Create maker
    maker = BulgeOldStellarMapMaker()

    # Set input
    bulges = FrameList(**args)

    # Run
    maker.run(bulges=bulges)

    # Return
    return maker.maps

# -----------------------------------------------------------------

class BulgeOldStellarMapMaker(Configurable):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(BulgeOldStellarMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The bulge frames
        self.bulges = None

        # The method name
        self.method_name = None

        # THe maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 4. Make the map of old stars
        self.make_map()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BulgeOldStellarMapMaker, self).setup(**kwargs)

        # Get the input
        self.bulges = kwargs.pop("bulges")

        # Get the method name
        self.method_name = kwargs.pop("method_name", None)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return self.bulges.filters

    # -----------------------------------------------------------------

    @property
    def has_method_name(self):

        """
        This function ...
        :return:
        """

        return self.method_name is not None

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of old stars ...")

        # Loop over the filters
        for fltr in self.filters:

            # Create copy
            bulge = self.bulges[fltr]

            # Normalize
            bulge.normalize()

            # Set name
            name = tostr(fltr, delimiter="_")

            # Add
            self.maps[name] = bulge

            # Set the origin
            self.origins[name] = [fltr]

            # Set methods
            if self.has_method_name: self.methods[name] = [self.method_name]

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
