#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.oldstars.total Contains the TotalOldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.basics.configurable import Configurable
from ....core.tools.stringify import tostr

# -----------------------------------------------------------------

def make_map(frame):

    """
    This function ...
    :param frame: 
    :return: 
    """

    pass

# -----------------------------------------------------------------

class TotalOldStellarMapMaker(Configurable):

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
        super(TotalOldStellarMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The frames
        self.frames = None

        # The method name
        self.method_name = None

        # The maps
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

        # 2. Make the map of old stars
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TotalOldStellarMapMaker, self).setup(**kwargs)

        # Get the frames
        self.frames = kwargs.pop("frames")

        # Get the already existing maps
        self.maps = kwargs.pop("maps", dict())

        # Get the method name
        self.method_name = kwargs.pop("method_name", None)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This fucntion ...
        :return: 
        """

        return self.frames.filters

    # -----------------------------------------------------------------

    @property
    def has_method_name(self):

        """
        This function ...
        :return:
        """

        return self.method_name is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of old stars ...")

        # Loop over the filters
        for fltr in self.filters:

            # Set name
            name = tostr(fltr, delimiter="_")

            # Set origin
            self.origins[name] = [fltr]

            # Set method
            if self.has_method_name: self.methods[name] = [self.method_name]

            # Check if already present
            if name in self.maps:
                log.success("The " + name + " total old stellar map is already created: not creating it again")
                continue

            # Set map
            total = self.frames[fltr]

            # Normalize
            total.normalize()

            # Add
            self.maps[name] = total

# -----------------------------------------------------------------
