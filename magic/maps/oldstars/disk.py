#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.oldstars.disk Contains the DiskOldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.basics.configurable import Configurable
from ...core.list import FrameList, NamedFrameList
from ....core.tools.stringify import tostr
from ...core.image import Image

# -----------------------------------------------------------------

def make_map(frame, bulge):

    """
    This function ...
    :param frame:
    :param bulge:
    :return: 
    """

    # Create the maker
    maker = DiskOldStellarMapMaker()

    # Set input
    frames = FrameList(frame)
    bulges = FrameList(bulge)

    # Run the maker
    maker.run(frames=frames, bulges=bulges)

    # Return the map
    return maker.single_map

# -----------------------------------------------------------------

class DiskOldStellarMapMaker(Configurable):

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
        super(DiskOldStellarMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The input
        self.frames = None
        self.bulges = None

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
        :return:
        """

        # Call the setup function of the base class
        super(DiskOldStellarMapMaker, self).setup(**kwargs)

        # Get input
        self.frames = kwargs.pop("frames")
        self.bulges = kwargs.pop("bulges")

        # Get already created maps
        self.maps = kwargs.pop("maps", dict())

        # Get method name
        self.method_name = kwargs.pop("method_name", None)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
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
        log.info("Making the maps of old stars ...")

        # Loop over the filters for which we have a bulge and a total image
        for fltr in self.filters:

            # Set name
            name = tostr(fltr, delimiter="_")

            # Set origin
            self.origins[name] = [fltr]

            # Set methods
            if self.has_method_name: self.methods[name] = [self.method_name]

            # Check if already present
            if name in self.maps:
                log.success("The " + name + " old stellar disk map is already created: not creating it again")
                continue

            # Get the observed frame and the bulge
            frame = self.frames[fltr]
            bulge = self.bulges[fltr]

            # REBIN TO THE SAME PIXELSCALE AND CONVOLVE TO SAME RESOLUTION
            frames = NamedFrameList(observation=frame, bulge=bulge)
            frames.convolve_and_rebin()

            # Subtract bulge from the IRAC I1 image
            minus_bulge = frames["observation"] - frames["bulge"]

            # Remove negatives, replace by zero
            negatives = minus_bulge.replace_negatives(0.0)

            # Normalize the old stellar map
            minus_bulge.normalize()

            # Create image
            image = Image()
            image.add_frame(minus_bulge, "disk")
            image.add_mask(negatives, "negatives")

            # Add
            self.maps[name] = image

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
