#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.composite Contains the CompositeRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion

# -----------------------------------------------------------------

class CompositeRegion(Region):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Set the region elements
        self.elements = args

        # Call the constructor of the base class
        super(CompositeRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelCompositeRegion(CompositeRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check that all elements are pixel regions
        for arg in args:
            if not isinstance(arg, PixelRegion): raise ValueError("All regions of the composite should be in pixel coordinates")

        # Call the constructor of the CompositeRegion base class
        CompositeRegion.__init__(self, *args, **kwargs)

# -----------------------------------------------------------------

class SkyCompositeRegion(CompositeRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check that all elements are sky regions
        for arg in args:
            if not isinstance(arg, SkyRegion): raise ValueError("All regions of the composite should be in sky coordinates")

        # Call the constructor of the CompositeRegion base class
        CompositeRegion.__init__(self, *args, **kwargs)

# -----------------------------------------------------------------

class PhysicalCompositeRegion(CompositeRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check that all elements are physical regions
        for arg in args:
            if not isinstance(arg, PhysicalRegion): raise ValueError("All regions of the composite should be in physical coordinates")

        # Call the constructor of the CompositeRegion base class
        CompositeRegion.__init__(self, *args, **kwargs)

# -----------------------------------------------------------------
