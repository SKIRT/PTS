#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.polygon Contains the PolygonRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate

# -----------------------------------------------------------------

class PolygonRegion(Region):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Set the points
        self.points = list(args)

        # Call the constructor of the base class
        super(PolygonRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelPolygonRegion(PolygonRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Verify the points
        for arg in args:
            if not isinstance(arg, PixelCoordinate): raise ValueError("Points should be of type 'PixelCoordinate'")

        # Call the constructor of the base class
        super(PixelPolygonRegion, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_point(self, point):

        """
        This function ...
        :return:
        """

        if not isinstance(point, PixelCoordinate): raise ValueError("Point should be of type 'PixelCoordinate'")

        # Add the point
        self.points.append(point)

# -----------------------------------------------------------------

class SkyPolygonRegion(PolygonRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Verify the points
        for arg in args:
            if not isinstance(arg, SkyCoordinate): raise ValueError("Points should be of type 'SkyCoordinate'")

        # Call the constructor of the base class
        super(SkyPolygonRegion, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_point(self, point):

        """
        This function ...
        :param point:
        :return:
        """

        if not isinstance(point, SkyCoordinate): raise ValueError("Point should be of type 'SkyCoordinate'")

        # Add the point
        self.points.append(point)

# -----------------------------------------------------------------

class PhysicalPolygonRegion(PolygonRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Verify the points
        for arg in args:
            if not isinstance(arg, PhysicalCoordinate): raise ValueError("Points should be of type 'PhysicalCoordinate'")

        # Call the constructor of the base class
        super(PhysicalPolygonRegion, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_point(self, point):

        """
        This function ...
        :param point:
        :return:
        """

        if not isinstance(point, PhysicalCoordinate): raise ValueError("Point should be of type 'PhysicalCoordinate'")

        # Add the point
        self.points.append(point)

# -----------------------------------------------------------------
