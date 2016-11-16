#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.rectangle Contains the RectangleRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Quantity
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .region import Region
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch

# -----------------------------------------------------------------

class RectangleRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the angle
        if angle is None: angle = Angle(0., "deg")
        elif not isinstance(angle, Angle): raise ValueError("Angle must be an Astropy Angle object")

        # Set the attributes
        self.center = center
        self.radius = radius
        self.angle = angle

        # Call the constructor of the base class
        super(RectangleRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelRectangleRegion(RectangleRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be a pixel coordinate")
        if not isinstance(radius, PixelStretch): raise ValueError("Radius must be a pixel stretch")

        # Call the constructor of the base class
        super(PixelRectangleRegion, self).__init__(center, radius, angle, **kwargs)

# -----------------------------------------------------------------

class SkyRectangleRegion(RectangleRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be a sky coordinate")
        if not isinstance(radius, SkyStretch): raise ValueError("Radius must be a sky stretch")

        # Call the constructor of the base class
        super(SkyRectangleRegion, self).__init__(center, radius, angle, **kwargs)

# -----------------------------------------------------------------

class PhysicalRectangleRegion(RectangleRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be a physical coordinate")
        if not isinstance(radius, PhysicalStretch): raise ValueError("Radius must be a physical stretch")

        # Call the constructor of the base class
        super(PhysicalRectangleRegion, self).__init__(center, radius, angle, **kwargs)

# -----------------------------------------------------------------
