#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.circle Contains the CircleRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import astronomical modules
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate

# -----------------------------------------------------------------

class CircleRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Set attributes
        self.center = center
        self.radius = radius

        # Call the constructor of the base class
        super(CircleRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    @property
    def circumference(self):

        """
        This function ...
        :return:
        """

        return 2. * math.pi * self.radius

    # -----------------------------------------------------------------

    @property
    def area(self):

        """
        This function ...
        :return:
        """

        return math.pi * self.radius ** 2

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :return:
        """

        return self.__class__(self.center, self.radius * value, self.meta)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :return:
        """

        self.radius *= value
        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :return:
        """

        return self.__class__(self.center, self.radius / value, self.meta)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :return:
        """

        # Divide and return
        self.radius /= value
        return self

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :return:
        """

        return self.__idiv__(value)

# -----------------------------------------------------------------

class PixelCircleRegion(CircleRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be a pixel coordinate")
        if not isinstance(radius, float): raise ValueError("Radius must be a float")

        # Call the constructor of the base class
        super(PixelCircleRegion, self).__init__(center, radius, **kwargs)

# -----------------------------------------------------------------

class SkyCircleRegion(CircleRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be a sky coordinate")
        if not isinstance(radius, Quantity): raise ValueError("Radius must be an angular quantity")

        # Call the constructor of the base class
        super(SkyCircleRegion, self).__init__(center, radius, **kwargs)

# -----------------------------------------------------------------

class PhysicalCircleRegion(CircleRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        This function ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be a physical coordinate")
        if not isinstance(radius, Quantity): raise ValueError("Radius must be a physical quantity of length")

        # Call the constructor of the base class
        super(PhysicalCircleRegion, self).__init__(center, radius, **kwargs)

# -----------------------------------------------------------------
