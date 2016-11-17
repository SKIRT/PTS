#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.text Contains the TextRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .rectangle import PixelRectangleRegion, SkyRectangleRegion
from ..basics.mask import Mask

# -----------------------------------------------------------------

class VectorRegion(Region):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        The constructor ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the angle
        if not isinstance(angle, Angle): raise ValueError("Angle must be a Astropy Angle object")

        # Set the attributes
        self.start = start
        self.length = length
        self.angle = angle

        # Call the constructor of the base class
        super(VectorRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    @property
    def axis1_center(self):

        """
        This property ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    @property
    def axis2_center(self):
        """
        This property ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

# -----------------------------------------------------------------

class PixelVectorRegion(VectorRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        if not isinstance(start, PixelCoordinate): raise ValueError("Start must be pixel coordinate")

        # Check the length
        if not isinstance(length, float): raise ValueError("Length must be float")

        # Call the constructor of VectorRegion class
        VectorRegion.__init__(self, start, length, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return PixelCoordinate(self.axis1_center, self.axis2_center)

    # -----------------------------------------------------------------

    @center.setter
    def center(self, value):

        """
        This function ...
        :return:
        """

        offset = value - self.center

        self.start += offset
        self.end += offset

# -----------------------------------------------------------------

class SkyVectorRegion(VectorRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        if not isinstance(start, SkyCoordinate): raise ValueError("Start must be sky coordinate")

        # Check the length
        if not isinstance(length, Quantity): raise ValueError("Length must be an angular quantity")

        # Call the constructor of VectorRegion class
        VectorRegion.__init__(self, start, length, angle, **kwargs)

# -----------------------------------------------------------------

class PhysicalVectorRegion(VectorRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        if not isinstance(start, PhysicalCoordinate): raise ValueError("Start must be physical coordinate")

        # Check the length
        if not isinstance(length, Quantity): raise ValueError("Length must be a physical quantity of length")

        # Call the constructor of VectorRegion class
        VectorRegion.__init__(self, start, length, angle, **kwargs)

# -----------------------------------------------------------------
