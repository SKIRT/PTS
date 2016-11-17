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

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Quantity
from photutils.geometry import elliptical_overlap_grid, circular_overlap_grid, rectangular_overlap_grid

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from ..basics.mask import Mask

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

    @property
    def rotated(self):

        """
        This function ...
        :return:
        """

        return not bool(self.angle == Angle(0.0, "deg"))

    # -----------------------------------------------------------------

    @property
    def diagonal(self):

        """
        This function ...
        :return:
        """

        return self.radius.norm

    # -----------------------------------------------------------------

    @property
    def _diagonal_angle(self):

        """
        This function ...
        :return:
        """

        return np.arctan2(self.radius.axis2, self.radius.axis1)

    # -----------------------------------------------------------------

    @property
    def corners(self):

        """
        This function ...
        :return:
        """

        return [self.corner11, self.corner12, self.corner21, self.corner22]

    # -----------------------------------------------------------------

    @property
    def axis1_min(self):

        """
        This function ...
        :return:
        """

        return min([corner.axis1 for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def axis1_max(self):

        """
        This property ...
        :return:
        """

        return max([corner.axis1 for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def axis2_min(self):

        """
        This property ...
        :return:
        """

        return min([corner.axis2 for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def axis2_max(self):

        """
        This function ...
        :return:
        """

        return max([corner.axis2 for corner in self.corners])

# -----------------------------------------------------------------

class PixelRectangleRegion(RectangleRegion, PixelRegion):

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
        if not isinstance(radius, PixelStretch):
            if not isinstance(radius, float): raise ValueError("Radius must be a pixel stretch or float (for a square)")
            else: radius = PixelStretch(radius, radius)

        # Call the constructor of the base class
        super(PixelRectangleRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def corner11(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x + self.radius.x, self.center.y - self.radius.y)
        else:
            angle_to_corner = 2.0 * np.pi - self._diagonal_angle + self.angle.to("radian").value
            return PixelCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def corner12(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x + self.radius.x, self.center.y + self.radius.y)
        else:
            angle_to_corner = self._diagonal_angle + self.angle.to("radian").value
            return PixelCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def corner21(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x - self.radius.x, self.center.y + self.radius.y)
        else:
            angle_to_corner = np.pi - self._diagonal_angle + self.angle.to("radian").value
            return PixelCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def corner22(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x - self.radius.x, self.center.y - self.radius.y)
        else:
            angle_to_corner = np.pi + self._diagonal_angle + self.angle.to("radian").value
            return PixelCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def unrotated_radius(self):

        """
        This function ...
        :return:
        """

        return PixelStretch(0.5 * (self.axis1_max - self.axis1_min), 0.5 * (self.axis2_max - self.axis2_min))

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        ## SIMPLE WAY

        # data = np.zeros((y_size, x_size))

        # Convert into integers
        # x_min = int(round(self.x_min))
        # x_max = int(round(self.x_max))
        # y_min = int(round(self.y_min))
        # y_max = int(round(self.y_max))

        # data[y_min:y_max, x_min:x_max] = 1

        # Return a new Mask object
        # return Mask(data)

        ## OTHER WAY

        rel_center = self.center

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        width = 2. * self.radius.x
        height = 2. * self.radius.y

        fraction = rectangular_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, width, height, self.angle.to("radian").value, 0, 1)

        # Return the mask
        return Mask(fraction)

# -----------------------------------------------------------------

class SkyRectangleRegion(RectangleRegion, SkyRegion):

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
        if not isinstance(radius, SkyStretch):
            if not isinstance(radius, Quantity): raise ValueError("Radius must be a sky stretch or an angular quantity (for a square)")
            else: radius = SkyStretch(radius, radius)

        # Call the constructor of the base class
        super(SkyRectangleRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def corner11(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return SkyCoordinate(self.center.ra + self.radius.ra, self.center.dec - self.radius.dec)
        else:
            angle_to_corner = 2.0 * np.pi - self._diagonal_angle + self.angle.to("radian").value
            return SkyCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def corner12(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return SkyCoordinate(self.center.ra + self.radius.ra, self.center.dec + self.radius.dec)
        else:
            angle_to_corner = self._diagonal_angle + self.angle.to("radian").value
            return SkyCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def corner21(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return SkyCoordinate(self.center.ra - self.radius.ra, self.center.dec + self.radius.dec)
        else:
            angle_to_corner = np.pi - self._diagonal_angle + self.angle.to("radian").value
            return SkyCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def corner22(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return SkyCoordinate(self.center.x - self.radius.x, self.center.y - self.radius.y)
        else:
            angle_to_corner = np.pi + self._diagonal_angle + self.angle.to("radian").value
            return SkyCoordinate(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def unrotated_radius(self):

        """
        This function ...
        :return:
        """

        return SkyStretch(0.5 * (self.axis1_max - self.axis1_min), 0.5 * (self.axis2_max - self.axis2_min))

# -----------------------------------------------------------------

class PhysicalRectangleRegion(RectangleRegion, PhysicalRegion):

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
        if not isinstance(radius, PhysicalStretch):
            if not isinstance(radius, Quantity): raise ValueError("Radius must be a physical stretch or a quantity of length (for a square)")
            else: radius = PhysicalStretch(radius, radius)

        # Call the constructor of the base class
        super(PhysicalRectangleRegion, self).__init__(center, radius, angle, **kwargs)

# -----------------------------------------------------------------
