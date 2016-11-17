#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.point Contains the PointRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Unit
from photutils.geometry import elliptical_overlap_grid, circular_overlap_grid, rectangular_overlap_grid

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .rectangle import PixelRectangleRegion, SkyRectangleRegion, PhysicalRectangleRegion
from ..basics.mask import Mask

# -----------------------------------------------------------------

class EllipseRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the angle
        if angle is None: angle = Angle(0., "deg")
        elif not isinstance(angle, Angle): raise ValueError("Angle must be an Astropy Angle object")

        # Set attributes
        self.center = center
        self.radius = radius
        self.angle = angle

        # Call the constructor of the base class
        super(EllipseRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelEllipseRegion(EllipseRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be a pixel coordinate")
        if not isinstance(radius, PixelStretch): raise ValueError("Radius must be a pixel stretch")

        # Call the constructor of the base class
        super(PixelEllipseRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def major(self):

        """
        This function ...
        :return:
        """

        return self.radius.x

    # -----------------------------------------------------------------

    @major.setter
    def major(self, value):

        """
        This function ...
        :return:
        """

        self.radius.x = value

    # -----------------------------------------------------------------

    @property
    def minor(self):

        """
        This function ...
        :return:
        """

        return self.radius.y

    # -----------------------------------------------------------------

    @minor.setter
    def minor(self, value):

        """
        This function ...
        :return:
        """

        self.radius.y = value

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        rel_center = self.center

        a = self.radius.x
        b = self.radius.y

        # theta in radians !
        theta = self.angle.radian

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        # Calculate the mask
        fraction = elliptical_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, a, b, theta, use_exact=0, subpixels=1)

        # xmin, xmax, ymin, ymax : float
        #    Extent of the grid in the x and y direction.
        # nx, ny : int
        #    Grid dimensions.
        # rx : float
        #    The semimajor axis of the ellipse.
        # ry : float
        #    The semiminor axis of the ellipse.
        # theta : float
        #    The position angle of the semimajor axis in radians (counterclockwise).
        # use_exact : 0 or 1
        #    If set to 1, calculates the exact overlap, while if set to 0, uses a
        #    subpixel sampling method with ``subpixel`` subpixels in each direction.
        # subpixels : int
        #    If ``use_exact`` is 0, each pixel is resampled by this factor in each
        #    dimension. Thus, each pixel is divided into ``subpixels ** 2``
        #    subpixels.

        return Mask(fraction)

# -----------------------------------------------------------------

class SkyEllipseRegion(EllipseRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be a sky coordinate")
        if not isinstance(radius, SkyStretch): raise ValueError("Radius must be a sky stretch")

        # Call the constructor of the base class
        super(SkyEllipseRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def major(self):

        """
        This function ...
        :return:
        """

        return self.radius.ra

    # -----------------------------------------------------------------

    @major.setter
    def major(self, value):

        """
        This function ...
        :return:
        """

        self.radius.ra = value

    # -----------------------------------------------------------------

    @property
    def minor(self):

        """
        This function ...
        :return:
        """

        return self.radius.dec

    # -----------------------------------------------------------------

    @minor.setter
    def minor(self, value):

        """
        This function ...
        :return:
        """

        self.radius.dec = value

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        x_radius = self.radius.ra.to("arcsec").value
        y_radius = self.radius.dec.to("arcsec").value

        a_projected_x = x_radius * math.cos(self.angle.radian)
        b_projected_x = y_radius * math.sin(self.angle.radian)
        a_projected_y = x_radius * math.sin(self.angle.radian)
        b_projected_y = y_radius * math.cos(self.angle.radian)

        box_x_radius = max(abs(a_projected_x), abs(b_projected_x)) * Unit("arcsec")
        box_y_radius = max(abs(a_projected_y), abs(b_projected_y)) * Unit("arcsec")

        radius = SkyStretch(box_x_radius, box_y_radius)

        # Return the bounding box
        return SkyRectangleRegion(self.center, radius)

# -----------------------------------------------------------------

class PhysicalEllipseRegion(EllipseRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be a physical coordinate")
        if not isinstance(radius, PhysicalStretch): raise ValueError("Radius must be a physical stretch")

        # Call the constructor of the base class
        super(PhysicalEllipseRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def major(self):

        """
        This function ...
        :return:
        """

        return self.radius.length1

    # -----------------------------------------------------------------

    @major.setter
    def major(self, value):

        """
        This function ...
        :return:
        """

        self.radius.length1 = value

    # -----------------------------------------------------------------

    @property
    def minor(self):

        """
        This function ...
        :return:
        """

        return self.radius.length2

    # -----------------------------------------------------------------

    @minor.setter
    def minor(self, value):

        """
        This function ...
        :return:
        """

        self.radius.length2 = value

# -----------------------------------------------------------------
