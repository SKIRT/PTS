#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.geometry Contains geometry related classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import astronomical modules
from astropy.coordinates import SkyCoord
from astropy.wcs import utils
import astropy.units as u

# Import the relevant PTS classes and modules
from .vector import Position, Extent

# -----------------------------------------------------------------

class SkyEllipse(object):

    """
    This class ...
    """

    def __init__(self, center, radius, angle):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        self.center = center
        self.radius = radius
        self.angle = angle

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, ellipse, wcs):

        """
        This function ...
        :param ellipse:
        :param wcs:
        :return:
        """

        center = SkyCoord.from_pixel(ellipse.center.x, ellipse.center.y, wcs, mode="wcs")

        ## GET THE PIXELSCALE
        result = utils.proj_plane_pixel_scales(wcs)
        # returns: A vector (ndarray) of projection plane increments corresponding to each pixel side (axis).
        # The units of the returned results are the same as the units of cdelt, crval, and cd for the celestial WCS
        # and can be obtained by inquiring the value of cunit property of the input WCS WCS object.
        x_pixelscale = result[0] * u.Unit("deg/pix")
        y_pixelscale = result[1] * u.Unit("deg/pix")
        #pixelscale = Extent(x_pixelscale, y_pixelscale)

        major = ellipse.major * u.Unit("pix") * x_pixelscale
        minor = ellipse.minor * u.Unit("pix") * y_pixelscale

        radius = Extent(major, minor)

        # Create a new SkyEllipse
        return cls(center, radius, ellipse.angle)

    # -----------------------------------------------------------------

    @property
    def major(self):

        """
        This function ...
        :return:
        """

        return self.radius.x

    # -----------------------------------------------------------------

    @property
    def minor(self):

        """
        This function ...
        :return:
        """

        return self.radius.y

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return SkyEllipse(self.center, self.radius * value, self.angle)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center, self.radius / value, self.angle)

    # -----------------------------------------------------------------

    def to_ellipse(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        pixel_center_x, pixel_center_y = self.center.to_pixel(wcs, origin=0, mode='wcs')
        center = Position(pixel_center_x, pixel_center_y)

        ## GET THE PIXELSCALE
        result = utils.proj_plane_pixel_scales(wcs)
        # returns: A vector (ndarray) of projection plane increments corresponding to each pixel side (axis).
        # The units of the returned results are the same as the units of cdelt, crval, and cd for the celestial WCS
        # and can be obtained by inquiring the value of cunit property of the input WCS WCS object.
        x_pixelscale = result[0] * u.Unit("deg/pix")
        y_pixelscale = result[1] * u.Unit("deg/pix")
        #pixelscale = Extent(x_pixelscale, y_pixelscale)

        major = (self.major / x_pixelscale).to("pix").value
        minor = (self.minor / y_pixelscale).to("pix").value

        radius = Extent(major, minor)

        # Create a new Ellipse and return it
        return Ellipse(center, radius, self.angle)

# -----------------------------------------------------------------

class Ellipse(object):
    
    """
    This class ...
    """
    
    def __init__(self, center, radius, angle):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        self.center = center
        self.radius = radius
        self.angle = angle

    # -----------------------------------------------------------------

    @property
    def major(self):

        """
        This function ...
        :return:
        """

        return self.radius.x if isinstance(self.radius, Extent) else self.radius

    # -----------------------------------------------------------------

    @property
    def minor(self):

        """
        This function ...
        :return:
        """

        return self.radius.y if isinstance(self.radius, Extent) else self.radius

    # -----------------------------------------------------------------

    @property
    def ellipticity(self):

        """
        This function ...
        :return:
        """

        return (self.major - self.minor) / self.major

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center, self.radius * value, self.angle)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center, self.radius / value, self.angle)

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center, self.radius / value, self.angle)

    # -----------------------------------------------------------------

    def __add__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center + value, self.radius, self.angle)

    # -----------------------------------------------------------------

    def __sub__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center - value, self.radius, self.angle)

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        x_radius = self.radius.x if isinstance(self.radius, Extent) else self.radius
        y_radius = self.radius.y if isinstance(self.radius, Extent) else self.radius

        a_projected_x = x_radius * math.cos(self.angle.radian)
        b_projected_x = y_radius * math.sin(self.angle.radian)
        a_projected_y = x_radius * math.sin(self.angle.radian)
        b_projected_y = y_radius * math.cos(self.angle.radian)

        box_x_radius = max(abs(a_projected_x), abs(b_projected_x))
        box_y_radius = max(abs(a_projected_y), abs(b_projected_y))

        radius = Extent(box_x_radius, box_y_radius)

        # Return the bounding box
        return Rectangle(self.center, radius)

# -----------------------------------------------------------------

class Rectangle(object):
    
    """
    This class ...
    """
    
    def __init__(self, center, radius):
        
        """
        The constructor ...
        :param center:
        :param radius:
        :return:
        """

        self.center = center
        self.radius = radius

    # -----------------------------------------------------------------

    def contains(self, position):

        """
        This function ...
        :return:
        """

        return self.x_min <= position.x <= self.x_max and self.y_min <= position.y <= self.y_max

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This function ...
        :return:
        """

        return self.center.x - self.radius.x

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This function ...
        :return:
        """

        return self.center.x + self.radius.x

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This function ...
        :return:
        """

        return self.center.y - self.radius.y

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This function ...
        :return:
        """

        return self.center.y + self.radius.y

    # -----------------------------------------------------------------

    @property
    def lower_right(self):

        """
        This function ...
        :return:
        """

        return Position(self.x_max, self.y_min)

    # -----------------------------------------------------------------

    @property
    def upper_right(self):

        """
        This function ...
        :return:
        """

        return Position(self.x_max, self.y_max)

    # -----------------------------------------------------------------

    @property
    def upper_left(self):

        """
        This function ...
        :return:
        """

        return Position(self.x_min, self.y_max)

    # -----------------------------------------------------------------
    
    @property
    def lower_left(self):

        """
        This function ...
        :return:
        """

        return Position(self.x_min, self.y_min)

    # -----------------------------------------------------------------

    @property
    def corners(self):

        """
        This function ...
        :return:
        """

        return [self.lower_right, self.upper_right, self.upper_left, self.lower_left]

# -----------------------------------------------------------------
