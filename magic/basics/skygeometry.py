#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.geometry Contains geometry related classes.

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
from .geometry import Line, Circle, Ellipse, Rectangle
from ..tools import coordinates

# -----------------------------------------------------------------

class SkyLine(object):

    """
    This function ...
    """

    def __init__(self, start, end):

        """
        This function ...
        :param start:
        :param end:
        :return:
        """

        self.start = start
        self.end = end

    # -----------------------------------------------------------------

    @property
    def length(self):

        """
        This function ...
        :return:
        """

        dec_start = self.start.dec.to("deg").value
        dec_end = self.end.dec.to("deg").value

        dec_center = 0.5 * (dec_start + dec_end)

        ra_start = self.start.ra.to("deg").value
        ra_end = self.end.ra.to("deg").value

        # Calculate the actual RA and DEC distance in degrees
        ra_distance = abs(coordinates.ra_distance(dec_center, ra_start, ra_end))
        dec_distance = abs(dec_end - dec_start)

        ra_span = ra_distance * u.Unit("deg")
        dec_span = dec_distance * u.Unit("deg")

        return math.sqrt(ra_span**2 + dec_span**2)

    # -----------------------------------------------------------------

    @classmethod
    def from_line(cls, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """



    # -----------------------------------------------------------------

    def to_line(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        start_x, start_y = self.start.to_pixel(wcs, mode="wcs", origin=0)
        end_x, end_y = self.end.to_pixel(wcs, mode="wcs", origin=0)

        # Return a new Line
        return Line(Position(start_x, start_y), Position(end_x, end_y))

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

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return SkyEllipse(self.center, self.radius / value, self.angle)

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

class SkyCircle(object):

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

    @classmethod
    def from_circle(cls, circle, wcs):

        """
        This function ...
        :param circle:
        :param wcs:
        :return:
        """

        center = SkyCoord.from_pixel(circle.center.x, circle.center.y, wcs, mode="wcs")

        # Get the pixelscale
        radius = circle.radius * u.Unit("pix") * wcs.xy_average_pixelscale

        # Create a new SkyCircle
        return cls(center, radius)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return SkyCircle(self.center, self.radius * value)

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return SkyCircle(self.center, self.radius / value)

    # -----------------------------------------------------------------

    def to_circle(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        pixel_center_x, pixel_center_y = self.center.to_pixel(wcs, origin=0, mode='wcs')
        center = Position(pixel_center_x, pixel_center_y)

        ## GET THE PIXELSCALE
        pixelscale = wcs.xy_average_pixelscale
        radius = (self.radius / pixelscale).to("pix").value

        # Create a new Circle and return it
        return Circle(center, radius)

# -----------------------------------------------------------------

class SkyRectangle(object):

    """
    This class
    """

    def __init__(self, center, radius):

        """
        This function ...
        :param center:
        :param radius:
        :return:
        """

        self.center = center
        self.radius = radius

    # -----------------------------------------------------------------

    @classmethod
    def from_rectangle(cls, rectangle, wcs):

        """
        This function ...
        :param rectangle:
        :param wcs:
        :return:
        """

        center = SkyCoord.from_pixel(rectangle.center.x, rectangle.center.y, wcs, mode="wcs")

        # Get the pixelscale
        radius = rectangle.radius * u.Unit("pix") * wcs.xy_average_pixelscale

        # Create a new SkyRectangle
        return cls(center, radius)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return SkyRectangle(self.center, self.radius * value)

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return SkyRectangle(self.center, self.radius / value)

    # -----------------------------------------------------------------

    def to_rectangle(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        pixel_center_x, pixel_center_y = self.center.to_pixel(wcs, origin=0, mode='wcs')
        center = Position(pixel_center_x, pixel_center_y)

        pixelscale = wcs.xy_average_pixelscale
        radius = (self.radius / pixelscale).to("pix").value

        # Return a new rectangle
        return Rectangle(center, radius)

# -----------------------------------------------------------------
