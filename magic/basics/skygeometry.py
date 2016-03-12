#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.geometry Contains geometry related classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import astronomical modules
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import utils
import astropy.units as u

# Import the relevant PTS classes and modules
from .vector import Extent
from .geometry import Coordinate, Line, Circle, Ellipse, Rectangle, Polygon
from ..tools import coordinates

# -----------------------------------------------------------------

class SkyCoordinate(SkyCoord):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SkyCoordinate, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def to_pixel(self, wcs, mode='wcs'):

        """
        This function ...
        :param wcs
        :param mode:
        :return:
        """

        x, y = super(SkyCoordinate, self).to_pixel(wcs, origin=0, mode=mode)
        return Coordinate(x, y)

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, coordinate, wcs, mode='wcs'):

        """
        This function ...
        :param coordinate:
        :param wcs:
        :param mode:
        :return:
        """

        return super(SkyCoordinate).from_pixel(coordinate.x, coordinate.y, wcs, origin=0, mode=mode)

    # -----------------------------------------------------------------

    def to_astropy(self):

        """
        This function ...
        :return:
        """

        return SkyCoord(ra=self.ra.to("deg").value, dec=self.dec.to("deg").value, unit="deg", frame="fk5")

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
    def from_pixel(cls, line, wcs):

        """
        This function ...
        :param line:
        :param wcs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        start = self.start.to_pixel(wcs)
        end = self.end.to_pixel(wcs)

        # Return a new Line
        return Line(start, end)

# -----------------------------------------------------------------

class SkyEllipse(object):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=0.0):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        if isinstance(angle, float) or isinstance(angle, int): angle = Angle(0.0, "deg")

        self.center = center
        self.radius = radius
        self.angle = angle

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, ellipse, wcs):

        """
        This function ...
        :param ellipse:
        :param wcs:
        :return:
        """

        center = SkyCoordinate.from_pixel(ellipse.center, wcs)

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

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        center = self.center.to_pixel(wcs)

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
    def from_pixel(cls, circle, wcs):

        """
        This function ...
        :param circle:
        :param wcs:
        :return:
        """

        center = SkyCoordinate.from_pixel(circle.center, wcs)

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

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        center = self.center.to_pixel(wcs)

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

    def __init__(self, center, radius, angle=0.0):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        if isinstance(angle, float) or isinstance(angle, int): angle = Angle(0.0, "deg")

        self.center = center
        self.radius = radius
        self.angle = angle

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, rectangle, wcs):

        """
        This function ...
        :param rectangle:
        :param wcs:
        :return:
        """

        center = SkyCoordinate.from_pixel(rectangle.center, wcs)

        # Get the pixelscale
        radius = rectangle.radius * u.Unit("pix") * wcs.xy_average_pixelscale

        # Create a new SkyRectangle
        return cls(center, radius, rectangle.angle)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return SkyRectangle(self.center, self.radius * value, self.angle)

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

        return SkyRectangle(self.center, self.radius / value, self.angle)

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        center = self.center.to_pixel(wcs)

        pixelscale = wcs.xy_average_pixelscale
        radius = (self.radius / pixelscale).to("pix").value

        # Return a new rectangle
        return Rectangle(center, radius, self.angle)

# -----------------------------------------------------------------

class SkyPolygon(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        self.points = []

    # -----------------------------------------------------------------

    def add_point(self, point):

        """
        This function ...
        :param point:
        :return:
        """

        self.points.append(point)

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, polygon, wcs):

        """
        This function ...
        :param polygon:
        :param wcs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Create a new polygon
        polygon = Polygon()

        # Loop over the points in this SkyPolygon
        for point in self.points:

            # Convert the coordinate to a pixel position
            position = point.to_pixel(wcs)

            # Add the pixel position to the points
            polygon.add_point(position)

        # Return the new polygon
        return polygon

# -----------------------------------------------------------------
