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
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .vector import Extent
from .geometry import Coordinate, Line, Circle, Ellipse, Rectangle, Polygon
from ..tools import coordinates

# -----------------------------------------------------------------

class SkyComposite(object):

    """
    This class ...
    """

    def __init__(self, base, exclude, meta=None):

        """
        This function ...
        :param base:
        :param exclude:
        :param meta:
        """

        self.base = base
        self.exclude = exclude

        self.meta = meta if meta is not None else dict()

    # -----------------------------------------------------------------

    @property
    def x(self):
        """
        This function ...
        :return:
        """

        return self.base.x

    # -----------------------------------------------------------------

    @property
    def y(self):
        """
        This function ...
        :return:
        """

        return self.base.y

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :param coordinate_system:
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta: suffix += " " + key + " = " + str(self.meta[key])
        else:
            suffix = ""

        pass

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :return:
        """

        return SkyComposite(self.base + extent, self.exclude + extent, self.meta)

    # -----------------------------------------------------------------

    def __iadd__(self, extent):

        """
        This function ...
        :return:
        """

        self.base += extent
        self.exclude += extent

        return self

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :return:
        """

        return SkyComposite(self.base - extent, self.exclude - extent, self.meta)

    # -----------------------------------------------------------------

    def __isub__(self, extent):

        """
        This function ...
        :return:
        """

        self.base -= extent
        self.exclude -= extent

        return self

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :return:
        """

        return SkyComposite(self.base * value, self.exclude * value, self.meta)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :return:
        """

        self.base *= value
        self.exclude *= value

        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :return:
        """

        return SkyComposite(self.base / value, self.exclude / value, self.meta)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :return:
        """

        self.base /= value
        self.exclude /= value

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

        self.__idiv__(value)
        return self

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

        meta = kwargs.pop("meta", None)

        # Call the constructor of the base class
        super(SkyCoordinate, self).__init__(*args, **kwargs)

        # Set meta information
        if meta is not None: self.meta = meta
        else: self.meta = dict()

    # -----------------------------------------------------------------

    def to_pixel(self, wcs, mode='wcs'):

        """
        This function ...
        :param wcs
        :param mode:
        :return:
        """

        x, y = super(SkyCoordinate, self).to_pixel(wcs, origin=0, mode=mode)
        return Coordinate(x, y, meta=self.meta)

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

        skycoordinate = super(SkyCoordinate, cls).from_pixel(coordinate.x, coordinate.y, wcs, origin=0, mode=mode)
        return cls(ra=skycoordinate.ra.deg, dec=skycoordinate.dec.deg, unit="deg", meta=coordinate.meta)

    # -----------------------------------------------------------------

    def to_astropy(self):

        """
        This function ...
        :return:
        """

        return SkyCoord(ra=self.ra.to("deg").value, dec=self.dec.to("deg").value, unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :param coordinate_system:
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Get the RA and DEC
        ra_deg = self.ra.to("deg").value
        dec_deg = self.dec.to("deg").value

        # Create and return the line
        if coordinate_system: line = "fk5;point({},{})".format(ra_deg, dec_deg) + suffix
        else: line = "point({},{})".format(ra_deg, dec_deg) + suffix
        return line

# -----------------------------------------------------------------

class SkyLine(object):

    """
    This function ...
    """

    def __init__(self, start, end, meta=None):

        """
        This function ...
        :param start:
        :param end:
        :return:
        """

        self.start = start
        self.end = end

        # Set meta information
        self.meta = meta if meta is not None else dict()

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

        ra_span = ra_distance * Unit("deg")
        dec_span = dec_distance * Unit("deg")

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

        start = line.start.to_sky(wcs)
        end = line.end.to_sky(wcs)
        return cls(start, end, meta=line.meta)

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
        return Line(start, end, meta=self.meta)

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Get the RA and DEC of the 'start' coordinate
        start_ra = self.start.ra.to("deg").value
        start_dec = self.start.dec.to("deg").value

        # Get the RA and DEC of the 'end' coordinate
        end_ra = self.end.ra.to("deg").value
        end_dec = self.end.dec.to("deg").value

        # Create and return the line
        if coordinate_system: line = "fk5;line({},{},{},{})".format(start_ra, start_dec, end_ra, end_dec) + suffix
        else: line = "line({},{},{},{})".format(start_ra, start_dec, end_ra, end_dec) + suffix
        return line

# -----------------------------------------------------------------

class SkyEllipse(object):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=0.0, meta=None):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        if isinstance(angle, float) or isinstance(angle, int): angle = Angle(angle, "deg")

        self.center = center
        self.radius = radius
        self.angle = angle

        # Set meta information
        self.meta = meta if meta is not None else dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, ellipse, wcs):

        """
        This function ...
        :param ellipse:
        :param wcs:
        :return:
        """

        center = SkyCoordinate.from_pixel(Coordinate(ellipse.center.x, ellipse.center.y), wcs)

        ## GET THE PIXELSCALE
        result = utils.proj_plane_pixel_scales(wcs)
        # returns: A vector (ndarray) of projection plane increments corresponding to each pixel side (axis).
        # The units of the returned results are the same as the units of cdelt, crval, and cd for the celestial WCS
        # and can be obtained by inquiring the value of cunit property of the input WCS WCS object.
        x_pixelscale = result[0] * Unit("deg/pix")
        y_pixelscale = result[1] * Unit("deg/pix")
        #pixelscale = Extent(x_pixelscale, y_pixelscale)

        major = ellipse.major * Unit("pix") * x_pixelscale
        minor = ellipse.minor * Unit("pix") * y_pixelscale

        radius = Extent(major, minor)

        # Create a new SkyEllipse
        return cls(center, radius, ellipse.angle, meta=ellipse.meta)

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
        return SkyRectangle(self.center, radius)

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
        x_pixelscale = result[0] * Unit("deg/pix")
        y_pixelscale = result[1] * Unit("deg/pix")
        #pixelscale = Extent(x_pixelscale, y_pixelscale)

        major = (self.major / x_pixelscale).to("pix").value
        minor = (self.minor / y_pixelscale).to("pix").value

        radius = Extent(major, minor)

        # Create a new Ellipse and return it
        return Ellipse(center, radius, self.angle, meta=self.meta)

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Get ellipse properties
        ra_deg = self.center.ra.to("deg").value
        dec_deg = self.center.dec.to("deg").value
        major = self.major.to("arcsec").value
        minor = self.minor.to("arcsec").value
        angle = self.angle.degree

        # Create and return the line
        if coordinate_system: line = "fk5;ellipse(%s,%s,%.2f\",%.2f\",%s)" % (ra_deg, dec_deg, major, minor, angle)
        else: line = "ellipse(%s,%s,%.2f\",%.2f\",%s)" % (ra_deg, dec_deg, major, minor, angle)
        line += suffix
        return line

# -----------------------------------------------------------------

class SkyCircle(object):

    """
    This class ...
    """

    def __init__(self, center, radius, meta=None):

        """
        The constructor ...
        :param center:
        :param radius:
        :return:
        """

        self.center = center
        self.radius = radius

        # Set meta information
        self.meta = meta if meta is not None else dict()

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
        radius = circle.radius * Unit("pix") * wcs.xy_average_pixelscale

        # Create a new SkyCircle
        return cls(center, radius, meta=circle.meta)

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
        return Circle(center, radius, meta=self.meta)

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Get circle properties
        ra_deg = self.center.ra.to("deg").value
        dec_deg = self.center.dec.to("deg").value
        radius = self.radius.to("arcsec").value

        # Create and return the line
        if coordinate_system: line = "fk5;circle(%s,%s,%.2f\")" % (ra_deg, dec_deg, radius)
        else:  line = "circle(%s,%s,%.2f\")" % (ra_deg, dec_deg, radius)
        line += suffix
        return line

# -----------------------------------------------------------------

class SkyRectangle(object):

    """
    This class
    """

    def __init__(self, center, radius, angle=0.0, meta=None):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        if isinstance(angle, float) or isinstance(angle, int): angle = Angle(angle, "deg")

        self.center = center
        self.radius = radius
        self.angle = angle

        # Set meta information
        self.meta = meta if meta is not None else dict()

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
        radius = rectangle.radius * Unit("pix") * wcs.xy_average_pixelscale

        # Create a new SkyRectangle
        return cls(center, radius, rectangle.angle, meta=rectangle.meta)

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
        radius_x = (self.radius.x / pixelscale).to("pix").value
        radius_y = (self.radius.y / pixelscale).to("pix").value

        # Return a new rectangle
        return Rectangle(center, Extent(radius_x, radius_y), self.angle, meta=self.meta)

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Get rectangle properties
        center_ra = self.center.ra.to("deg").value
        center_dec = self.center.dec.to("deg").value
        width = 2.0 * self.radius.to("arcsec").value
        height = 2.0 * self.radius.y.to("arcsec").value
        angle = self.angle.degree

        # Create and return the line
        if coordinate_system: line = "fk5;box(%s,%s,%.2f\",%.2f\",%s)" % (center_ra, center_dec, width, height, angle)
        else: line = "box(%s,%s,%.2f\",%.2f\",%s)" % (center_ra, center_dec, width, height, angle)
        line += suffix
        return line

# -----------------------------------------------------------------

class SkyPolygon(object):

    """
    This class ...
    """

    def __init__(self, meta=None):

        """
        This function ...
        :return:
        """

        self.points = []

        # Set meta information
        self.meta = meta if meta is not None else dict()

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

        # Create the SkyPolygon object
        skypolygon = cls(meta=polygon.meta)

        # Loop over the points in the polygon
        for point in polygon.points:

            # Convert the coordinate to a SkyCoordinate
            coordinate = point.to_sky(wcs)

            # Add the SkyCoordinate to the points of the new SkyPolygon
            skypolygon.add_point(coordinate)

        # Return the new sky polygon
        return skypolygon

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Create a new polygon
        polygon = Polygon(meta=self.meta)

        # Loop over the points in this SkyPolygon
        for point in self.points:

            # Convert the coordinate to a pixel position
            position = point.to_pixel(wcs)

            # Add the pixel position to the points
            polygon.add_point(position)

        # Return the new polygon
        return polygon

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Initialize line
        if coordinate_system: line = "fk5;polygon("
        else: line = "polygon("

        lines = []

        # Add the points to the line
        for point in self.points:

            ra = point.ra.to("deg").value
            dec = point.dec.to("deg").value
            lines.append("{},{}".format(ra, dec))

        # Finish line
        line += ",".join(lines) + ")" + suffix

        # Return the line
        return line

# -----------------------------------------------------------------
