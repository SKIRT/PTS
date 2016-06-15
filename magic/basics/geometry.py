#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.geometry Contains geometry related classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np

# Import astronomical modules
from astropy.coordinates import Angle

# Function to create mask from ellipse
from photutils.geometry import elliptical_overlap_grid, circular_overlap_grid, rectangular_overlap_grid
from PIL import Image, ImageDraw

# Import the relevant PTS classes and modules
from .vector import Position, Extent
from .mask import Mask

# -----------------------------------------------------------------

class Composite(object):

    """
    This function ...
    """

    def __init__(self, base, exclude, meta=None):

        """
        This function ...
        :return:
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

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        base = self.base.to_mask(x_size, y_size)
        exclude = self.exclude.to_mask(x_size, y_size)

        # Return the mask
        return base * np.logical_not(exclude)

    # -----------------------------------------------------------------

    def to_region_string(self):

        """
        This function ...
        :return:
        """

        line1 = self.base.to_region_string()
        line2 = self.exclude.to_region_string()

        return line1 + "\n" + line2.replace("image;", "image;-")

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :return:
        """

        return Composite(self.base + extent, self.exclude + extent, self.meta)

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

        return Composite(self.base - extent, self.exclude - extent, self.meta)

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

        return Composite(self.base * value, self.exclude * value, self.meta)

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

        return Composite(self.base / value, self.exclude / value, self.meta)

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

class Coordinate(Position):

    """
    This class ...
    """

    def __init__(self, x, y, meta=None):

        """
        The constructor ...
        :param x:
        :param y:
        :param meta:
        :return:
        """

        # Call the constructor of the base class
        super(Coordinate, self).__init__(x, y)

        # Set the meta information
        self.meta = meta if meta is not None else dict()

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        return Mask.empty(x_size, y_size)

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module scope
        from .skygeometry import SkyCoordinate

        # Return a new SkyCoordinate
        return SkyCoordinate.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    def to_region_string(self):

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

        # Create and return the line
        line = "image;point({},{})".format(self.x+1, self.y+1) + suffix
        return line

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        return Rectangle(Position(self.x, self.y), Extent(1.0,1.0))

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Coordinate(self.x + extent.x, self.y + extent.y, self.meta)

    # -----------------------------------------------------------------

    def __iadd__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        self.x += extent.x
        self.y += extent.y

        return self

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Coordinate(self.x - extent.x, self.y - extent.y, self.meta)

    # -----------------------------------------------------------------

    def __isub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        self.x -= extent.x
        self.y -= extent.y

        return self

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Coordinate(self.x, self.y, self.meta)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Coordinate(self.x, self.y, self.meta)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__idiv__(value)

# -----------------------------------------------------------------

class Line(object):

    """
    This class ...
    """

    def __init__(self, start, end, meta=None):

        """
        This function ...
        :param start:
        :param end:
        :return:
        """

        # The start and end coordinate
        self.start = start
        self.end = end

        # Set meta information
        self.meta = meta if meta is not None else dict()

    # -----------------------------------------------------------------

    @property
    def x(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.start.x + self.end.x)

    # -----------------------------------------------------------------

    @property
    def y(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.start.y + self.end.y)

    # -----------------------------------------------------------------

    @property
    def length(self):

        """
        This function ...
        :return:
        """

        return Extent(self.end.x - self.start.x, self.end.y - self.start.y).norm

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This property ...
        :return:
        """

        return min(self.start.x, self.end.x)

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This property ...
        :return:
        """

        return max(self.start.x, self.end.x)

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This property ...
        :return:
        """

        return min(self.start.y, self.end.y)

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This property ...
        :return:
        """

        return max(self.start.x, self.start.y)

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        return Mask.empty(x_size, y_size)

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module scope
        from .skygeometry import SkyLine

        # Return a new SkyLine
        return SkyLine.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    def to_region_string(self):

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

        # Get the line properties
        start = self.start
        end = self.end

        # Create and return the line
        line = "image;line({},{},{},{})".format(start.x+1, start.y+1, end.x+1, end.y+1) + suffix
        return line

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        x_min = self.x_min
        x_max = self.x_max
        y_min = self.y_min
        y_max = self.y_max

        center = Position(0.5*(x_min+x_max), 0.5*(y_min+y_max))
        radius = Extent(0.5*(x_max-x_min), 0.5*(y_max-y_min))

        return Rectangle(center, radius)

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :return:
        """

        return Line(self.start + extent, self.end + extent, self.meta)

    # -----------------------------------------------------------------

    def __iadd__(self, extent):

        """
        This function ...
        :return:
        """

        self.start += extent
        self.end += extent

        return self

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Line(self.start - extent, self.end - extent, self.meta)

    # -----------------------------------------------------------------

    def __isub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        self.start -= extent
        self.end -= extent

        return self

# -----------------------------------------------------------------

class Circle(object):

    """
    This class ...
    """

    def __init__(self, center, radius, meta=None):

        """
        This function ...
        :param center:
        :param radius:
        :return:
        """

        # Set basic properties
        self.center = center
        self.radius = radius

        # Set the meta information
        self.meta = meta if meta is not None else dict()

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

        return math.pi * self.radius**2

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :return:
        """

        return Circle(self.center, self.radius * value, self.meta)

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

        return Circle(self.center, self.radius / value, self.meta)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :return:
        """

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

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        rel_center = self.center

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        fraction = circular_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, self.radius, use_exact=0, subpixels=1)

        # Return a new mask
        return Mask(fraction)

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module level
        from .skygeometry import SkyCircle

        # Return a new SkyCircle
        return SkyCircle.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    def to_region_string(self):

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
        center = self.center
        radius = self.radius

        # Create and return the line
        line = "image;circle({},{},{})".format(center.x+1, center.y+1, radius) + suffix
        return line

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Return the bounding box
        return Rectangle(self.center, Extent(self.radius, self.radius))

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Circle(self.center + extent, self.radius, self.meta)

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Circle(self.center - extent, self.radius, self.meta)

# -----------------------------------------------------------------

class Ellipse(object):
    
    """
    This class ...
    """
    
    def __init__(self, center, radius, angle=0.0, meta=None):

        """
        The constructor ...
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

    @property
    def ellipticity(self):

        """
        This function ...
        :return:
        """

        return (self.major - self.minor) / self.major

    # -----------------------------------------------------------------

    @property
    def circumference(self):

        """
        This function ...
        :return:
        """

        h = (self.radius.x - self.radius.y)**2 / (self.radius.x + self.radius.y)**2

        # Approximation
        return math.pi * (self.radius.x + self.radius.y) * (1. + 3. * h / (10. + math.sqrt(4.-3.*h)))

    # -----------------------------------------------------------------

    @property
    def area(self):

        """
        This function ...
        :return:
        """

        return math.pi * self.radius.x * self.radius.y

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center, self.radius * value, self.angle)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.radius *= value

        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center, self.radius / value, self.angle)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.radius /= value

        return self

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center, self.radius / value, self.angle)

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.radius /= value

        return self

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Ellipse(self.center + extent, self.radius, self.angle, self.meta)

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param value:
        :return:
        """

        return Ellipse(self.center - extent, self.radius, self.angle, self.meta)

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

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        rel_center = self.center

        a = self.radius.x if isinstance(self.radius, Extent) else self.radius
        b = self.radius.y if isinstance(self.radius, Extent) else self.radius

        # theta in radians !
        theta = self.angle.radian

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        # Calculate the mask
        fraction = elliptical_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, a, b, theta, use_exact=0, subpixels=1)

        #xmin, xmax, ymin, ymax : float
        #    Extent of the grid in the x and y direction.
        #nx, ny : int
        #    Grid dimensions.
        #rx : float
        #    The semimajor axis of the ellipse.
        #ry : float
        #    The semiminor axis of the ellipse.
        #theta : float
        #    The position angle of the semimajor axis in radians (counterclockwise).
        #use_exact : 0 or 1
        #    If set to 1, calculates the exact overlap, while if set to 0, uses a
        #    subpixel sampling method with ``subpixel`` subpixels in each direction.
        #subpixels : int
        #    If ``use_exact`` is 0, each pixel is resampled by this factor in each
        #    dimension. Thus, each pixel is divided into ``subpixels ** 2``
        #    subpixels.

        return Mask(fraction)

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module level
        from .skygeometry import SkyEllipse

        # Return a new SkyEllipse
        return SkyEllipse.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    def to_region_string(self):

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
        center = self.center
        major = self.major
        minor = self.minor
        angle = self.angle.degree

        # Create and return the line
        line = "image;ellipse({},{},{},{},{})".format(center.x+1, center.y+1, major, minor, angle) + suffix
        return line

# -----------------------------------------------------------------

class Rectangle(object):
    
    """
    This class ...
    """
    
    def __init__(self, center, radius, angle=0.0, meta=None):
        
        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        if isinstance(angle, float) or isinstance(angle, int): angle = Angle(angle, "deg")

        self.center = center
        self.radius = radius
        self.angle = angle

        self.meta = meta if meta is not None else dict()

    # -----------------------------------------------------------------

    @property
    def rotated(self):

        """
        This function ...
        :return:
        """

        return not bool(self.angle == Angle(0.0,"deg"))

    # -----------------------------------------------------------------

    def contains(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        if not self.rotated: return self.x_min <= position.x <= self.x_max and self.y_min <= position.y <= self.y_max
        else:

            # NOTE: does not work properly yet!

            # Source: http://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle
            # AM = position - self.upper_left
            # AB = self.upper_right - self.upper_left
            # AD = self.lower_left - self.upper_left
            am = position - self.upper_left
            ab = self.upper_right - self.upper_left
            ad = self.lower_left - self.upper_left

            return (0 < am.dot(ab) < ab.dot(ab)) and (0 < am.dot(ad) < ad.dot(ad))

    # -----------------------------------------------------------------

    @property
    def circumference(self):

        """
        This function ...
        :return:
        """

        return 4. * (self.radius.x + self.radius.y)

    # -----------------------------------------------------------------

    @property
    def area(self):

        """
        This function ...
        :return:
        """

        return 4. * self.radius.x * self.radius.y

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This function ...
        :return:
        """

        return min([corner.x for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This function ...
        :return:
        """

        return max([corner.x for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This function ...
        :return:
        """

        return min([corner.y for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This function ...
        :return:
        """

        return max([corner.y for corner in self.corners])

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

        return np.arctan2(self.radius.y, self.radius.x)

    # -----------------------------------------------------------------

    @property
    def lower_right(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return Position(self.center.x + self.radius.x, self.center.y - self.radius.y)
        else:

            angle_to_corner = 2.0 * np.pi - self._diagonal_angle + self.angle.to("radian").value
            return Position(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def upper_right(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return Position(self.center.x + self.radius.x, self.center.y + self.radius.y)
        else:

            angle_to_corner = self._diagonal_angle + self.angle.to("radian").value
            return Position(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def upper_left(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return Position(self.center.x - self.radius.x, self.center.y + self.radius.y)
        else:

            angle_to_corner = np.pi - self._diagonal_angle + self.angle.to("radian").value
            return Position(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------
    
    @property
    def lower_left(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return Position(self.center.x - self.radius.x, self.center.y - self.radius.y)
        else:

            angle_to_corner = np.pi + self._diagonal_angle + self.angle.to("radian").value
            return Position(self.diagonal * np.cos(angle_to_corner), self.diagonal * np.sin(angle_to_corner))

    # -----------------------------------------------------------------

    @property
    def corners(self):

        """
        This function ...
        :return:
        """

        return [self.lower_right, self.upper_right, self.upper_left, self.lower_left]

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        ## SIMPLE WAY

        #data = np.zeros((y_size, x_size))

        # Convert into integers
        #x_min = int(round(self.x_min))
        #x_max = int(round(self.x_max))
        #y_min = int(round(self.y_min))
        #y_max = int(round(self.y_max))

        #data[y_min:y_max, x_min:x_max] = 1

        # Return a new Mask object
        #return Mask(data)

        ## OTHER WAY

        rel_center = self.center

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        width = 2. * self.radius.x
        height = 2. * self.radius.y

        fraction = rectangular_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, width, height, self.angle.to("radian").value, 0, 1)
        return Mask(fraction)

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module level
        from .skygeometry import SkyRectangle

        # Return a new SkyRectangle
        return SkyRectangle.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    def to_region_string(self):

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
        center = self.center
        radius = self.radius
        width = 2.0 * radius.x
        height = 2.0 * radius.y
        angle = self.angle.degree

        # Create and return the line
        line = "image;box({},{},{},{},{})".format(center.x+1, center.y+1, width, height, angle) + suffix
        return line

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        return Rectangle(self.center, self.radius)

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Rectangle(self.center + extent, self.radius, self.angle, self.meta)

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Rectangle(self.center - extent, self.radius, self.angle, self.meta)

# -----------------------------------------------------------------

class Polygon(object):

    """
    This function ...
    """

    def __init__(self, meta=None):

        """
        The constructor ...
        :return:
        """

        # Initilize a list to contain the corner points
        self.points = []

        # Set the meta information
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

    @property
    def circumference(self):

        """
        This function ...
        :return:
        """

        total = 0.0

        # Calculate the circumference
        for i in range(1, len(self.points)): total += (self.points[i]-self.points[i-1]).norm
        total += (self.points[0] - self.points[len(self.points)-1]).norm

        # Return the result
        return total

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This property ...
        :return:
        """

        return min([point.x for point in self.points])

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This property ...
        :return:
        """

        return max([point.x for point in self.points])

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This function ...
        :return:
        """

        return min([point.y for point in self.points])

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This function ...
        :return:
        """

        return max([point.y for point in self.points])

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        width = x_size
        height = y_size

        #points = [(x1,y1),(x2,y2),...] or [x1,y1,x2,y2,...]
        points = []
        for point in self.points:

            #points.append((point.x, point.y))
            points.append(float(point.x))
            points.append(float(point.y))

        img = Image.new('L', (width, height), 0)
        ImageDraw.Draw(img).polygon(points, outline=1, fill=1)

        return Mask(np.array(img))

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module level
        from .skygeometry import SkyPolygon

        # Return a new SkyPolygon
        return SkyPolygon.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    def to_region_string(self):

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
        line = "image;polygon("

        # Add the points to the line
        for point in self.points: line += "{},{}".format(point.x+1, point.y+1)

        # Finish line
        line += ")" + suffix

        # Return the line
        return line

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        x_min = self.x_min
        x_max = self.x_max
        y_min = self.y_min
        y_max = self.y_max

        center = Position(0.5*(x_min+x_max), 0.5*(y_min+y_max))
        radius = Extent(0.5*(x_max - x_min), 0.5*(y_max - y_min))

        return Rectangle(center, radius)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This property returns the center coordinate of the polygon, calcualted as the mean x value with the mean y value
        :return:
        """

        x = np.array([point.x for point in self.points])
        y = np.array([point.y for point in self.points])

        return Coordinate(np.mean(x), np.mean(y))

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :param extent
        :return:
        """

        # Create a new polygon
        polygon = Polygon(self.meta)

        for point in self.points: polygon.add_point(point + extent)

        return polygon

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        # Create a new polygon
        polygon = Polygon(self.meta)

        for point in self.points: polygon.add_point(point - extent)

        return polygon

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new = Polygon(meta=self.meta)

        center = self.center

        for point in self.points:

            point_to_center = point - self.center # is Extent with x and y
            new_point = Coordinate(center.x + point_to_center.x * value, center.y + point_to_center.y * value)
            new.add_point(new_point)

        # Return the new polygon
        return new

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        center = self.center

        for point in self.points:

            point_to_center = point - self.center  # is Extent with x and y
            point.x = center.x + point_to_center.x * value
            point.y = center.y + point_to_center.y * value

        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__mul__(1./value)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__imul__(1./value)

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__idiv__(value)

# -----------------------------------------------------------------
