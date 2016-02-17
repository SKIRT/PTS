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
import numpy as np

# Function to create mask from ellipse
from photutils.geometry import elliptical_overlap_grid, circular_overlap_grid, rectangular_overlap_grid

# Import the relevant PTS classes and modules
from .vector import Position, Extent
from .mask import Mask

# -----------------------------------------------------------------

class Line(object):

    """
    This class ...
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

        return Extent(self.end.x - self.start.x, self.end.y - self.start.y).norm

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

class Circle(object):

    """
    This class ...
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
        :param position:
        :return:
        """

        return self.x_min <= position.x <= self.x_max and self.y_min <= position.y <= self.y_max

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

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        ## NAIVE WAY

        data = np.zeros((y_size, x_size))

        # Convert into integers
        x_min = int(round(self.x_min))
        x_max = int(round(self.x_max))
        y_min = int(round(self.y_min))
        y_max = int(round(self.y_max))

        data[y_min:y_max, x_min:x_max] = 1

        # Return a new Mask object
        return Mask(data)

        ## OTHER WAY

        #rel_center = rectangle.center

        #x_min = - rel_center.x
        #x_max = x_size - rel_center.x
        #y_min = - rel_center.y
        #y_max = y_size - rel_center.y

        #width = 2. * rectangle.radius.x
        #height = 2. * rectangle.radius.y

        #fraction = rectangular_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, width, height, use_exact=0, subpixels=1)

        #return cls(fraction)

# -----------------------------------------------------------------
