#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.vector Contains the Vector class and the subclasses Position and Extent.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import the relevant PTS classes and modules
from ...core.tools import types

# -----------------------------------------------------------------

class Vector(object):

    """
    This class ...
    """

    def __init__(self, x, y):

        """
        The constructor ...
        """

        # Set the attributes
        self.x = x
        self.y = y

    # -----------------------------------------------------------------

    @property
    def norm(self):
        return math.sqrt(self.x**2 + self.y**2)

    # -----------------------------------------------------------------

    @property
    def angle(self):
        return math.atan2(self.x, self.y)

    # -----------------------------------------------------------------

    @property
    def polar(self):
        return self.norm, self.angle

    # -----------------------------------------------------------------

    @property
    def cartesian(self):
        return self.x, self.y

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        """

        return self.__class__.__name__ + '(x={0}, y={1})'.format(self.x, self.y)

     # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function
        """

        return '<' + self.__class__.__name__ + ' x={0}, y={1}>'.format(self.x, self.y)

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        from .mask import Mask
        return Mask.empty(x_size, y_size)

    # -----------------------------------------------------------------

    def dot(self, vector):

        """
        This function calculates the dot product with another vector
        :param vector:
        :return:
        """

        return self.x * vector.x + self.y * vector.y

# -----------------------------------------------------------------

class Vector3D(object):

    """
    This class ...
    """

    def __init__(self, x, y, z):

        """
        The constructor ...
        :param x:
        :param y:
        :param z:
        """

        # Set the attributes
        self.x = x
        self.y = y
        self.z = z

    # -----------------------------------------------------------------

    @property
    def norm(self):

        """
        This function ...
        """

        return math.sqrt(self.x ** 2 + self.y ** 2)

# -----------------------------------------------------------------

class Position(Vector):

    """
    This class ...
    """

    def __init__(self, x, y):

        """
        The constructor ...
        """

        # Call the constructor of the Vector base class
        super(Position, self).__init__(x, y)

    # -----------------------------------------------------------------

    def __sub__(self, position):

        """
        This function subtracts two positions to obtain an Extent
        """

        return Extent(self.x - position.x, self.y - position.y)

# -----------------------------------------------------------------

class Position3D(Vector3D):

    """
    This classs ...
    """

    def __init__(self, x, y, z):

        """
        This function ...
        :param self:
        :param x:
        :param y:
        :param z:
        :return:
        """

        # Call the constructor of the Vector3D base class
        super(Position3D, self).__init__(x, y, z)

# -----------------------------------------------------------------

class Extent(Vector):

    """
    This class ...
    """

    def __init__(self, x, y):

        """
        The constructor ...
        :param x:
        :param y:
        :return:
        """

        # Call the constructor of the Vector base class
        super(Extent, self).__init__(x, y)

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :param y:
        :return:
        """

        return self.__class__(self.x + extent.x, self.y + extent.y)

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return self.__class__(self.x - extent.x, self.y - extent.y)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        #print(self.x, type(self.x))
        #print(self.y, type(self.y))
        #print(value, type(value))
        return self.__class__(self.x * value, self.y * value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__class__(self.x / value, self.y / value)

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__class__(self.x / value, self.y / value)

# -----------------------------------------------------------------

class IntegerExtent(Extent):

    """
    This function ...
    """

    def __init__(self, x, y):

        """
        The constructor ..
        :param x:
        :param y:
        """

        # Call the constructor of the base class
        super(IntegerExtent, self).__init__(x, y)

        # Check
        if not types.is_integer_type(self.x): raise ValueError("x and y must be integers")
        if not types.is_integer_type(self.y): raise ValueError("x and y must be integers")

# -----------------------------------------------------------------

class RealExtent(Extent):

    """
    This function ...
    """

    def __init__(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        """

        # Call the constructor of the base class
        super(RealExtent, self).__init__(x, y)

        # Check
        if not types.is_real_type(self.x): raise ValueError("x and y must be real values")
        if not types.is_real_type(self.y): raise ValueError("x and y must be real values")

# -----------------------------------------------------------------

class AngleExtent(Extent):

    """
    This class ...
    """

    def __init__(self, x, y):

        """
        The constructor ...
        :param x:
        :param y:
        """

        # Call the constructor of the base calss
        super(AngleExtent, self).__init__(x, y)

        # Convert quantities to angles
        from ...core.units.quantity import quantity_to_angle
        if types.is_quantity(self.x): self.x = quantity_to_angle(self.x)
        if types.is_quantity(self.y): self.y = quantity_to_angle(self.y)

        # Check
        if not types.is_angle(self.x): raise ValueError("Arguments must be angles")
        if not types.is_angle(self.y): raise ValueError("Arguments must be angles")

# -----------------------------------------------------------------

class QuantityExtent(Extent):

    """
    This function ...
    """

    def __init__(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        """

        # Call the constructor of the base class
        super(QuantityExtent, self).__init__(x, y)

        # Check
        if not types.is_quantity(self.x): raise ValueError("x and y must be quantities")
        if not types.is_quantity(self.y): raise ValueError("x and y must be quantities")

        # Check whether same physical type
        from ...core.units.quantity import same_physical_type
        if not same_physical_type(self.x, self.y): raise ValueError("x and y must represent the same physical type")

# -----------------------------------------------------------------

class Pixel(Vector):

    """
    This class ...
    """

    def __init__(self, x, y):

        """
        The constructor ...
        """

        # Check the arguments
        if not types.is_integer_type(x): raise ValueError("Arguments must be integer numbers")
        if not types.is_integer_type(y): raise ValueError("Arguments must be integer numbers")

        # Call the constructor of the base class
        super(Pixel, self).__init__(x, y)

    # -----------------------------------------------------------------

    def exists_in(self, frame):

        """
        This function ...
        :param frame: 
        :return: 
        """

        # Check whether the position can be within the mask considering its dimensions
        if self.x < 0 or self.y < 0 or self.x >= frame.xsize or self.y >= frame.ysize: return False
        else: return True

    # -----------------------------------------------------------------

    @classmethod
    def for_coordinate(cls, coordinate, round_first=False):

        """
        This function ...
        :param coordinate:
        :param round_first:
        :return:
        """

        if round_first: return cls(int(round(coordinate.x)), int(round(coordinate.y)))
        else: return cls(int(coordinate.x), int(coordinate.y))

    # -----------------------------------------------------------------

    @classmethod
    def center_of(cls, data):

        """
        This function ...
        :param data:
        :return:
        """

        # Get the shape
        ny, nx = data.shape

        # Get center
        x = int(math.round((nx - 1.) / 2.))
        y = int(math.round((ny - 1.) / 2.))

        # Create
        return cls(x, y)

# -----------------------------------------------------------------

class PixelShape(tuple):

    """
    This function ...
    """

    def __new__(cls, y, x):

        """
        This function ...
        :param a:
        :param b:
        :return:
        """

        # Call the constructor of the base class
        return super(PixelShape, cls).__new__(cls, tuple([y, x]))

    # -----------------------------------------------------------------

    @classmethod
    def square(cls, npixels):

        """
        This function ...
        :param npixels:
        :return:
        """

        return cls(npixels, npixels)

    # -----------------------------------------------------------------

    @classmethod
    def from_tuple(cls, shape):

        """
        This function ...
        :param shape:
        :return:
        """

        return cls(x=shape[1], y=shape[0])

    # -----------------------------------------------------------------

    @classmethod
    def from_xy(cls, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        return cls(x=x, y=y)

    # -----------------------------------------------------------------

    @classmethod
    def from_xy_tuple(cls, xy):

        """
        This function ...
        :param xy:
        :return:
        """

        return cls(x=xy[0], y=xy[1])

    # -----------------------------------------------------------------

    @classmethod
    def from_yx(cls, y, x):

        """
        This function ...
        :param y:
        :param x:
        :return:
        """

        return cls(x=x, y=y)

    # -----------------------------------------------------------------

    @classmethod
    def from_yx_tuple(cls, yx):

        """
        This functino ...
        :param yx:
        :return:
        """

        return cls(x=yx[1], y=xy[0])

    # -----------------------------------------------------------------

    @property
    def x(self):

        """
        This function ...
        :return:
        """

        return self[1]

    # -----------------------------------------------------------------

    @property
    def y(self):

        """
        This function ...
        :return:
        """

        return self[0]

    # -----------------------------------------------------------------

    @property
    def nx(self):

        """
        This function ...
        :return:
        """

        return self.x

    # -----------------------------------------------------------------

    @property
    def ny(self):

        """
        This function ...
        :return:
        """

        return self.y

    # -----------------------------------------------------------------

    @property
    def nxpixels(self):

        """
        This function ...
        :return:
        """

        return self.nx

    # -----------------------------------------------------------------

    @property
    def nypixels(self):

        """
        This function ...
        :return:
        """

        return self.ny

    # -----------------------------------------------------------------

    @property
    def xpixels(self):

        """
        This function ...
        :return:
        """

        return self.nx

    # -----------------------------------------------------------------

    @property
    def ypixels(self):

        """
        This function ...
        :return:
        """

        return self.ny

    # -----------------------------------------------------------------

    @property
    def xy(self):

        """
        This function ...
        :return:
        """

        return self.nx * self.ny

    # -----------------------------------------------------------------

    @property
    def nxy(self):

        """
        This function ...
        :return:
        """

        return self.xy

    # -----------------------------------------------------------------

    @property
    def ntotal(self):

        """
        This function ...
        :return:
        """

        return self.nxy

    # -----------------------------------------------------------------

    @property
    def ntotalpixels(self):

        """
        This function ...
        :return:
        """

        return self.nxy

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return "(nx=" + str(self.nx) + ", ny=" + str(self.ny) + ")"

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return "PixelShape(x=" + str(self.nx) + ", y=" + str(self.ny) + ")"

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :return:
        """

        return self.x == other.x and self.y == other.y

    # -----------------------------------------------------------------

    def __lt__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self.ntotal < other.ntotal

    # -----------------------------------------------------------------

    def __gt__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self.ntotal > other.ntotal

    # -----------------------------------------------------------------

    def __le__(self, other):

        """
        This function ...
        :return:
        """

        return self.ntotal <= other.ntotal

    # -----------------------------------------------------------------

    def __ge__(self, other):

        """
        This function ...
        :return:
        """

        return self.ntotal >= other.ntotal

    # -----------------------------------------------------------------

    def as_tuple(self):

        """
        This function returns the shape as a regular tuple
        :return:
        """

        return (self.y, self.x,)

# -----------------------------------------------------------------
