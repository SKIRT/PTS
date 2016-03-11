#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.vector Contains the Vector class and the subclasses Position and Extent.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

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

        """
        This function ...
        """

        return math.sqrt(self.x**2 + self.y**2)

    # -----------------------------------------------------------------

    @property
    def angle(self):

        """
        This function ...
        """

        return math.atan2(self.x, self.y)

    # -----------------------------------------------------------------

    def polar(self):

        """
        This function returns the polar coordinates (radius, angle) of this vector
        """

        return self.norm, self.angle

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

        return Extent(self.x + extent.x, self.y + extent.y)

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Extent(self.x - extent.x, self.y - extent.y)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Extent(self.x * value, self.y * value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Extent(self.x / value, self.y / value)

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Extent(self.x / value, self.y / value)

# -----------------------------------------------------------------