#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import math
import numpy as np

# Import astronomical units
from astropy import units as u

# *****************************************************************

class Vector(object):

    """
    This class ...
    """

    # Define the symbol for displaying as string
    symbol = "V"

    # *****************************************************************

    def __init__(self, x, y):

        """
        The constructor ...
        """

        # Set the attributes
        self.x = x
        self.y = y

    # *****************************************************************

    @property
    def norm(self):

        """
        This function ...
        """

        return math.sqrt(self.x**2 + self.y**2)

    # *****************************************************************

    @property
    def angle(self):

        """
        This function ...
        """

        return math.atan2(self.x, self.y)

    # *****************************************************************

    def polar(self):

        """
        This function returns the polar coordinates (radius, angle) of this vector
        """

        return self.norm, self.angle

    # *****************************************************************

    def __str__(self):

        """
        This function ...
        """

        return self.symbol + '(x={0}, y={1})'.format(self.x, self.y)

     # *****************************************************************

    def __repr__(self):

        """
        This function
        """

        return '<' + self.__class__.__name__ + ' x={0}, y={1}>'.format(self.x, self.y)

# *****************************************************************

class Position(Vector):

    """
    This class ...
    """

    # Define the symbol for displaying as string
    symbol = "P"

    # *****************************************************************

    def __init__(self, x, y):

        """
        The constructor ...
        """

        # Call the constructor of the Vector base class
        super(Position, self).__init__(x, y)

    # *****************************************************************

    def __sub__(self, position):

        """
        This function subtracts two positions to obtain an Extent
        """

        return Extent(self.x - position.x, self.y - position.y)

# *****************************************************************

class Extent(Vector):

    """
    This class ...
    """

    # Define the symbol for displaying as string
    symbol = "E"

    # *****************************************************************

    def __init__(self, x, y):

        """
        The constructor ...
        :param x:
        :param y:
        :return:
        """

        # Call the constructor of the Vector base class
        super(Extent, self).__init__(x, y)

    # *****************************************************************

    def __add__(self, extent):

        """
        This function ...
        :param y:
        :return:
        """

        return Extent(self.x + extent.x, self.y + extent.y)

    # *****************************************************************

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return Extent(self.x - extent.x, self.y - extent.y)

    # *****************************************************************

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Extent(self.x * value, self.y * value)

    # *****************************************************************

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Extent(self.x / value, self.y / value)

# *****************************************************************
