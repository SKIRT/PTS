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

# *****************************************************************

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
