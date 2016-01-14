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

# Import the relevant PTS classes and modules
from .vector import Position, Extent

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
