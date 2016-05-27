#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.instruments Contains the SimpleInstrument class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class SEDInstrument(object):

    """
    This class ...
    """

    def __init__(self, distance, inclination, azimuth, position_angle):

        """
        This function ...
        :param distance:
        :param inclination:
        :param azimuth:
        :param position_angle:
        :return:
        """

        self.distance = distance
        self.inclination = inclination
        self.azimuth = azimuth
        self.position_angle = position_angle

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        return cls(projection.distance, projection.azimuth, projection.position_angle)

# -----------------------------------------------------------------

class FrameInstrument(object):

    """
    This class ...
    """

    def __init__(self, distance, inclination, azimuth, position_angle, field_x, field_y, pixels_x, pixels_y, center_x, center_y):

        """
        This function ...
        :param distance:
        :param inclination:
        :param azimuth:
        :param position_angle:
        :param field_x:
        :param field_y:
        :param pixels_x:
        :param pixels_y:
        :param center_x:
        :param center_y:
        """

        self.distance = distance
        self.inclination = inclination
        self.azimuth = azimuth
        self.position_angle = position_angle
        self.field_x = field_x
        self.field_y = field_y
        self.pixels_x = pixels_x
        self.pixels_y = pixels_y
        self.center_x = center_x
        self.center_y = center_y

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):
        """
        This function ...
        :param projection:
        :return:
        """

        return cls(projection.distance, projection.azimuth, projection.position_angle, projection.field_x_physical,
                   projection.field_y_physical, projection.pixels_x, projection.pixels_y, projection.center_x,
                   projection.center_y)

# -----------------------------------------------------------------

class SimpleInstrument(object):

    """
    This class ...
    """

    def __init__(self, distance, inclination, azimuth, position_angle, field_x, field_y, pixels_x, pixels_y, center_x, center_y):

        """
        This function ...
        :param distance:
        :param inclination:
        :param azimuth:
        :param position_angle:
        :param field_x:
        :param field_y:
        :param pixels_x:
        :param pixels_y:
        :param center_x:
        :param center_y:
        :return:
        """

        self.distance = distance
        self.inclination = inclination
        self.azimuth = azimuth
        self.position_angle = position_angle
        self.field_x = field_x
        self.field_y = field_y
        self.pixels_x = pixels_x
        self.pixels_y = pixels_y
        self.center_x = center_x
        self.center_y = center_y

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        return cls(projection.distance, projection.azimuth, projection.position_angle, projection.field_x_physical,
                   projection.field_y_physical, projection.pixels_x, projection.pixels_y, projection.center_x, projection.center_y)

# -----------------------------------------------------------------

class FullInstrument(object):

    """
    This class ...
    """

    def __init__(self, distance, inclination, azimuth, position_angle, field_x, field_y, pixels_x, pixels_y, center_x, center_y):

        """
        This function ...
        :param distance:
        :param inclination:
        :param azimuth:
        :param position_angle:
        :param field_x:
        :param field_y:
        :param pixels_x:
        :param pixels_y:
        :param center_x:
        :param center_y:
        :return:
        """

        self.distance = distance
        self.inclination = inclination
        self.azimuth = azimuth
        self.position_angle = position_angle
        self.field_x = field_x
        self.field_y = field_y
        self.pixels_x = pixels_x
        self.pixels_y = pixels_y
        self.center_x = center_x
        self.center_y = center_y

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        return cls(projection.distance, projection.azimuth, projection.position_angle, projection.field_x_physical,
                   projection.field_y_physical, projection.pixels_x, projection.pixels_y, projection.center_x,
                   projection.center_y)

# -----------------------------------------------------------------
