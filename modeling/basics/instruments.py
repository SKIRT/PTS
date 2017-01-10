#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.instruments Contains the SEDInstrument, FrameInstrument, SimpleInstrument and
#  FullInstrument classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

def load_instrument(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Get the first line of the file
    with open(path, 'r') as f: first_line = f.readline()

    # Create the appropriate model
    if "SEDInstrument" in first_line: return SEDInstrument.from_file(path)
    elif "FrameInstrument" in first_line: return FrameInstrument.from_file(path)
    elif "SimpleInstrument" in first_line: return SimpleInstrument.from_file(path)
    elif "FullInstrument" in first_line: return FullInstrument.from_file(path)
    else: raise ValueError("Unrecognized instrument file")

# -----------------------------------------------------------------

class Instrument(SimplePropertyComposite):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

# -----------------------------------------------------------------

class SEDInstrument(Instrument):

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

        # Call the constructor of the base class
        super(SEDInstrument, self).__init__()

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

        return cls(projection.distance, projection.inclination, projection.azimuth, projection.position_angle)

# -----------------------------------------------------------------

class FrameInstrument(Instrument):

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

        # Call the constructor of the base class
        super(FrameInstrument, self).__init__()

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

        return cls(projection.distance, projection.inclination, projection.azimuth, projection.position_angle,
                   projection.field_x_physical, projection.field_y_physical, projection.pixels_x, projection.pixels_y,
                   projection.center_x, projection.center_y)

# -----------------------------------------------------------------

class SimpleInstrument(Instrument):

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

        # Call the constructor of the base class
        super(SimpleInstrument, self).__init__()

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

        return cls(projection.distance, projection.inclination, projection.azimuth, projection.position_angle, projection.field_x_physical,
                   projection.field_y_physical, projection.pixels_x, projection.pixels_y, projection.center_x, projection.center_y)

# -----------------------------------------------------------------

class FullInstrument(Instrument):

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

        # Call the constructor of the base class
        super(FullInstrument, self).__init__()

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

        return cls(projection.distance, projection.inclination, projection.azimuth, projection.position_angle,
                   projection.field_x_physical, projection.field_y_physical, projection.pixels_x, projection.pixels_y,
                   projection.center_x, projection.center_y)

# -----------------------------------------------------------------
