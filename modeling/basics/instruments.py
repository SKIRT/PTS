#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.instruments Contains the SEDInstrument, FrameInstrument SimpleInstrument and FullInstruemnt classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ...core.tools import parsing
from ...core.basics.configuration import stringify_not_list
from ...core.tools.logging import log

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

class Instrument(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading instrument from " + path + " ...")

        properties = dict()

        with open(path, 'r') as instrumentfile:

            for line in instrumentfile:

                if "Type:" in line: continue

                name, rest = line.split(": ")
                value, dtype = rest.split("[")
                dtype = dtype.split("]")[0]

                # Set the property value
                properties[name] = getattr(parsing, dtype)(value)

        # Create the class instance
        return cls(**properties)

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the instrument to " + path + " ...")

        # Write the properties
        with open(path, 'w') as instrumentfile:

            # Print the type
            print("Type:", self.__class__.__name__, file=instrumentfile)

            # Loop over the variables
            for name in vars(self):

                dtype, value = stringify_not_list(getattr(self, name))
                print(name + ":", value + " [" + dtype + "]", file=instrumentfile)

    # -----------------------------------------------------------------

    def copy(self):
        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

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
