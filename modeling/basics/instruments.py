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

    def __init__(self):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(Instrument, self).__init__()

        # Define properties
        self.add_property("distance", "quantity", "distance", None)
        self.add_property("inclination", "angle", "inclination", None)
        self.add_property("azimuth", "angle", "azimuth", None)
        self.add_property("position_angle", "angle", "position angle", None)

# -----------------------------------------------------------------

class SEDInstrument(Instrument):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param **kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SEDInstrument, self).__init__()

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        return cls(distance=projection.distance, inclination=projection.inclination, azimuth=projection.azimuth, position_angle=projection.position_angle)

# -----------------------------------------------------------------

class FrameInstrument(Instrument):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FrameInstrument, self).__init__()

        # Define properties
        self.add_property("field_x", "quantity", "x field")
        self.add_property("field_y", "quantity", "y field")
        self.add_property("pixels_x", "positive_integer", "x pixels")
        self.add_property("pixels_y", "positive_integer", "y pixels")
        self.add_property("center_x", "quantity", "x center")
        self.add_property("center_y", "quantity", "y center")

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        return cls(distance=projection.distance, inclination=projection.inclination, azimuth=projection.azimuth,
                   position_angle=projection.position_angle, field_x=projection.field_x_physical,
                   field_y=projection.field_y_physical, pixels_x=projection.pixels_x, pixels_y=projection.pixels_y,
                   center_x=projection.center_x, center_y=projection.center_y)

# -----------------------------------------------------------------

class SimpleInstrument(Instrument):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SimpleInstrument, self).__init__()

        # Define properties
        self.add_property("field_x", "quantity", "x field")
        self.add_property("field_y", "quantity", "y field")
        self.add_property("pixels_x", "positive_integer", "x pixels")
        self.add_property("pixels_y", "positive_integer", "y pixels")
        self.add_property("center_x", "quantity", "x center")
        self.add_property("center_y", "quantity", "y center")

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        return cls(distance=projection.distance, inclination=projection.inclination, azimuth=projection.azimuth,
                   position_angle=projection.position_angle, field_x=projection.field_x_physical,
                   field_y=projection.field_y_physical, pixels_x=projection.pixels_x, pixels_y=projection.pixels_y,
                   center_x=projection.center_x, center_y=projection.center_y)

# -----------------------------------------------------------------

class FullInstrument(Instrument):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FullInstrument, self).__init__()

        # Define properties
        self.add_property("field_x", "quantity", "x field")
        self.add_property("field_y", "quantity", "y field")
        self.add_property("pixels_x", "positive_integer", "x pixels")
        self.add_property("pixels_y", "positive_integer", "y pixels")
        self.add_property("center_x", "quantity", "x center")
        self.add_property("center_y", "quantity", "y center")
        self.add_property("scattering_levels", "integer", "number of scattering levels", 0) # default=0

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        return cls(distance=projection.distance, inclination=projection.inclination, azimuth=projection.azimuth,
                   position_angle=projection.position_angle, field_x=projection.field_x_physical,
                   field_y=projection.field_y_physical, pixels_x=projection.pixels_x, pixels_y=projection.pixels_y,
                   center_x=projection.center_x, center_y=projection.center_y)

# -----------------------------------------------------------------
