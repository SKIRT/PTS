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
from abc import ABCMeta, abstractproperty

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite
from .projection import GalaxyProjection
from ...core.tools.utils import abstractclassmethod

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

    @abstractclassmethod
    def from_projection(self, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @classmethod
    def from_deprojection(cls, deprojection):

        """
        This function ...
        :param deprojection:
        :return:
        """

        return cls.from_projection(deprojection.projection)

    # -----------------------------------------------------------------

    @classmethod
    def from_deprojection_faceon(cls, deprojection):

        """
        This function ...
        :param deprojection:
        :return:
        """

        return cls.from_projection(deprojection.faceon_projection)

    # -----------------------------------------------------------------

    @classmethod
    def from_deprojection_edgeon(cls, deprojection):

        """
        This function ...
        :param deprojection:
        :return:
        """

        return cls.from_projection(deprojection.edgeon_projection)

    # -----------------------------------------------------------------

    @classmethod
    def from_wcs(cls, wcs, center, distance, inclination, azimuth, position_angle):

        """
        This function ...
        :param wcs:
        :param center:
        :param distance:
        :param inclination:
        :param azimuth:
        :param position_angle:
        :return:
        """

        # Create projection
        # wcs, center, distance, inclination, azimuth, position_angle
        projection = GalaxyProjection.from_wcs(wcs, center, distance, inclination, azimuth, position_angle)

        # Create and return
        return cls.from_projection(projection)

    # -----------------------------------------------------------------

    @abstractproperty
    def npixels(self):

        """
        This function ...
        :return:
        """

        pass

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
    def from_properties(cls, distance, inclination, position_angle, azimuth=None):

        """
        This function ...
        :param distance:
        :param inclination:
        :param position_angle:
        :param azimuth:
        :return:
        """

        # Set azimuth
        if azimuth is None: azimuth = Angle(0.0, "deg")

        # Create and return the instrument
        return cls(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle)

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

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return 1

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

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return self.pixels_x * self.pixels_y

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

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return self.pixels_x * self.pixels_y

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
        self.add_property("counts", "boolean", "write photon counts", False)

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection, scattering_levels=0, counts=False):

        """
        This function ...
        :param projection:
        :param scattering_levels:
        :param counts:
        :return:
        """

        return cls(distance=projection.distance, inclination=projection.inclination, azimuth=projection.azimuth,
                   position_angle=projection.position_angle, field_x=projection.field_x_physical,
                   field_y=projection.field_y_physical, pixels_x=projection.pixels_x, pixels_y=projection.pixels_y,
                   center_x=projection.center_x, center_y=projection.center_y, scattering_levels=scattering_levels, counts=counts)

    # -----------------------------------------------------------------

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return self.pixels_x * self.pixels_y

# -----------------------------------------------------------------

class FullSEDInstrument(Instrument):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FullSEDInstrument, self).__init__()

        # The number of scattering levels
        self.add_property("scattering_levels", "integer", "number of scattering levels", 0)

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_properties(cls, distance, inclination, position_angle, azimuth=None, scattering_levels=0):

        """
        This function ...
        :param distance:
        :param inclination:
        :param position_angle:
        :param azimuth:
        :param scattering_levels:
        :return:
        """

        # Set azimuth
        if azimuth is None: azimuth = Angle(0.0, "deg")

        # Create and return the instrument
        return cls(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                   scattering_levels=scattering_levels)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection, scattering_levels=0):

        """
        This function ...
        :param projection:
        :param scattering_levels:
        :return:
        """

        return cls(distance=projection.distance, inclination=projection.inclination, azimuth=projection.azimuth,
                   position_angle=projection.position_angle, scattering_levels=scattering_levels)

    # -----------------------------------------------------------------

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return 1

# -----------------------------------------------------------------

class MultiFrameInstrument(Instrument):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MultiFrameInstrument, self).__init__()

        # Add basic properties
        self.add_property("write_total", "boolean", "write total", False)
        self.add_property("write_stellar_components", "boolean", "write stellar components", False)

        # Add the frames property
        self.add_property("frames", "instrument_frame_list", "list of instrument frames", [])

        # Add the frames
        for frame in args: self.add_frame(frame)

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented for multi frame instrument")

    # -----------------------------------------------------------------

    def add_frame(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Add the frame
        self.frames.append(frame)

    # -----------------------------------------------------------------

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        raise ValueError("The number of pixels cannot be determined")

# -----------------------------------------------------------------

class InstrumentFrame(SimplePropertyComposite):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(InstrumentFrame, self).__init__()

        # Define properties
        self.add_property("pixels_x", "positive_integer", "x pixels", None)
        self.add_property("pixels_y", "positive_integer", "y pixels", None)
        self.add_property("field_x", "quantity", "x field", None)
        self.add_property("field_y", "quantity", "y field", None)

        # Set properties
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return self.pixels_x * self.pixels_y

# -----------------------------------------------------------------
