#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.projection Contains the GalaxyProjection, FaceOnProjection and EdgeOnProjection classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Unit, dimensionless_angles

# Import the relevant PTS classes and modules
from ...magic.basics.pixelscale import Pixelscale
from ...magic.basics.vector import Position

# -----------------------------------------------------------------

# SKIRT:  incl.  azimuth PA
# XY-plane	0	 0	    90
# XZ-plane	90	 -90	0
# YZ-plane	90	 0	    0

# -----------------------------------------------------------------

class GalaxyProjection(object):

    """
    This function ...
    """

    def __init__(self, distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y):

        """
        This function ...
        :param distance:
        :param inclination:
        :param azimuth:
        :param position_angle:
        :param pixels_x:
        :param pixels_y:
        :param center_x:
        :param center_y:
        :param field_x:
        :param field_y:
        """

        self.distance = distance
        self.inclination = Angle(inclination, "deg")
        self.azimuth = Angle(azimuth, "deg")
        self.position_angle = Angle(position_angle, "deg")
        self.pixels_x = pixels_x
        self.pixels_y = pixels_y
        self.center_x = center_x
        self.center_y = center_y
        self.field_x_physical = field_x
        self.field_y_physical = field_y

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        # From field of view to pixel scale
        pixelscale_x  = self.field_x_physical / self.pixels_x
        pixelscale_y = self.field_y_physical / self.pixels_y

        # Create and return the pixelscale
        return Pixelscale(pixelscale_x, pixelscale_y)

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

        # Get derived properties
        pixels_x, pixels_y, center_x, center_y, field_x, field_y = get_relevant_wcs_properties(wcs, center, distance)

        # Create and return a new class instance
        return cls(distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_deprojection(cls, deprojection, distance, azimuth):

        """
        This function ...
        :param deprojection:
        :param distance:
        :param azimuth:
        :return:
        """

        pixels_x = deprojection.x_size
        pixels_y = deprojection.y_size

        pixelscale = deprojection.pixelscale

        # In pixel
        x_center = deprojection.x_center
        y_center = deprojection.y_center

        # To physical
        center_x = x_center * pixelscale
        center_y = y_center * pixelscale

        # Detemrine field of view
        field_x = pixelscale * pixels_x
        field_y = pixelscale * pixels_y

        inclination = deprojection.inclination
        position_angle = deprojection.position_angle

        # Create and return a new class instance
        return cls(distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_instrument(cls, instrument, default_pixels_x=None, default_pixels_y=None, default_field_x=None,
                        default_field_y=None):

        """
        This function ...
        :param instrument:
        :param default_pixels_x:
        :param default_pixels_y:
        :param default_field_x:
        :param default_field_y:
        :return:
        """

        # distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y

        # Each instrument
        distance = instrument.distance
        inclination = instrument.inclination
        azimuth = instrument.azimuth
        position_angle = instrument.position_angle

        from .instruments import SEDInstrument, FrameInstrument, SimpleInstrument, FullInstrument

        # SED instrument
        if isinstance(instrument, SEDInstrument):

            pixels_x = default_pixels_x
            pixels_y = default_pixels_y

            # Put the galaxy center at the center of the instrument
            center_x = 0.5 * (pixels_x + 1)
            center_y = 0.5 * (pixels_y + 1)

            # Determine the pixelscale
            pixelscale_physical_x = default_field_x / pixels_x
            pixelscale_physical_y = default_field_y / pixels_y

            # Calculate center in physical units
            center_x = center_x * pixelscale_physical_x
            center_y = center_y * pixelscale_physical_y

            # Set field
            field_x = default_field_x
            field_y = default_field_y

        # Instrument with pixels
        elif isinstance(instrument, FrameInstrument) or isinstance(instrument, SimpleInstrument) or isinstance(instrument, FullInstrument):

            pixels_x = instrument.pixels_x
            pixels_y = instrument.pixels_y
            center_x = instrument.center_x
            center_y = instrument.center_y
            field_x = instrument.field_x
            field_y = instrument.field_y

        # Not recognized
        else: raise ValueError("Unrecognized instrument object")

        # Create and return a new class instance
        return cls(distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        distance = None
        inclination = None
        azimuth = None
        position_angle = None
        pixels_x = None
        pixels_y = None
        center_x = None
        center_y = None
        field_x = None
        field_y = None

        # Read the projection file
        with open(path, 'r') as projection_file:

            # Loop over all lines in the file
            for line in projection_file:

                splitted = line.split(": ")

                if splitted[0] == "Distance": distance = get_quantity(splitted[1])
                elif splitted[0] == "Inclination": inclination = get_angle(splitted[1])
                elif splitted[0] == "Azimuth": azimuth = get_angle(splitted[1])
                elif splitted[0] == "Position angle": position_angle = get_angle(splitted[1])
                elif splitted[0] == "Pixels x": pixels_x = int(splitted[1])
                elif splitted[0] == "Pixels y": pixels_y = int(splitted[1])
                elif splitted[0] == "Center x": center_x = get_quantity(splitted[1])
                elif splitted[0] == "Center y": center_y = get_quantity(splitted[1])
                elif splitted[0] == "Field x": field_x = get_quantity(splitted[1])
                elif splitted[0] == "Field y": field_y = get_quantity(splitted[1])

        # Create and return a new class instance
        projection = cls(distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

        # Set the path
        projection.path = path

        # Return the projection
        return projection

    # -----------------------------------------------------------------

    def to_wcs(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether the path is valid
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create the projection file
        with open(path, 'w') as projection_file:

            print("Distance:", str(self.distance), file=projection_file)
            print("Inclination:", str(self.inclination.to("deg").value) + " deg", file=projection_file)
            print("Azimuth:", str(self.azimuth.to("deg").value) + " deg", file=projection_file)
            print("Position angle:", str(self.position_angle.to("deg").value) + " deg", file=projection_file)
            print("Pixels x:", str(self.pixels_x), file=projection_file)
            print("Pixels y:", str(self.pixels_y), file=projection_file)
            print("Center x:", str(self.center_x), file=projection_file)
            print("Center y:", str(self.center_y), file=projection_file)
            print("Field x:", str(self.field_x_physical), file=projection_file)
            print("Field y:", str(self.field_y_physical), file=projection_file)

        # Update the path
        self.path = path

# -----------------------------------------------------------------

class FaceOnProjection(GalaxyProjection):

    """
    This class ...
    """

    def __init__(self, distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(FaceOnProjection, self).__init__(distance, 0.0, 0.0, 90., pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_wcs(cls, wcs, center, distance):

        """
        This function ...
        :param wcs:
        :param center:
        :param distance:
        :return:
        """

        # Get derived properties
        pixels_x, pixels_y, center_x, center_y, field_x, field_y = get_relevant_wcs_properties(wcs, center, distance)

        # Call the constructor
        return cls(distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_deprojection(cls, deprojection, distance):

        """
        This function ...
        :param deprojection:
        :param distance:
        :return:
        """

        pixels_x = deprojection.x_size
        pixels_y = deprojection.y_size

        pixelscale = deprojection.pixelscale

        # In pixel
        x_center = deprojection.x_center
        y_center = deprojection.y_center

        # To physical
        center_x = x_center * pixelscale
        center_y = y_center * pixelscale

        # Detemrine field of view
        field_x = pixelscale * pixels_x
        field_y = pixelscale * pixels_y

        # Create and return
        return cls(distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        # Create and return
        return cls(projection.distance, projection.pixels_x, projection.pixels_y, projection.center_x,
                   projection.center_y, projection.field_x_physical, projection.field_y_physical)

# -----------------------------------------------------------------

class EdgeOnProjection(GalaxyProjection):

    """
    This class ...
    """

    def __init__(self, distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y):

        """
        The constructor ...
        :param distance:
        :param pixels_x:
        :param pixels_y:
        :param center_x:
        :param center_y:
        :param field_x:
        :param field_y:
        """

        # Call the constructor of the base class
        super(EdgeOnProjection, self).__init__(distance, 90., 0., 0., pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_wcs(cls, wcs, center, distance):

        """
        This function ...
        :param wcs:
        :param center:
        :param distance:
        :return:
        """

        # Get derived properties
        pixels_x, pixels_y, center_x, center_y, field_x, field_y = get_relevant_wcs_properties(wcs, center, distance)

        # Call the constructor
        return cls(distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_deprojection(cls, deprojection, distance):

        """
        This function ...
        :param deprojection:
        :param distance:
        :return:
        """

        pixels_x = deprojection.x_size
        pixels_y = deprojection.y_size

        pixelscale = deprojection.pixelscale

        # In pixel
        x_center = deprojection.x_center
        y_center = deprojection.y_center

        # To physical
        center_x = x_center * pixelscale
        center_y = y_center * pixelscale

        # Detemrine field of view
        field_x = pixelscale * pixels_x
        field_y = pixelscale * pixels_y

        # Call the constructor
        return cls(distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        # Create and return
        return cls(projection.distance, projection.pixels_x, projection.pixels_y, projection.center_x,
                   projection.center_y, projection.field_x_physical, projection.field_y_physical)

# -----------------------------------------------------------------

def get_quantity(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create a quantity object and return it
    if unit is not None: value = value * Unit(unit)
    return value

# -----------------------------------------------------------------

def get_angle(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    if "d" in entry and "m" in entry and "s" in entry: # dms format

        value = float(entry.split("d")[0]) + float(entry.split("d")[1].split("m")[0])/60. + float(entry.split("m")[1].split("s")[0])/3600.
        unit = "deg"

    else:

        splitted = entry.split()
        value = float(splitted[0])
        try: unit = splitted[1]
        except IndexError: unit = default_unit

    # Create an Angle object and return it
    if unit is not None: value = Angle(value, unit)
    return value

# -----------------------------------------------------------------

def get_relevant_wcs_properties(wcs, center, distance):

    """
    This function ...
    :param wcs:
    :param center:
    :param distance:
    :return:
    """

    # PIXEL SIZE
    pixels_x = wcs.xsize
    pixels_y = wcs.ysize

    # CENTER PIXEL
    pixel_center = center.to_pixel(wcs)
    center = Position(0.5*pixels_x - pixel_center.x - 0.5, 0.5*pixels_y - pixel_center.y - 0.5)
    center_x = center.x
    center_y = center.y
    center_x = (center_x * wcs.pixelscale.x.to("deg") * distance).to("pc", equivalencies=dimensionless_angles())
    center_y = (center_y * wcs.pixelscale.y.to("deg") * distance).to("pc", equivalencies=dimensionless_angles())

    # FIELD OF VIEW
    field_x_angular = wcs.pixelscale.x.to("deg") * pixels_x
    field_y_angular = wcs.pixelscale.y.to("deg") * pixels_y
    field_x_physical = (field_x_angular * distance).to("pc", equivalencies=dimensionless_angles())
    field_y_physical = (field_y_angular * distance).to("pc", equivalencies=dimensionless_angles())

    return pixels_x, pixels_y, center_x, center_y, field_x_physical, field_y_physical

# -----------------------------------------------------------------
