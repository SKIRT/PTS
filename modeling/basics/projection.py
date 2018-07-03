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
from ...magic.basics.pixelscale import Pixelscale, PhysicalPixelscale
from ...magic.basics.vector import Position
from ...magic.basics.vector import Extent
from ...magic.basics.stretch import PhysicalExtent
from ...magic.basics.coordinate import PixelCoordinate, PhysicalCoordinate
from ...core.tools import types
from ...core.basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

# SKIRT:  incl.  azimuth PA
# XY-plane	0	 0	    90
# XZ-plane	90	 -90	0
# YZ-plane	90	 0	    0

# -----------------------------------------------------------------

def load_projection(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Get the first line of the file
    with open(path, 'r') as f: first_line = f.readline()

    # Create the appropriate model
    if "GalaxyProjection" in first_line: return GalaxyProjection.from_file(path)
    elif "FaceOnProjection" in first_line: return FaceOnProjection.from_file(path)
    elif "EdgeOnProjection" in first_line: return EdgeOnProjection.from_file(path)
    else: return GalaxyProjection.from_old_file(path)
    #raise ValueError("Unrecognized instrument file")

# -----------------------------------------------------------------

class GalaxyProjection(SimplePropertyComposite):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GalaxyProjection, self).__init__()

        # Define properties
        self.add_property("distance", "quantity", "distance", None)
        self.add_property("inclination", "angle", "inclination", None)
        self.add_property("azimuth", "angle", "azimuth", None)
        self.add_property("position_angle", "angle", "position angle", None)
        self.add_property("pixels_x", "positive_integer", "number of pixels along the horizontal axis", None)
        self.add_property("pixels_y", "positive_integer", "number of pixels along the vertical axis", None)
        self.add_property("center_x", "length_quantity", "x position of the galaxy center in the frame", None)
        self.add_property("center_y", "length_quantity", "y position of the galaxy center in the frame", None)
        self.add_property("field_x", "length_quantity", "x field of view", None)
        self.add_property("field_y", "length_quantity", "y field of view", None)

        # Set values
        self.set_properties(kwargs)

    # -----------------------------------------------------------------

    # compatibility
    @classmethod
    def from_file(cls, path, remote=None):

        """
        Thisf unction ...
        :param path:
        :param remote:
        :return:
        """

        try: return super(GalaxyProjection, cls).from_file(path, remote=remote)
        except ValueError: return cls.from_old_file(path)

    # -----------------------------------------------------------------

    @classmethod
    def from_old_file(cls, path):

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

        # Check if all defined
        assert distance is not None
        assert inclination is not None
        assert azimuth is not None
        assert position_angle is not None
        assert pixels_x is not None
        assert pixels_y is not None
        assert center_x is not None
        assert center_y is not None
        assert field_x is not None
        assert field_y is not None

        # Create and return
        projection = cls(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                   pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x, center_y=center_y, field_x=field_x, field_y=field_y)

        # Set the path
        projection._path = path

        # Return the projection
        return projection

    # -----------------------------------------------------------------

    # compatibility
    @property
    def field_x_physical(self):

        """
        This function ...
        :return:
        """

        return self.field_x

    # -----------------------------------------------------------------

    # compatibility
    @property
    def field_y_physical(self):

        """
        This function ...
        :return:
        """

        return self.field_y

    # -----------------------------------------------------------------

    @classmethod
    def prompt(cls, **properties):

        """
        This function ...
        :param properties:
        :return:
        """

        from ...core.basics.configuration import prompt_variable

        # Distance
        if properties.get("distance", None) is None: distance = prompt_variable("distance", "length_quantity", "distance")
        else: distance = properties.pop("distance")

        # Inclination
        if properties.get("inclination", None) is None: inclination = prompt_variable("inclination", "angle", "inclination angle")
        else: inclination = properties.pop("inclination")

        # Azimuth
        if properties.get("azimuth", None) is None: azimuth = prompt_variable("azimuth", "angle", "azimuth angle", "0 deg", convert_default=True)
        else: azimuth = properties.pop("azimuth")

        # Position angle
        if properties.get("position_angle", None) is None: position_angle = prompt_variable("position_angle", "angle", "position angle")
        else: position_angle = properties.get("position_angle")

        # Pixel size
        if properties.get("npixels", None) is None:
            #pixels_x = prompt_integer("pixels_x", "number of x pixels")
            #pixels_y = prompt_integer("pixels_y", "number of y pixels")
            npixels = prompt_variable("npixels", "pixelshape", "number of pixels")
        else: npixels = properties.pop("npixels")

        # Field
        if properties.get("field", None) is None:
            #field_x = prompt_variable("field_x", "length_quantity", "x field of view")
            #field_y = prompt_variable("field_y", "length_quantity", "y field of view")
            field = prompt_variable("field", "physical_extent", "field of view")
        else: field = properties.pop("field")

        # Determine physical pixelscale
        pixelscale = PhysicalPixelscale.from_shape_and_field(npixels, field)

        # Center
        if properties.get("center", None) is None:
            #center_x = prompt_variable("center_x", "length_quantity", "x center")
            #center_y = prompt_variable("center_y", "length_quantity", "y center")
            center = prompt_variable("center", "physical_or_pixel_coordinate", "center (in pixel or length coordinates)")
        else: center = properties.pop("center")
        if isinstance(center, PixelCoordinate): center = PhysicalCoordinate(center.x * pixelscale.x, center.y * pixelscale.y)

        # Create and return
        return cls(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                   pixels_x=npixels.x, pixels_y=npixels.y, center_x=center.x, center_y=center.y, field_x=field.x, field_y=field.y)

    # -----------------------------------------------------------------

    @property
    def physical_pixelscale(self):

        """
        This function
        :return:
        """

        # From field of view to pixel scale
        pixelscale_x = self.field_x_physical / self.pixels_x
        pixelscale_y = self.field_y_physical / self.pixels_y

        # Create xy extent
        return PhysicalPixelscale(pixelscale_x, pixelscale_y)

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        physical = self.physical_pixelscale
        pixelscale_x_angular = (physical.x / self.distance).to("arcsec", equivalencies=dimensionless_angles())
        pixelscale_y_angular = (physical.y / self.distance).to("arcsec", equivalencies=dimensionless_angles())

        # Create and return the pixelscale
        return Pixelscale(pixelscale_x_angular, pixelscale_y_angular)

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
        return cls(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                   pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x, center_y=center_y, field_x=field_x, field_y=field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_deprojection(cls, deprojection, distance=None, azimuth=0.0):

        """
        This function ...
        :param deprojection:
        :param distance:
        :param azimuth:
        :return:
        """

        from ...magic.basics.vector import PixelShape
        from ...magic.basics.coordinate import PixelCoordinate

        pixels_x = deprojection.x_size
        pixels_y = deprojection.y_size
        npixels = PixelShape.from_xy(pixels_x, pixels_y)

        pixelscale = deprojection.pixelscale

        # In pixel
        x_center = deprojection.x_center
        y_center = deprojection.y_center

        # To physical
        #center_x = x_center * pixelscale
        #center_y = y_center * pixelscale

        # Determine field of view
        field_x = pixelscale * pixels_x
        field_y = pixelscale * pixels_y
        field = PhysicalExtent(field_x, field_y)

        # Get physical center
        center = PixelCoordinate(x=x_center, y=y_center)
        center_physical = get_physical_center(field, npixels, center)

        inclination = deprojection.inclination
        position_angle = deprojection.position_angle

        # Get distance
        if distance is None: distance = deprojection.galaxy_distance

        # Create and return a new class instance
        return cls(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                   pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_physical.x, center_y=center_physical.y,
                   field_x=field.x, field_y=field.y)

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
        return cls(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                   pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x, center_y=center_y, field_x=field_x, field_y=field_y)

    # -----------------------------------------------------------------

    def to_wcs(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

class FaceOnProjection(GalaxyProjection):

    """
    This class ...
    """

    #def __init__(self, distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y):
    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        # inclination, azimuth, position_angle
        #super(FaceOnProjection, self).__init__(distance, 0.0, 0.0, 90., pixels_x, pixels_y, center_x, center_y, field_x, field_y)

        # Call the constructor of the base class
        super(FaceOnProjection, self).__init__()

        from ...core.units.quantity import zero_angle, right_angle

        # Checks
        if "inclination" in kwargs and kwargs.pop("inclination") != zero_angle(): raise ValueError("Inclination is not zero")
        if "azimuth" in kwargs and kwargs.pop("azimuth") != zero_angle(): raise ValueError("Azimuth is not zero")
        if "position_angle" in kwargs and kwargs.pop("position_angle") != right_angle(): raise ValueError("Position angle is not 90 degrees")

        # Set fixed properties
        self.set_fixed("inclination", value=zero_angle())
        self.set_fixed("azimuth", value=zero_angle())
        self.set_fixed("position_angle", value=right_angle())

        # Set values
        self.set_properties(kwargs)

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
        return cls(distance=distance, pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x, center_y=center_y, field_x=field_x, field_y=field_y)

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
        return cls(distance=distance, pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x, center_y=center_y, field_x=field_x, field_y=field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        # Create and return
        return cls(distance=projection.distance, pixels_x=projection.pixels_x, pixels_y=projection.pixels_y,
                   center_x=projection.center_x, center_y=projection.center_y, field_x=projection.field_x_physical,
                   field_y=projection.field_y_physical)

# -----------------------------------------------------------------

class EdgeOnProjection(GalaxyProjection):

    """
    This class ...
    """

    #def __init__(self, distance, pixels_x, pixels_y, center_x, center_y, field_x, field_y):
    def __init__(self, **kwargs):

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
        #super(EdgeOnProjection, self).__init__(distance, 90., 0., 0., pixels_x, pixels_y, center_x, center_y, field_x, field_y)

        # Call the constructor of the base class
        super(EdgeOnProjection, self).__init__()

        from ...core.units.quantity import zero_angle, right_angle

        # Checks
        if "inclination" in kwargs and kwargs.pop("inclination") != right_angle(): raise ValueError("Inclination is not 90 degrees")
        if "azimuth" in kwargs and kwargs.pop("azimuth") != zero_angle(): raise ValueError("Azimuth is not zero")
        if "position_angle" in kwargs and kwargs.pop("position_angle") != zero_angle(): raise ValueError("Position angle is not zero")

        # Set fixed properties
        self.set_fixed("inclination", value=right_angle())
        self.set_fixed("azimuth", value=zero_angle())
        self.set_fixed("position_angle", value=zero_angle())

        # Set values
        self.set_properties(kwargs)

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
        return cls(distance=distance, pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x, center_y=center_y, field_x=field_x, field_y=field_y)

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
        return cls(distance=distance, pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x, center_y=center_y, field_x=field_x, field_y=field_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection):

        """
        This function ...
        :param projection:
        :return:
        """

        # Create and return
        return cls(distance=projection.distance, pixels_x=projection.pixels_x, pixels_y=projection.pixels_y,
                   center_x=projection.center_x, center_y=projection.center_y, field_x=projection.field_x_physical,
                   field_y=projection.field_y_physical)

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

    # Return the properties
    return pixels_x, pixels_y, center_x, center_y, field_x_physical, field_y_physical

# -----------------------------------------------------------------

def get_npixels(npixels):

    """
    This function ...
    :param npixels:
    :return:
    """

    from ...magic.basics.vector import PixelShape
    from ...magic.basics.vector import IntegerExtent

    # Set npixels
    if types.is_integer_type(npixels): npixels = PixelShape.square(npixels)
    elif isinstance(npixels, PixelShape): pass
    elif isinstance(npixels, IntegerExtent): npixels = PixelShape.from_xy(npixels.x, npixels.y)
    else: raise ValueError("Don't know what to do with npixels of type " + str(type(npixels)))

    return npixels

# -----------------------------------------------------------------

def get_field(pixelscale, npixels, distance):

    """
    This function ...
    :param pixelscale:
    :param npixels:
    :param distance:
    :return:
    """

    from ...magic.basics.pixelscale import Pixelscale
    from astropy.units import Quantity

    # Get pixelscale instance
    if isinstance(pixelscale, Quantity): pixelscale = Pixelscale(pixelscale)
    elif isinstance(pixelscale, Pixelscale): pass
    else: raise ValueError("Don't know what to do with pixelscale of type " + str(type(pixelscale)))

    # Determine physical pixelscale
    phys_pixelscale = pixelscale.to_physical(distance)

    # Determine field of view
    field = PhysicalExtent(phys_pixelscale.x * npixels.x, phys_pixelscale.y * npixels.y)

    # Reutnr the field of view
    return field

# -----------------------------------------------------------------

def get_center(npixels):

    """
    This function ...
    :param npixels:
    :return:
    """

    from ...magic.basics.vector import PixelShape
    if types.is_integer_type(npixels): npixels = PixelShape.square(npixels)
    from ...magic.basics.coordinate import PixelCoordinate
    center = PixelCoordinate(x=0.5*npixels.x, y=0.5*npixels.y)
    return center

# -----------------------------------------------------------------

def get_physical_center(field, npixels, center):

    """
    Thisf unction ...
    :param field:
    :param npixels:
    :param center: pixel center
    :return:
    """

    from ...magic.basics.coordinate import PhysicalCoordinate, PixelCoordinate

    # Physical scales per pixel
    x_scale = field.x / npixels.x
    y_scale = field.y / npixels.y

    # Get center in physical coordinates
    center = PixelCoordinate(0.5 * (npixels.x-1) - center.x, 0.5 * (npixels.y-1) - center.y)

    center_x = center.x * x_scale
    center_y = center.y * y_scale

    # Return the physical coordinate
    return PhysicalCoordinate(center_x, center_y)

# -----------------------------------------------------------------
