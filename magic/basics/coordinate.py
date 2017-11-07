#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.coordinate Contains the Coordinate class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import SkyCoord
from astropy.units import dimensionless_angles, Unit

# Import the relevant PTS classes and modules
from .vector import Position, Position3D

# -----------------------------------------------------------------

class Coordinate(object):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """
            
        # Set the meta information
        self.meta = kwargs.pop("meta", dict())
        
# -----------------------------------------------------------------

class PixelCoordinate(Position, Coordinate):

    """
    This class ...
    """

    def __init__(self, x, y, **kwargs):

        """
        The constructor ...
        :param x:
        :param y:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        Position.__init__(self, x, y)
        Coordinate.__init__(self, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_sky(cls, coordinate, wcs, mode='wcs'):

        """
        This function ...
        :param coordinate:
        :param wcs:
        :param mode:
        :return:
        """

        standard = SkyCoord(ra=coordinate.ra, dec=coordinate.dec)
        x, y = standard.to_pixel(wcs, origin=0, mode=mode)
        #x, y = super(SkyCoordinate, coordinate).to_pixel(wcs, origin=0, mode=mode)
        return cls(float(x), float(y), meta=coordinate.meta)

    # -----------------------------------------------------------------

    @property
    def axis1(self):

        """
        This property ...
        :return:
        """

        return self.x

    # -----------------------------------------------------------------

    @property
    def axis2(self):

        """
        This function ...
        :return:
        """

        return self.y

    # -----------------------------------------------------------------

    def to_sky(self, wcs, mode='wcs'):

        """
        This function ...
        :param wcs:
        :param mode:
        :return:
        """

        # Return a new SkyCoordinate
        return SkyCoordinate.from_pixel(self, wcs, mode)

# -----------------------------------------------------------------

class SkyCoordinate(SkyCoord, Coordinate):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the Coordinate base class
        Coordinate.__init__(self, **kwargs)

        # Remove keyword arguments that are not recognized by the SkyCoord class
        allowed_keys = ["frame", "unit", "obstime", "equinox", "representation", "copy", "ra", "dec", "l", "b", "x", "y", "z", "w", "u", "v"]
        to_be_removed_keys = []
        for key in kwargs:
            if key not in allowed_keys: to_be_removed_keys.append(key)
        for key in to_be_removed_keys: del kwargs[key]

        # Call the constructor of the base class
        SkyCoord.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_astropy(cls, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        return cls(ra=coordinate.ra, dec=coordinate.dec)

    # -----------------------------------------------------------------

    @property
    def axis1(self):

        """
        This property ...
        :return:
        """

        return self.ra

    # -----------------------------------------------------------------

    @property
    def axis2(self):

        """
        This function ...
        :return:
        """

        return self.dec

    # -----------------------------------------------------------------

    def to_pixel(self, wcs, mode='wcs'):

        """
        This function ...
        :param wcs
        :param mode:
        :return:
        """

        # Make pixel coordinate
        return PixelCoordinate.from_sky(self, wcs, mode)

    # -----------------------------------------------------------------

    @property
    def ra_min(self):

        """
        This function ...
        :return:
        """

        return self.ra

    # -----------------------------------------------------------------

    @property
    def ra_max(self):

        """
        This function ...
        :return:
        """

        return self.ra

    # -----------------------------------------------------------------

    @property
    def dec_min(self):

        """
        This function ...
        :return:
        """

        return self.dec

    # -----------------------------------------------------------------

    @property
    def dec_max(self):

        """
        This function ...
        :return:
        """

        return self.dec

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, coordinate, wcs, mode='wcs'):

        """
        This function ...
        :param coordinate:
        :param wcs:
        :param mode:
        :return:
        """

        skycoordinate = super(SkyCoordinate, cls).from_pixel(coordinate.x, coordinate.y, wcs, origin=0, mode=mode)
        return cls(ra=skycoordinate.ra.deg, dec=skycoordinate.dec.deg, unit="deg")

    # -----------------------------------------------------------------

    def to_physical(self, wcs, distance):

        """
        This function ...
        :param wcs:
        :param distance:
        :return:
        """

        return PhysicalCoordinate.from_sky(self, wcs, distance)

    # -----------------------------------------------------------------

    def to_astropy(self):

        """
        This function ...
        :return:
        """

        return SkyCoord(ra=self.ra.to("deg").value, dec=self.dec.to("deg").value, unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :param coordinate_system:
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Get the RA and DEC
        ra_deg = self.ra.to("deg").value
        dec_deg = self.dec.to("deg").value

        # Create and return the line
        if coordinate_system: line = "fk5;point({},{})".format(ra_deg, dec_deg) + suffix
        else: line = "point({},{})".format(ra_deg, dec_deg) + suffix
        return line

# -----------------------------------------------------------------

class PhysicalCoordinate(Position, Coordinate):

    """
    This class ...
    """

    def __init__(self, axis1, axis2, **kwargs):

        """
        The constructor ...
        :param axis1:
        :param axis2:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base classes
        Position.__init__(self, axis1, axis2)
        Coordinate.__init__(self, **kwargs)

        # Check whether of right type
        from ...core.tools import types
        if not types.is_length_quantity(self.x): raise ValueError("Arguments must be length quantities")
        if not types.is_length_quantity(self.y): raise ValueError("Arguments must be length quantities")

    # -----------------------------------------------------------------

    @classmethod
    def zero(cls, unit="pc"):

        """
        Thisj function ...
        :return:
        """

        # Create zero length
        zero_length = 0. * Unit(unit)

        # Create and return
        return cls(zero_length, zero_length)

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, coordinate, wcs, distance, from_center=False):

        """
        This function ...
        :param coordinate:
        :param wcs:
        :param distance:
        :param from_center:
        :return:
        """

        if from_center:
            # PIXEL SIZE
            pixels_x = wcs.xsize
            pixels_y = wcs.ysize
            center = Position(0.5 * pixels_x - coordinate.x - 0.5, 0.5 * pixels_y - coordinate.y - 0.5)
            center_x = center.x
            center_y = center.y
        else:
            center_x = coordinate.x
            center_y = coordinate.y

        # Calculate x and y coordinate in physical length unit
        center_x = (center_x * wcs.pixelscale.x.to("deg") * distance).to("pc", equivalencies=dimensionless_angles())
        center_y = (center_y * wcs.pixelscale.y.to("deg") * distance).to("pc", equivalencies=dimensionless_angles())

        # Create and return the coordinate
        return cls(center_x, center_y)

    # -----------------------------------------------------------------

    @classmethod
    def from_sky(cls, coordinate, wcs, distance, from_center=False, mode="wcs"):

        """
        This function ...
        :param coordinate:
        :param wcs:
        :param distance:
        :param from_center:
        :param mode:
        :return:
        """

        # Convert to pixel coordinate
        pixel = coordinate.to_pixel(wcs, mode=mode)

        # From pixel coordinate
        return cls.from_pixel(pixel, wcs, distance, from_center=from_center)

    # -----------------------------------------------------------------

    @property
    def axis1(self):

        """
        This property ...
        :return:
        """

        return self.x

    # -----------------------------------------------------------------

    @property
    def axis2(self):

        """
        This function ...
        :return:
        """

        return self.y

# -----------------------------------------------------------------

class PhysicalCoordinate3D(Position3D, Coordinate):

    """
    This class ...
    """

    def __init__(self, axis1, axis2, axis3, **kwargs):

        """
        The constructor ...
        :param axis1:
        :param axis2:
        :param axis3:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base classes
        Position3D.__init__(self, axis1, axis2, axis3)
        Coordinate.__init__(self, **kwargs)

        # Check whether of right type
        from ...core.tools import types
        if not types.is_length_quantity(self.x): raise ValueError("Arguments must be length quantities")
        if not types.is_length_quantity(self.y): raise ValueError("Arguments must be length quantities")
        if not types.is_length_quantity(self.z): raise ValueError("Arguments must be length quantities")

# -----------------------------------------------------------------
