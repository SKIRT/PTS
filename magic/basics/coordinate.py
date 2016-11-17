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

# Import the relevant PTS classes and modules
from .vector import Position, Extent

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

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Return a new SkyCoordinate
        return SkyCoordinate.from_pixel(self, wcs)

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

        # Call the constructor of the base class
        SkyCoord.__init__(self, *args, **kwargs)

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

        x, y = super(SkyCoordinate, self).to_pixel(wcs, origin=0, mode=mode)
        return PixelCoordinate(float(x), float(y), meta=self.meta)

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
        return cls(ra=skycoordinate.ra.deg, dec=skycoordinate.dec.deg, unit="deg", meta=coordinate.meta)

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

class PhysicalCoordinate(Coordinate):

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

        self.axis1 = axis1
        self.axis2 = axis2

        # Call the Coordinate constructor
        Coordinate.__init__(self, **kwargs)

# -----------------------------------------------------------------
