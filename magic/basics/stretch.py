#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.stretch Contains the Stretch class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .vector import Extent

# -----------------------------------------------------------------

class SkyExtent(object):

    """
    This class ...
    """

    def __init__(self, ra, dec):

        """
        The constructor ...
        :param ra:
        :param dec:
        """

        self.ra = ra
        self.dec = dec

# -----------------------------------------------------------------

class PhysicalExtent(object):

    """
    This class ...
    """

    def __init__(self, length1, length2):

        """
        The constructor ...
        :param length1:
        :param length2:
        """

        self.length1 = length1
        self.length2 = length2

# -----------------------------------------------------------------

class Stretch(object):

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

class PixelStretch(Extent, Stretch):

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

        # Call the constructor of the base classes
        Extent.__init__(self, x, y)
        Stretch.__init__(self, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_sky(cls, sky_stretch, wcs, mode="wcs"):

        """
        This function ...
        :param sky_stretch:
        :param wcs:
        :param mode:
        :return:
        """

        x_pixels = (sky_stretch.ra / wcs.pixelscale.x.to("deg/pix")).to("pix").value
        y_pixels = (sky_stretch.dec / wcs.pixelscale.y.to("deg/pix")).to("pix").value

        return cls(x_pixels, y_pixels)

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Return a new SkyCoordinate
        return SkyStretch.from_pixel(self, wcs)

# -----------------------------------------------------------------

class SkyStretch(SkyExtent, Stretch):

    """
    This class ...
    """

    def __init__(self, ra, dec, **kwargs):

        """
        The constructor ...
        :param ra, dec:
        :param kwargs:
        :return:
        """

        # Call the constructors of the base classes
        SkyExtent.__init__(self, ra, dec)
        Stretch.__init__(self, **kwargs)

    # -----------------------------------------------------------------

    def to_pixel(self, wcs, mode='wcs'):

        """
        This function ...
        :param wcs:
        :param mode:
        :return:
        """

        return PixelStretch.from_sky(self, wcs, mode)

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, stretch, wcs, mode='wcs'):

        """
        This function ...
        :param stretch:
        :param wcs:
        :param mode:
        :return:
        """

        x_angular = wcs.pixelscale.x.to("deg/pix") * stretch.x * Unit("pix")
        y_angular = wcs.pixelscale.y.to("deg/pix") * stretch.y * Unit("pix")

        return cls(x_angular, y_angular)

# -----------------------------------------------------------------

class PhysicalStretch(PhysicalExtent, Stretch):

    """
    This class ...
    """

    def __init__(self, length1, length2, **kwargs):

        """
        The constructor ...
        :param length1:
        :param length2:
        :param kwargs:
        """

        # Call the constructors of the base classes
        PhysicalExtent.__init__(self, length1, length2)
        Stretch.__init__(self, **kwargs)

# -----------------------------------------------------------------
