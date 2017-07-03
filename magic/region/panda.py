#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.panda Contains the PandaRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate

# -----------------------------------------------------------------

class PandaRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, start_angle, stop_angle, nangle, inner, outer, nradius, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Check the angle
        #if not isinstance(angle, Angle): raise ValueError("Angle must be a Astropy Angle object")

        # Set the attributes
        self.center = center
        self.start_angle = start_angle
        self.stop_angle = stop_angle
        self.nangle = nangle
        self.inner = inner
        self.outer = outer
        self.nradius = nradius

        # Call the constructor of the base class
        super(PandaRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelPandaRegion(PandaRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, start_angle, stop_angle, nangle, inner, outer, nradius, **kwargs):

        """
        This function ...
        """

        # Check the start coordinate
        #if not isinstance(start, PixelCoordinate): raise ValueError("Start must be pixel coordinate")

        # Check the length
        #if not isinstance(length, float): raise ValueError("Length must be float")

        # Call the constructor of VectorRegion class
        PandaRegion.__init__(self, center, start_angle, stop_angle, nangle, inner, outer, nradius, **kwargs)

# -----------------------------------------------------------------

class SkyPandaRegion(PandaRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, start_angle, stop_angle, nangle, inner, outer, nradius, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        #if not isinstance(start, SkyCoordinate): raise ValueError("Start must be sky coordinate")

        # Check the length
        #if not isinstance(length, Quantity): raise ValueError("Length must be an angular quantity")

        # Call the constructor of VectorRegion class
        PandaRegion.__init__(self, center, start_angle, stop_angle, nangle, inner, outer, nradius, **kwargs)

# -----------------------------------------------------------------

class PhysicalPandaRegion(PandaRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, center, start_angle, stop_angle, nangle, inner, outer, nradius, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        #if not isinstance(start, PhysicalCoordinate): raise ValueError("Start must be physical coordinate")

        # Check the length
        #if not isinstance(length, Quantity): raise ValueError("Length must be a physical quantity of length")

        # Call the constructor of VectorRegion class
        PandaRegion.__init__(self, center, start_angle, stop_angle, nangle, inner, outer, nradius, **kwargs)

# -----------------------------------------------------------------
