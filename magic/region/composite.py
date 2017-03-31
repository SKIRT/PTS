#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.composite Contains the CompositeRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion

# -----------------------------------------------------------------

class CompositeRegion(Region):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Set the region elements
        self.elements = args

        # Get the angle
        self.angle = kwargs.pop("angle", Angle(0.0, "deg"))

        # Check the angle
        if not isinstance(self.angle, Angle): raise ValueError("Angle must be an Astropy Angle object")

        # Call the constructor of the base class
        super(CompositeRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.elements)

    # -----------------------------------------------------------------

    @property
    def axis1_center(self):

        """
        This property ...
        :return:
        """

        return sum([element.center.axis1 for element in self.elements]) / float(len(self))

    # -----------------------------------------------------------------

    @property
    def axis2_center(self):

        """
        This property ...
        :return:
        """

        return sum([element.center.axis2 for element in self.elements]) / float(len(self))

# -----------------------------------------------------------------

class PixelCompositeRegion(CompositeRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check that all elements are pixel regions
        for arg in args:
            if not isinstance(arg, PixelRegion): raise ValueError("All regions of the composite should be in pixel coordinates")

        # Call the constructor of the CompositeRegion base class
        CompositeRegion.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return PixelCoordinate(self.axis1_center, self.axis2_center)

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        from ..core.mask import Mask
        masks = [element.to_mask(x_size, y_size) for element in self.elements]
        return Mask.union(*masks)

# -----------------------------------------------------------------

class SkyCompositeRegion(CompositeRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check that all elements are sky regions
        for arg in args:
            if not isinstance(arg, SkyRegion): raise ValueError("All regions of the composite should be in sky coordinates")

        # Call the constructor of the CompositeRegion base class
        CompositeRegion.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return SkyCoordinate(self.axis1_center, self.axis2_center)

# -----------------------------------------------------------------

class PhysicalCompositeRegion(CompositeRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check that all elements are physical regions
        for arg in args:
            if not isinstance(arg, PhysicalRegion): raise ValueError("All regions of the composite should be in physical coordinates")

        # Call the constructor of the CompositeRegion base class
        CompositeRegion.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return PhysicalCoordinate(self.axis1_center, self.axis2_center)

# -----------------------------------------------------------------
