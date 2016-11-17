#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.region Contains the (abstract) Region class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy

# -----------------------------------------------------------------

class Region(object):

    """
    This class ...
    """

    default_extension = "reg"

    # -----------------------------------------------------------------

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Set the label
        self.label = kwargs.pop("label", None)

        # Set the 'exclude' flag
        self.include = kwargs.pop("include", True)

        # Set the appearance info
        self.appearance = kwargs.pop("appearance", dict())

        # Set the meta information
        self.meta = kwargs.pop("meta", dict())

    # -----------------------------------------------------------------

    @property
    def has_label(self):

        """
        This function ...
        :return:
        """

        return self.label is not None

    # -----------------------------------------------------------------

    @property
    def has_appearance(self):

        """
        This function ...
        :return:
        """

        return len(self.appearance) > 0

    # -----------------------------------------------------------------

    @property
    def has_meta(self):

        """
        This function ...
        :return:
        """

        return len(self.meta) > 0

    # -----------------------------------------------------------------

    @property
    def has_info(self):

        """
        This property ...
        :return:
        """

        return self.has_label or self.has_appearance or self.has_meta

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

# -----------------------------------------------------------------

class PixelRegion(Region):

    """
    This class ...
    """

    def __add__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        copy = self.copy()
        copy += other
        return copy

    # -----------------------------------------------------------------

    def __iadd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        self.center += other
        return self

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        copy = self.copy()
        copy -= other
        return copy

    # -----------------------------------------------------------------

    def __isub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        self.center -= other
        return self

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This property ...
        :return:
        """

        return self.axis1_min

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This function ...
        :return:
        """

        return self.axis1_max

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This function ...
        :return:
        """

        return self.axis2_min

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This function ...
        :return:
        """

        return self.axis2_max

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Import the relevant PTS classes and modules
        from .rectangle import PixelRectangleRegion

        # Create the rectangle region and return it
        return PixelRectangleRegion(self.center, self.unrotated_radius)

# -----------------------------------------------------------------

class SkyRegion(Region):

    """
    This class ...
    """

    @property
    def ra_min(self):

        """
        This property ...
        :return:
        """

        return self.axis1_min

    # -----------------------------------------------------------------

    @property
    def ra_max(self):

        """
        This property ...
        :return:
        """

        return self.axis1_max

    # -----------------------------------------------------------------

    @property
    def dec_min(self):

        """
        This property ...
        :return:
        """

        return self.axis2_min

    # -----------------------------------------------------------------

    @property
    def dec_max(self):

        """
        This function ...
        :return:
        """

        return self.axis2_max

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Import the relevant PTS classes and modules
        from .rectangle import SkyRectangleRegion

        # Create the rectangle region and return it
        return SkyRectangleRegion(self.center, self.radius)

# -----------------------------------------------------------------

class PhysicalRegion(Region):

    """
    This class ...
    """

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Import the relevant PTS classes and modules
        from .rectangle import PhysicalRectangleRegion

        # Create the rectangle region and return it
        return PhysicalRectangleRegion(self.center, self.radius)

# -----------------------------------------------------------------
