#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.point Contains the PointRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .rectangle import PixelRectangleRegion, SkyRectangleRegion
from ..basics.mask import Mask

# -----------------------------------------------------------------

class PointRegion(Region):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param axis1: value along axis 1
        :param axis2: value along axis 2
        :param kwargs:
        """

        # Call the constructor of the base class
        #super(PointRegion, self).__init__(**kwargs) # DOESN'T WORK: TypeError: __init__() takes exactly 3 arguments (1 given) ??
        Region.__init__(self, **kwargs)

# -----------------------------------------------------------------

class PixelPointRegion(PointRegion, PixelCoordinate, PixelRegion):

    """
    This class ...
    """

    def __init__(self, x, y, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        """

        # Call the constructor of PointRegion class
        PointRegion.__init__(self, **kwargs)

        # Call the constructor of the PixelCoordinate class
        PixelCoordinate.__init__(self, x, y, **kwargs)

    # -----------------------------------------------------------------

    @property
    def axis1(self):

        """
        This function ...
        :return:
        """

        return self.x

    # -----------------------------------------------------------------

    @axis1.setter
    def axis1(self, value):

        """
        This function ...
        :return:
        """

        self.x = value

    # -----------------------------------------------------------------

    @property
    def axis2(self):

        """
        This function ...
        :return:
        """

        return self.y

    # -----------------------------------------------------------------

    @axis2.setter
    def axis2(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.y = value

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        return Mask.empty(x_size, y_size)

    # -----------------------------------------------------------------

    def to_mpl_patch(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def to_region_string(self, coordinate_system=True):

        """
        This function ...
        :return:
        """

        # Create the suffix
        if len(self.meta) > 0:
            suffix = " #"
            for key in self.meta:
                if key == "text": suffix += " " + key + " = {" + str(self.meta[key]) + "}"
                else: suffix += " " + key + " = " + str(self.meta[key])
        else: suffix = ""

        # Create and return the line
        if coordinate_system: line = "image;point({},{})".format(self.x+1, self.y+1) + suffix
        else: line = "point({},{})".format(self.x+1, self.y+1) + suffix
        return line

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        return PixelRectangleRegion(PixelCoordinate(self.x, self.y), PixelStretch(1.0, 1.0))

    # -----------------------------------------------------------------

    def __add__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return self.__class__(self.x + extent.x, self.y + extent.y, self.meta)

    # -----------------------------------------------------------------

    def __iadd__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        # Add x and y
        self.x += extent.x
        self.y += extent.y
        return self

    # -----------------------------------------------------------------

    def __sub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        return self.__class__(self.x - extent.x, self.y - extent.y, self.meta)

    # -----------------------------------------------------------------

    def __isub__(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        self.x -= extent.x
        self.y -= extent.y

        return self

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__class__(self.x, self.y, self.meta)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__class__(self.x, self.y, self.meta)


    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__idiv__(value)

# -----------------------------------------------------------------

class SkyPointRegion(PointRegion, SkyCoordinate, SkyRegion):

    """
    This class ...
    """

    def __init__(self, ra, dec, **kwargs):

        """
        This function ...
        :param ra:
        :param dec:
        :param kwargs:
        """

        # Call the constructor of PointRegion class
        PointRegion.__init__(self, **kwargs)

        # Call the constructor of the PixelCoordinate class
        SkyCoordinate.__init__(ra, dec, **kwargs)

    # -----------------------------------------------------------------

    @property
    def axis1(self):

        """
        This function ...
        :return:
        """

        return self.ra

    # -----------------------------------------------------------------

    @axis1.setter
    def axis1(self, value):

        """
        This function ...
        :return:
        """

        self.ra = value

    # -----------------------------------------------------------------

    @property
    def axis2(self):

        """
        This function ...
        :return:
        """

        return self.dec

    # -----------------------------------------------------------------

    @axis2.setter
    def axis2(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.dec = value

# -----------------------------------------------------------------

class PhysicalPointRegion(PointRegion, PhysicalCoordinate):

    """
    This class ...
    """

    def __init__(self, axis1, axis2, **kwargs):

        """
        This function ...
        :param axis1:
        :param axis2:
        :param kwargs:
        """

        # Call the constructor of PointRegion class
        PointRegion.__init__(self, **kwargs)

        # Call the constructor of the PixelCoordinate class
        PhysicalCoordinate.__init__(axis1, axis2, **kwargs)

# -----------------------------------------------------------------
