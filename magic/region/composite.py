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

    @property
    def axis1_min(self):

        """
        This function ...
        :return: 
        """

        the_min = None
        for element in self.elements:
            if the_min is None or element.axis1_min < the_min: the_min = element.axis1_min
        return the_min

    # -----------------------------------------------------------------

    @property
    def axis1_max(self):

        """
        This function ...
        :return: 
        """

        the_max = None
        for element in self.elements:
            if the_max is None or element.axis1_max > the_max: the_max = element.axis1_max
        return the_max

    # -----------------------------------------------------------------

    @property
    def axis2_min(self):

        """
        This function ...
        :return: 
        """

        the_min = None
        for element in self.elements:
            if the_min is None or element.axis2_min < the_min: the_min = element.axis2_min
        return the_min

    # -----------------------------------------------------------------

    @property
    def axis2_max(self):

        """
        This function ...
        :return: 
        """

        the_max = None
        for element in self.elements:
            if the_max is None or element.axis2_max > the_max: the_max = element.axis2_max
        return the_max

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new = self.copy()
        new *= value
        return new

    # -----------------------------------------------------------------

    def __rmul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new = self.copy()
        new *= value
        return new

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for element in self.elements: element *= value
        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new = self.copy()
        new /= value
        return new

    # -----------------------------------------------------------------

    def __rdiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        raise SyntaxError("You cannot divide by a region")

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for element in self.elements: element /= value
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

    def __rtruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        raise SyntaxError("You cannot divide by a region")

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__idiv__(value)

    # -----------------------------------------------------------------

    def to_mpl_patches(self):

        """
        This function ...
        :return:
        """

        patches = []
        for element in self.elements:
            patch = element.to_mpl_patch()
            patches.append(patch)
        return patches

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

        include = []
        exclude = []

        for element in self.elements:

            if element.include: include.append(element.to_mask(x_size, y_size))
            else: exclude.append(element.to_mask(x_size, y_size, invert=True))

        if len(exclude) > 0: return Mask.intersection(Mask.union(*include), Mask.union(*exclude))
        else: return Mask.union(*include)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        prefix = "" if self.include else "-"

        # Get properties
        x = self.center.x
        y = self.center.y
        ang = self.angle

        # Create string
        composite_string = prefix + ds9_strings['composite'].format(**locals())
        composite_string = add_info(composite_string, composite)

        # Add the strings for the composite elements
        output = composite_string + "\n" + " ||\n".join([regular_to_string(element, ds9_strings, frame, radunit, fmt) for element in composite.elements])

        # Return the string
        return output

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

    def __str__(self):

        """
        This function ...
        :return:
        """

        prefix = "" if composite.include else "-"

        # Get properties
        x = float(composite.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(composite.center.transform_to(frame).spherical.lat.to('deg').value)
        ang = composite.angle

        # Create string
        composite_string = prefix + ds9_strings['composite'].format(**locals())
        composite_string = add_info(composite_string, composite)

        # Add the strings for the composite elements
        output = composite_string + "\n" + " ||\n".join(
            [regular_to_string(element, ds9_strings, frame, radunit, fmt) for element in composite.elements])

        # Return the string
        return output

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
