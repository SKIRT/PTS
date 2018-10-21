#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.circle Contains the CircleRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import astronomical modules
from astropy.units import Quantity, Unit
from photutils.geometry import circular_overlap_grid
from astropy.coordinates import frame_transform_graph

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.mask import Mask
from .region import add_info, make_circle_template, coordsys_name_mapping

# -----------------------------------------------------------------

class CircleRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Set attributes
        self.center = center
        self.radius = radius

        # Call the constructor of the base class
        super(CircleRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    @property
    def unrotated_radius(self):

        """
        This function ...
        :return:
        """

        return self.radius

    # -----------------------------------------------------------------

    @property
    def circumference(self):

        """
        This function ...
        :return:
        """

        return 2. * math.pi * self.radius

    # -----------------------------------------------------------------

    @property
    def area(self):

        """
        This function ...
        :return:
        """

        return math.pi * self.radius ** 2

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :return:
        """

        return self.__class__(self.center, self.radius * value, meta=self.meta)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :return:
        """

        self.radius *= value
        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :return:
        """

        return self.__class__(self.center, self.radius / value, self.meta)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :return:
        """

        # Divide and return
        self.radius /= value
        return self

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :return:
        """

        return self.__idiv__(value)

# -----------------------------------------------------------------

class PixelCircleRegion(CircleRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be a pixel coordinate")
        if not isinstance(radius, float): raise ValueError("Radius must be a float")

        # Call the constructor of the base class
        super(PixelCircleRegion, self).__init__(center, radius, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_sky(cls, region, wcs):

        """
        This function ...
        :param region:
        :param wcs:
        :return:
        """

        # Convert the center coordinate
        center = PixelCoordinate.from_sky(region.center, wcs)

        # Convert the radius
        radius = (region.radius / wcs.average_pixelscale).to("").value

        # Create the pixel circle region
        return cls(center, radius, meta=region.meta, label=region.label, include=region.include, appearance=region.appearance)

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        return SkyCircleRegion.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size, invert=False):

        """
        This function ...
        :param x_size:
        :param y_size:
        :param invert:
        :return:
        """

        rel_center = self.center

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        fraction = circular_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, self.radius, use_exact=0, subpixels=1)

        #from ..tools import plotting
        #plotting.plot_mask(fraction)

        # Return a new mask
        #return Mask(fraction)
        mask = Mask(fraction)
        if invert: return mask.inverse()
        else: return mask

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.include: prefix = ""
        else: prefix = "-"

        coordsys = 'fk5'
        fmt = '.4f'
        radunit = 'deg'

        if radunit == 'arcsec':
            if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
            else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else: radunitstr = ''

        x = self.center.x
        y = self.center.y
        r = self.radius

        string = prefix + make_circle_template(fmt, radunitstr).format(**locals())
        string = add_info(string, self)
        return string

    # -----------------------------------------------------------------

    @property
    def axis1_min(self):
        return self.bounding_box.x_min

    # -----------------------------------------------------------------

    @property
    def axis1_max(self):
        return self.bounding_box.x_max

    # -----------------------------------------------------------------

    @property
    def axis2_min(self):
        return self.bounding_box.y_min

    # -----------------------------------------------------------------

    @property
    def axis2_max(self):
        return self.bounding_box.y_max

# -----------------------------------------------------------------

class SkyCircleRegion(CircleRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be a sky coordinate")
        if not isinstance(radius, Quantity): raise ValueError("Radius must be an angular quantity")

        # Call the constructor of the base class
        super(SkyCircleRegion, self).__init__(center, radius, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, region, wcs):

        """
        This function ...
        :param region:
        :param wcs:
        :return:
        """

        # Convert the center coordinate to sky coordinate
        center = SkyCoordinate.from_pixel(region.center, wcs)

        # Convert the radius
        radius = region.radius * wcs.average_pixelscale

        # Create a new SkyCircleRegion
        return cls(center, radius, meta=region.meta)

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        return PixelCircleRegion.from_sky(self, wcs)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.include: prefix = ""
        else: prefix = "-"

        coordsys = 'fk5'
        fmt = '.4f'
        radunit = 'deg'

        if radunit == 'arcsec':
            if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
            else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else: radunitstr = ''

        # convert coordsys string to coordsys object
        if coordsys in coordsys_name_mapping: frame = frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
        else: frame = None  # for pixel/image/physical frames

        x = float(self.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(self.center.transform_to(frame).spherical.lat.to('deg').value)
        r = float(self.radius.to(radunit).value)

        string = prefix + make_circle_template(fmt, radunitstr).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class PhysicalCircleRegion(CircleRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, **kwargs):

        """
        This function ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be a physical coordinate")
        if not isinstance(radius, Quantity): raise ValueError("Radius must be a physical quantity of length")

        # Call the constructor of the base class
        super(PhysicalCircleRegion, self).__init__(center, radius, **kwargs)

# -----------------------------------------------------------------
