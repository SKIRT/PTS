#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.rectangle Contains the RectangleRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Quantity
from photutils.geometry import elliptical_overlap_grid, circular_overlap_grid, rectangular_overlap_grid
from astropy.coordinates import frame_transform_graph
from astropy.wcs import utils

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from ..basics.mask import Mask
from .region import add_info, make_rectangle_template, coordsys_name_mapping
from ...core.units.parsing import parse_unit as u
from ..tools import coordinates

# -----------------------------------------------------------------

class RectangleRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the angle
        if angle is None: angle = Angle(0., "deg")
        elif not isinstance(angle, Angle): raise ValueError("Angle must be an Astropy Angle object")

        # Set the attributes
        self.center = center
        self.radius = radius
        self.angle = angle

        # Call the constructor of the base class
        super(RectangleRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    @property
    def width(self):
        return self.radius.axis1

    # -----------------------------------------------------------------

    @property
    def height(self):
        return self.radius.axis2

    # -----------------------------------------------------------------

    @property
    def rotated(self):
        return not bool(self.angle == Angle(0.0, "deg"))

    # -----------------------------------------------------------------

    @property
    def diagonal(self):
        return self.radius.norm

    # -----------------------------------------------------------------

    @property
    def _diagonal_angle(self):

        """
        This function ...
        :return:
        """

        return np.arctan2(self.radius.axis2, self.radius.axis1)

    # -----------------------------------------------------------------

    @property
    def corners(self):
        return [self.corner11, self.corner12, self.corner21, self.corner22]

    # -----------------------------------------------------------------

    @property
    def axis1_min(self):
        return min([corner.axis1 for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def axis1_max(self):
        return max([corner.axis1 for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def axis2_min(self):
        return min([corner.axis2 for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def axis2_max(self):
        return max([corner.axis2 for corner in self.corners])

    # -----------------------------------------------------------------

    @property
    def upper_left(self):

        """
        Upper left means lowest axis1 coordinate, highest axis2 coordinate
        :return:
        """

        most_upper_left = None

        for corner in self.corners:

            if most_upper_left is None or (corner.axis1 <= most_upper_left.axis1 and corner.axis2 >= most_upper_left.axis2):
                most_upper_left = corner

        return most_upper_left

    # -----------------------------------------------------------------

    @property
    def upper_right(self):

        """
        This function ...
        :return:
        """

        most_upper_right = None

        for corner in self.corners:

            if most_upper_right is None or (corner.axis1 >= most_upper_right.axis1 and corner.axis2 >= most_upper_right.axis2):
                most_upper_right = corner

        return most_upper_right

    # -----------------------------------------------------------------

    @property
    def lower_left(self):

        """
        This function ...
        :return:
        """

        most_lower_left = None

        for corner in self.corners:

            if most_lower_left is None or (corner.axis1 <= most_lower_left.axis1 and corner.axis2 <= most_lower_left.axis2):
                most_lower_left = corner

        return most_lower_left

    # -----------------------------------------------------------------

    @property
    def lower_right(self):

        """
        This function ...
        :return:
        """

        most_lower_right = None

        for corner in self.corners:

            if most_lower_right is None or (corner.axis1 >= most_lower_right.axis1 and corner.axis2 <= most_lower_right.axis2):
                most_lower_right = corner

        return most_lower_right

    # -----------------------------------------------------------------

    def contains(self, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        # TODO: IS THIS CORRECT???

        if coordinate.axis1 > self.axis1_max or coordinate.axis1 < self.axis1_min: return False
        if coordinate.axis2 > self.axis2_max or coordinate.axis2 < self.axis2_min: return False
        return True

    # -----------------------------------------------------------------

    @property
    def axial_ratio(self):
        return self.width / self.height

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

        self.radius *= value
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

        self.radius /= value
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

class PixelRectangleRegion(RectangleRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be a pixel coordinate")
        if not isinstance(radius, PixelStretch):
            if not isinstance(radius, float): raise ValueError("Radius must be a pixel stretch or float (for a square)")
            else: radius = PixelStretch(radius, radius)

        # Call the constructor of the base class
        super(PixelRectangleRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_sky(cls, region, wcs):

        """
        This function ...
        :param region:
        :param wcs:
        :return:
        """

        # Get center pixel coordinate
        center = region.center.to_pixel(wcs)

        ## GET THE PIXELSCALE
        result = utils.proj_plane_pixel_scales(wcs)
        # returns: A vector (ndarray) of projection plane increments corresponding to each pixel side (axis).
        # The units of the returned results are the same as the units of cdelt, crval, and cd for the celestial WCS
        # and can be obtained by inquiring the value of cunit property of the input WCS WCS object.
        x_pixelscale = result[0] * u("deg")
        y_pixelscale = result[1] * u("deg")

        #semimajor = (region.semimajor / x_pixelscale).to("").value
        #semiminor = (region.semiminor / y_pixelscale).to("").value

        width = (region.width / x_pixelscale).to("").value
        height = (region.height / y_pixelscale).to("").value

        radius = PixelStretch(width, height)

        # Convert angle
        # Set the angle
        angle = region.angle
        if angle is not None:
            try:
                orientation = wcs.standard_orientation_angle
            except ValueError:
                orientation = wcs.orientation_angle
            # Add the orientation angle (w.r.t. standard E-W and S-N projection on the x and y axes) to the position angle
            # that is expressed in the standard way
            # return self.pa + orientation
            angle = angle + orientation
        else: angle = Angle(0.0, "deg")

        # Create a new PixelRectangleRegion
        return cls(center, radius, angle, meta=region.meta, label=region.label, include=region.include, appearance=region.appearance)

    # -----------------------------------------------------------------

    @property
    def corner11(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x + self.radius.x, self.center.y - self.radius.y)
        else:
            angle_to_corner = 2.0 * np.pi - self._diagonal_angle + self.angle.to("radian").value
            x_to_corner = self.diagonal * np.cos(angle_to_corner)
            y_to_corner = self.diagonal * np.sin(angle_to_corner)
            return PixelCoordinate(self.center.x + x_to_corner, self.center.y + y_to_corner)

    # -----------------------------------------------------------------

    @property
    def corner12(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x + self.radius.x, self.center.y + self.radius.y)
        else:
            angle_to_corner = self._diagonal_angle + self.angle.to("radian").value
            x_to_corner = self.diagonal * np.cos(angle_to_corner)
            y_to_corner = self.diagonal * np.sin(angle_to_corner)
            return PixelCoordinate(self.center.x + x_to_corner, self.center.y + y_to_corner)

    # -----------------------------------------------------------------

    @property
    def corner21(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x - self.radius.x, self.center.y + self.radius.y)
        else:
            angle_to_corner = np.pi - self._diagonal_angle + self.angle.to("radian").value
            x_to_corner = self.diagonal * np.cos(angle_to_corner)
            y_to_corner = self.diagonal * np.sin(angle_to_corner)
            return PixelCoordinate(self.center.x + x_to_corner, self.center.y + y_to_corner)

    # -----------------------------------------------------------------

    @property
    def corner22(self):

        """
        This function ...
        :return:
        """

        if not self.rotated: return PixelCoordinate(self.center.x - self.radius.x, self.center.y - self.radius.y)
        else:
            angle_to_corner = np.pi + self._diagonal_angle + self.angle.to("radian").value
            x_to_corner = self.diagonal * np.cos(angle_to_corner)
            y_to_corner = self.diagonal * np.sin(angle_to_corner)
            return PixelCoordinate(self.center.x + x_to_corner, self.center.y + y_to_corner)

    # -----------------------------------------------------------------

    @property
    def unrotated_radius(self):

        """
        This function ...
        :return:
        """

        return PixelStretch(0.5 * (self.axis1_max - self.axis1_min), 0.5 * (self.axis2_max - self.axis2_min))

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        ## SIMPLE WAY

        # data = np.zeros((y_size, x_size))

        # Convert into integers
        # x_min = int(round(self.x_min))
        # x_max = int(round(self.x_max))
        # y_min = int(round(self.y_min))
        # y_max = int(round(self.y_max))

        # data[y_min:y_max, x_min:x_max] = 1

        # Return a new Mask object
        # return Mask(data)

        ## OTHER WAY

        rel_center = self.center

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        width = 2. * self.radius.x
        height = 2. * self.radius.y

        fraction = rectangular_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, width, height, self.angle.to("radian").value, 0, 1)

        # Return the mask
        return Mask(fraction)

    # -----------------------------------------------------------------

    def contains(self, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        if not self.rotated: return self.x_min <= coordinate.x <= self.x_max and self.y_min <= coordinate.y <= self.y_max
        else:

            # NOTE: does not work properly yet!

            # Source: http://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle
            # AM = position - self.upper_left
            # AB = self.upper_right - self.upper_left
            # AD = self.lower_left - self.upper_left
            am = coordinate - self.upper_left
            ab = self.upper_right - self.upper_left
            ad = self.lower_left - self.upper_left

            return (0 < am.dot(ab) < ab.dot(ab)) and (0 < am.dot(ad) < ad.dot(ad))

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
        d1 = 2.0 * self.radius.x
        d2 = 2.0 * self.radius.y
        ang = self.angle.to("deg").value

        string = prefix + make_rectangle_template(fmt, radunitstr).format(**locals())
        string = add_info(string, self)
        return string

    # -----------------------------------------------------------------

    def to_mpl_patch(self):

        """
        This function ...
        :return:
        """

        from matplotlib.patches import Rectangle as plt_Rectangle
        raise NotImplementedError("Not yet implemented")

# -----------------------------------------------------------------

class SkyRectangleRegion(RectangleRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be a sky coordinate")
        if not isinstance(radius, SkyStretch):
            if not isinstance(radius, Quantity): raise ValueError("Radius must be a sky stretch or an angular quantity (for a square)")
            else: radius = SkyStretch(radius, radius)

        # Call the constructor of the base class
        super(SkyRectangleRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def corner11(self):

        """
        This function ...
        :return:
        """

        if not self.rotated:
            ra_min, ra_max = coordinates.ra_around(self.center.ra, self.radius.ra, self.center.dec)
            return SkyCoordinate(ra_max, self.center.dec - self.radius.dec)
        else:
            angle_to_corner = 2.0 * np.pi - self._diagonal_angle + self.angle.to("radian").value
            ra_to_corner = self.diagonal * np.cos(angle_to_corner)
            dec_to_corner = self.diagonal * np.sin(angle_to_corner)
            ra_min, ra_max = coordinates.ra_around(self.center.ra, ra_to_corner, self.center.dec)
            return SkyCoordinate(ra_max, self.center.dec + dec_to_corner)

    # -----------------------------------------------------------------

    @property
    def corner12(self):

        """
        This function ...
        :return:
        """

        if not self.rotated:
            ra_min, ra_max = coordinates.ra_around(self.center.ra, self.radius.ra, self.center.dec)
            return SkyCoordinate(ra_max, self.center.dec + self.radius.dec)
        else:
            angle_to_corner = self._diagonal_angle + self.angle.to("radian").value
            ra_to_corner = self.diagonal * np.cos(angle_to_corner)
            dec_to_corner = self.diagonal * np.sin(angle_to_corner)
            ra_min, ra_max = coordinates.ra_around(self.center.ra, ra_to_corner, self.center.dec)
            return SkyCoordinate(ra_max, dec_to_corner)

    # -----------------------------------------------------------------

    @property
    def corner21(self):

        """
        This function ...
        :return:
        """

        if not self.rotated:
            ra_min, ra_max = coordinates.ra_around(self.center.ra, self.radius.ra, self.center.dec)
            return SkyCoordinate(ra_min, self.center.dec + self.radius.dec)
        else:
            angle_to_corner = np.pi - self._diagonal_angle + self.angle.to("radian").value
            ra_to_corner = self.diagonal * np.cos(angle_to_corner)
            dec_to_corner = self.diagonal * np.sin(angle_to_corner)
            ra_min, ra_max = coordinates.ra_around(self.center.ra, ra_to_corner, self.center.dec)
            return SkyCoordinate(ra_min, self.center.dec + dec_to_corner)

    # -----------------------------------------------------------------

    @property
    def corner22(self):

        """
        This function ...
        :return:
        """

        if not self.rotated:
            ra_min, ra_max = coordinates.ra_around(self.center.ra, self.radius.ra, self.center.dec)
            return SkyCoordinate(ra_min, self.center.dec - self.radius.dec)
        else:
            angle_to_corner = np.pi + self._diagonal_angle + self.angle.to("radian").value
            ra_to_corner = self.diagonal * np.cos(angle_to_corner)
            dec_to_corner = self.diagonal * np.sin(angle_to_corner)
            ra_min, ra_max = coordinates.ra_around(self.center.ra, ra_to_corner, self.center.dec)
            return SkyCoordinate(ra_min, self.center.dec + dec_to_corner)

    # -----------------------------------------------------------------

    @property
    def unrotated_radius(self):

        """
        This function ...
        :return:
        """

        return SkyStretch(0.5 * (self.axis1_max - self.axis1_min), 0.5 * (self.axis2_max - self.axis2_min))

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
        d1 = 2.0 * float(self.radius.ra.to(radunit).value)
        d2 = 2.0 * float(self.radius.dec.to(radunit).value)
        ang = float(self.angle.to('deg').value)

        string = prefix + make_rectangle_template(fmt, radunitstr).format(**locals())
        string = add_info(string, self)
        return string

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Create a pixel ellipse region
        return PixelRectangleRegion.from_sky(self, wcs)

# -----------------------------------------------------------------

class PhysicalRectangleRegion(RectangleRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be a physical coordinate")
        if not isinstance(radius, PhysicalStretch):
            if not isinstance(radius, Quantity): raise ValueError("Radius must be a physical stretch or a quantity of length (for a square)")
            else: radius = PhysicalStretch(radius, radius)

        # Call the constructor of the base class
        super(PhysicalRectangleRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def corner11(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    @property
    def corner12(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    @property
    def corner21(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    @property
    def corner22(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

# -----------------------------------------------------------------
