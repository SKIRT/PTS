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

# Import standard modules
import math
import warnings

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.wcs import utils
from photutils.geometry import elliptical_overlap_grid, circular_overlap_grid, rectangular_overlap_grid
from astropy.coordinates import frame_transform_graph

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from ...core.units.parsing import parse_unit as u
from ..core.mask import Mask
from .region import add_info, make_ellipse_template, coordsys_name_mapping
from ..tools import coordinates

# -----------------------------------------------------------------

class EllipseRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :param kwargs:
        """


        # Check the angle
        if angle is None: angle = Angle(0., "deg")
        elif angle == 0: angle = Angle(0., "deg")
        elif not isinstance(angle, Angle): raise ValueError("Angle must be an Astropy Angle object")

        # Check whether axis1 > axis2
        if radius.axis1 < radius.axis2:
            #raise ValueError("Semimajor axis length (axis1) must be larger than semiminor axis length (axis2). Consider rotation.")
            axis2 = radius.axis2
            radius.axis2 = radius.axis1
            radius.axis1 = axis2
            angle += Angle(90, "deg")

        # Set attributes
        self.center = center
        self.radius = radius
        self.angle = angle

        # Call the constructor of the base class
        super(EllipseRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    @property
    def major(self):

        """
        This function ...
        :return:
        """

        return 2.0 * self.radius.axis1

    # -----------------------------------------------------------------

    @property
    def semimajor(self):

        """
        This function ...
        :return:
        """

        return self.radius.axis1

    # -----------------------------------------------------------------

    @semimajor.setter
    def semimajor(self, value):

        """
        This function ...
        :return:
        """

        self.radius.axis1 = value

    # -----------------------------------------------------------------

    @property
    def minor(self):

        """
        This function ...
        :return:
        """

        return 2.0 * self.radius.axis2

    # -----------------------------------------------------------------

    @property
    def semiminor(self):

        """
        This function ...
        :return:
        """

        return self.radius.axis2

    # -----------------------------------------------------------------

    @semiminor.setter
    def semiminor(self, value):

        """
        This function ...
        :return:
        """

        self.radius.axis2 = value

    # -----------------------------------------------------------------

    @property
    def axial_ratio(self):

        """
        This function ...
        :return:
        """

        return self.semiminor / self.semimajor

    # -----------------------------------------------------------------

    @property
    def ellipticity(self):

        """
        This function ...
        :return:
        """

        return 1. - self.axial_ratio

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

class PixelEllipseRegion(EllipseRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be a pixel coordinate")
        if not isinstance(radius, PixelStretch): raise ValueError("Radius must be a pixel stretch")

        # Call the constructor of the base class
        super(PixelEllipseRegion, self).__init__(center, radius, angle, **kwargs)

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
        # pixelscale = Extent(x_pixelscale, y_pixelscale)

        #print(region.semimajor)
        #print(region.semiminor)

        #print(x_pixelscale)
        #print(y_pixelscale)

        semimajor = (region.semimajor / x_pixelscale).to("").value
        semiminor = (region.semiminor / y_pixelscale).to("").value

        radius = PixelStretch(semimajor, semiminor)

        # Convert angle
        # Set the angle
        angle = region.angle
        if angle is not None:
            try: orientation = wcs.standard_orientation_angle
            except ValueError: orientation = wcs.orientation_angle
            # Add the orientation angle (w.r.t. standard E-W and S-N projection on the x and y axes) to the position angle
            # that is expressed in the standard way
            #return self.pa + orientation
            angle = angle + orientation
        else: angle = Angle(0.0, "deg")

        # Create a new PixelEllipse
        return cls(center, radius, angle, meta=region.meta, label=region.label, include=region.include, appearance=region.appearance)

    # -----------------------------------------------------------------

    @property
    def unrotated_radius(self):

        """
        This function ...
        :return:
        """

        x_radius = self.radius.x
        y_radius = self.radius.y

        a_projected_x = x_radius * math.cos(self.angle.radian)
        b_projected_x = y_radius * math.sin(self.angle.radian)
        a_projected_y = x_radius * math.sin(self.angle.radian)
        b_projected_y = y_radius * math.cos(self.angle.radian)

        box_x_radius = max(abs(a_projected_x), abs(b_projected_x))
        box_y_radius = max(abs(a_projected_y), abs(b_projected_y))

        radius = PixelStretch(box_x_radius, box_y_radius)

        # Return the radius
        return radius

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

    def to_mask(self, x_size, y_size, invert=False):

        """
        This function ...
        :param x_size:
        :param y_size:
        :param invert:
        :return:
        """

        rel_center = self.center

        a = self.radius.x
        b = self.radius.y

        # theta in radians !
        theta = self.angle.radian

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        # Calculate the mask
        fraction = elliptical_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, a, b, theta, use_exact=0, subpixels=1)

        # xmin, xmax, ymin, ymax : float
        #    Extent of the grid in the x and y direction.
        # nx, ny : int
        #    Grid dimensions.
        # rx : float
        #    The semimajor axis of the ellipse.
        # ry : float
        #    The semiminor axis of the ellipse.
        # theta : float
        #    The position angle of the semimajor axis in radians (counterclockwise).
        # use_exact : 0 or 1
        #    If set to 1, calculates the exact overlap, while if set to 0, uses a
        #    subpixel sampling method with ``subpixel`` subpixels in each direction.
        # subpixels : int
        #    If ``use_exact`` is 0, each pixel is resampled by this factor in each
        #    dimension. Thus, each pixel is divided into ``subpixels ** 2``
        #    subpixels.

        mask = Mask(fraction)

        # Return
        if invert: return mask.inverse()
        else: return mask

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        return SkyEllipseRegion.from_pixel(self, wcs)

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
        r1 = self.semimajor
        r2 = self.semiminor
        ang = self.angle.to("deg").value

        string = prefix + make_ellipse_template(fmt, radunitstr).format(**locals())
        string = add_info(string, self)
        return string

    # -----------------------------------------------------------------

    def to_mpl_patch(self):

        """
        This function ...
        :return:
        """

        from matplotlib.patches import Ellipse as plt_Ellipse
        return plt_Ellipse((self.center.x, self.center.y), 2.0 * self.radius.x, 2.0 * self.radius.y, self.angle.to("deg").value, edgecolor="red", facecolor="none", lw=5)

# -----------------------------------------------------------------

class SkyEllipseRegion(EllipseRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be a sky coordinate")
        if not isinstance(radius, SkyStretch): raise ValueError("Radius must be a sky stretch")

        # Call the constructor of the base class
        super(SkyEllipseRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, region, wcs):

        """
        This function ...
        :param region:
        :param wcs:
        :return:
        """

        # Get center sky coordinate
        center = region.center.to_sky(wcs)

        ## GET THE PIXELSCALE
        result = utils.proj_plane_pixel_scales(wcs)
        # returns: A vector (ndarray) of projection plane increments corresponding to each pixel side (axis).
        # The units of the returned results are the same as the units of cdelt, crval, and cd for the celestial WCS
        # and can be obtained by inquiring the value of cunit property of the input WCS WCS object.
        x_pixelscale = result[0] * u("deg")
        y_pixelscale = result[1] * u("deg")
        # pixelscale = Extent(x_pixelscale, y_pixelscale)

        semimajor = region.semimajor * x_pixelscale
        semiminor = region.semiminor * y_pixelscale

        radius = SkyStretch(semimajor, semiminor)

        # Create a new SkyEllipse
        return cls(center, radius, region.angle, meta=region.meta)

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Create a pixel ellipse region
        return PixelEllipseRegion.from_sky(self, wcs)

    # -----------------------------------------------------------------

    def to_physical(self, distance, wcs=None):

        """
        This function ...
        :param distance:
        :param wcs:
        :return:
        """

        return PhysicalEllipseRegion.from_sky(self, distance, wcs=wcs)

    # -----------------------------------------------------------------

    @property
    def unrotated_radius(self):

        """
        This function ...
        :return:
        """

        x_radius = self.radius.ra.to("arcsec").value
        y_radius = self.radius.dec.to("arcsec").value

        a_projected_x = x_radius * math.cos(self.angle.radian)
        b_projected_x = y_radius * math.sin(self.angle.radian)
        a_projected_y = x_radius * math.sin(self.angle.radian)
        b_projected_y = y_radius * math.cos(self.angle.radian)

        box_x_radius = max(abs(a_projected_x), abs(b_projected_x)) * u("arcsec")
        box_y_radius = max(abs(a_projected_y), abs(b_projected_y)) * u("arcsec")

        radius = SkyStretch(box_x_radius, box_y_radius)

        return radius

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
        #radunit = 'deg'
        radunit = "arcsec"

        if radunit == 'arcsec':
            if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
            else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else: radunitstr = ''

        # convert coordsys string to coordsys object
        if coordsys in coordsys_name_mapping: frame = frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
        else: frame = None  # for pixel/image/physical frames

        x = float(self.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(self.center.transform_to(frame).spherical.lat.to('deg').value)
        ra = coordinates.degrees_to_hms(ra=x)
        dec = coordinates.degrees_to_hms(dec=y)
        x = ra
        y = dec
        r1 = float(self.semimajor.to(radunit).value)
        r2 = float(self.semiminor.to(radunit).value)
        ang = float(self.angle.to('deg').value)

        #preprefix = coordsys + ";"
        preprefix = ""

        string = preprefix + prefix + make_ellipse_template(fmt, radunitstr, hmsdms=True).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class PhysicalEllipseRegion(EllipseRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, center, radius, angle=None, **kwargs):

        """
        The constructor ...
        :param center:
        :param radius:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be a physical coordinate")
        if not isinstance(radius, PhysicalStretch): raise ValueError("Radius must be a physical stretch")

        # Call the constructor of the base class
        super(PhysicalEllipseRegion, self).__init__(center, radius, angle, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_sky(cls, region, distance, wcs=None):

        """
        This function ...
        :param region:
        :param distance:
        :param wcs:
        :return:
        """

        # Convert radius
        radius = region.radius.to_physical(distance)

        # Convert center??
        if wcs is not None: center = region.center.to_physical()

        # No coordinate system
        else:
            warnings.warn("Center coordinate is lost during conversion from sky ellipse region to physical ellipse region when coordinate system is not passed")
            center = PhysicalCoordinate.zero()

        # Call the constructor
        return cls(center, radius, angle=region.angle)

# -----------------------------------------------------------------
