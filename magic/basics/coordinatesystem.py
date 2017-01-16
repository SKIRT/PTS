#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.coordinatesystem Contains the CoordinateSystem class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import copy
import numpy as np
import warnings

# Import astronomical modules
from astropy import wcs
from astropy.wcs import utils
from astropy.units import Unit
from astropy.io import fits
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .coordinate import PixelCoordinate, SkyCoordinate
from ..region.rectangle import SkyRectangleRegion
from .stretch import SkyStretch
from ..tools import coordinates
from .pixelscale import Pixelscale

# -----------------------------------------------------------------

class CoordinateSystem(wcs.WCS):

    """
    This class ...
    """

    def __init__(self, header=None):

        """
        This function ...
        :return:
        """

        # Check header
        coordinates_keywords = ["PC1_1", "PC1_2", "PC2_1", "PC2_2", "PC1", "PC2", "CDELT1", "CDELT2", "CD1_1", "CD1_2", "CD2_2", "CD2_1"]
        if header is not None:
            for keyword in coordinates_keywords:
                if keyword in header: break
            else: raise ValueError("Header does not contain coordinate information")

        # Call the constructor of the base class
        super(CoordinateSystem, self).__init__(header)

        # Set the number of pixels in both axis directions
        if header is not None:

            if not hasattr(self, "naxis1"):

                try: self.naxis1 = header["NAXIS1"]
                except TypeError: self.naxis1 = self._naxis1

            if not hasattr(self, "naxis2"):

                try: self.naxis2 = header["NAXIS2"]
                except TypeError: self.naxis2 = self._naxis2

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get the header and flatten it (remove references to third axis)
        header = fits.getheader(path)
        header["NAXIS"] = 2
        if "NAXIS3" in header: del header["NAXIS3"]
        for key in header:
            if "PLANE" in key: del header[key]

        self = cls(header)

        # Set the path
        self.path = path

        # Return
        return self

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        # DEFAULT COPY FUNCTION OF WCS IS NOT A REAL DEEP COPY!! THIS HAS CAUSED ME BUGS THAT GAVE ME HEADACHES

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return self.ysize, self.xsize

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.naxis1

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.naxis2

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        result = utils.proj_plane_pixel_scales(self)

        # returns: A vector (ndarray) of projection plane increments corresponding to each pixel side (axis).
        # The units of the returned results are the same as the units of cdelt, crval, and cd for the celestial WCS
        # and can be obtained by inquiring the value of cunit property of the input WCS WCS object.

        x_pixelscale = result[0] * Unit("deg/pix")
        y_pixelscale = result[1] * Unit("deg/pix")

        # Return the pixel scale as an extent
        return Pixelscale(x_pixelscale, y_pixelscale)

    # -----------------------------------------------------------------

    @property
    def average_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.pixelscale.average

    # -----------------------------------------------------------------

    @property
    def pixelarea(self):

        """
        This function ...
        :return:
        """

        square_pixel = 1.0 * Unit("pix2")
        return (self.pixelscale.x * self.pixelscale.y * square_pixel).to("sr")

    # -----------------------------------------------------------------

    @property
    def center_pixel(self):

        """
        This property ...
        :return:
        """

        x, y = self.wcs.crpix
        return PixelCoordinate(x, y)

    # -----------------------------------------------------------------

    @property
    def center_sky(self):

        """
        This property ...
        :return:
        """

        # Return the coordinate of the center in sky coordinates
        return SkyCoordinate.from_pixel(self.center_pixel, self)

    # -----------------------------------------------------------------

    @property
    def is_distorted(self):

        """
        This function ...
        :return:
        """

        return utils.is_proj_plane_distorted(self)

    # -----------------------------------------------------------------

    @property
    def is_flipped(self):

        """
        This function ...
        :return:
        """

        return (self.pixel_scale_matrix[0,0] >= 0) == (self.pixel_scale_matrix[1,1] >= 0)

    # -----------------------------------------------------------------

    def __eq__(self, other_wcs):

        """
        This function ...
        :param other_wcs:
        :return:
        """

        try:

            # Check whether the CRPIX is equal
            if not self.wcs.crpix[0] == other_wcs.wcs.crpix[0]: return False
            if not self.wcs.crpix[1] == other_wcs.wcs.crpix[1]: return False

            # Check whether the pixel scale is equal
            for element_self, element_other in zip(list(self.pixel_scale_matrix.flatten()), list(other_wcs.pixel_scale_matrix.flatten())):
                #print(element_self, element_other)
                if not element_self == element_other: return False

            # Check whether the CRVAL is equal
            if not self.wcs.crval[0] == other_wcs.wcs.crval[0]: return False
            if not self.wcs.crval[1] == other_wcs.wcs.crval[1]: return False

            # Check whether the number of axis is equal
            if not self.naxis == other_wcs.naxis: return False

            # Check whether the axis sizes are equal
            if not self.naxis1 == other_wcs.naxis1: return False
            if not self.naxis2 == other_wcs.naxis2: return False

            # If all of the above tests succeeded, the coordinate systems may be considered equal
            return True

        except AttributeError: return False

    # -----------------------------------------------------------------

    def to_pixel(self, shape):

        """
        This function ...
        :param shape:
        :return:
        """

        return shape.to_pixel(self)

    # -----------------------------------------------------------------

    def to_sky(self, shape):

        """
        This function ...
        :param shape:
        :return:
        """

        return shape.to_sky(self)

    # -----------------------------------------------------------------

    @property
    def coordinate_range(self):

        """
        This function ...
        :return:
        """

        coor1 = self.wcs_pix2world(0.0, 0.0, 0)
        coor2 = self.wcs_pix2world(self.xsize - 1.0, self.ysize - 1.0, 0)

        co1 = SkyCoordinate(ra=float(coor1[0]), dec=float(coor1[1]), unit="deg", frame='fk5')
        co2 = SkyCoordinate(ra=float(coor2[0]), dec=float(coor2[1]), unit="deg", frame='fk5')

        #print("co1=", co1.to_string('hmsdms'))
        #print("co2=", co2.to_string('hmsdms'))

        ra_range = [co1.ra.value, co2.ra.value]
        dec_range = [co1.dec.value, co2.dec.value]

        # Determine the center in RA and DEC (in degrees)
        ra_center = 0.5 * (ra_range[0] + ra_range[1])
        dec_center = 0.5 * (dec_range[0] + dec_range[1])

        # New
        dec_begin = dec_range[0]
        dec_end = dec_range[1]
        ra_begin = ra_range[0]
        ra_end = ra_range[1]

        # Calculate the actual RA and DEC distance in degrees
        ra_distance = abs(coordinates.ra_distance(dec_center, ra_begin, ra_end))
        dec_distance = abs(dec_end - dec_begin)

        # Calculate the pixel scale of this image in degrees
        x_pixelscale_deg = self.pixelscale.x.to("deg").value
        y_pixelscale_deg = self.pixelscale.y.to("deg").value

        # Get the center pixel
        ref_pix = self.wcs.crpix
        ref_world = self.wcs.crval

        center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame="fk5")

        # Get the orientation of the coordinate system
        try: orientation = self.orientation
        except ValueError:
            largest_distance = max(ra_distance, dec_distance) * Unit("deg")
            return center, largest_distance, largest_distance

        if "x" in orientation[0] and "y" in orientation[1]: # RA axis = x axis and DEC axis = y axis

            size_ra_deg = self.xsize * x_pixelscale_deg
            size_dec_deg = self.ysize * y_pixelscale_deg

        elif "y" in orientation[0] and "x" in orientation[1]: # RA axis = y axis and DEC axis = x axis

            size_ra_deg = self.ysize * y_pixelscale_deg
            size_dec_deg = self.xsize * x_pixelscale_deg

        else: raise ValueError("Invalid coordinate system orientation:" + str(orientation))

        # Check whether the two different ways of calculating the RA width result in approximately the same value
        #assert np.isclose(ra_distance, size_ra_deg, rtol=0.05), "The coordinate system and pixel scale do not match: ra_distance=" + str(ra_distance) + ",size_ra_deg=" + str(size_ra_deg)
        #assert np.isclose(dec_distance, size_dec_deg, rtol=0.05), "The coordinate system and pixel scale do not match: dec_distance=" + str(dec_distance) + ",size_dec_deg=" + str(size_dec_deg)

        if not np.isclose(ra_distance, size_ra_deg, rtol=0.05):
            warnings.warn("The coordinate system and pixel scale do not match: ra_distance = " + str(ra_distance) + ", size_ra_deg = " + str(size_ra_deg))
        if not np.isclose(dec_distance, size_dec_deg, rtol=0.05):
            warnings.warn("The coordinate system and pixel scale do not match: dec_distance = " + str(dec_distance) + ", size_dec_deg = " + str(size_dec_deg))

        # Create RA and DEC span as quantities
        ra_span = ra_distance * Unit("deg")
        dec_span = dec_distance * Unit("deg")

        # Return the center coordinate and the RA and DEC span
        return center, ra_span, dec_span

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Get coordinate range
        center, ra_span, dec_span = self.coordinate_range

        #ra = center.ra.to(unit).value
        #dec = center.dec.to(unit).value

        #ra_span = ra_span.to(unit).value
        #dec_span = dec_span.to(unit).value

        # Create rectangle
        #center = Position(ra, dec)
        #radius = Extent(0.5 * ra_span, 0.5 * dec_span)
        #box = Rectangle(center, radius)

        radius = SkyStretch(0.5 * ra_span, 0.5 * dec_span)

        box = SkyRectangleRegion(center, radius)

        # Return the box
        return box

    # -----------------------------------------------------------------

    @property
    def orientation(self):

        """
        This property ...
        :return:
        """

        # Calculate the number of pixels in the x and y direction that corresponds one arcmin
        pix_x = (1. / self.pixelscale.x * Unit("arcmin")).to("pix").value
        pix_y = (1. / self.pixelscale.y * Unit("arcmin")).to("pix").value

        #print(pix_x, pix_y)

        # Get the center pixel of the coordinate system
        center = self.center_pixel
        center_coordinate = self.center_sky

        # Calculate the new pixel position
        new_x_pixel = PixelCoordinate(center.x + pix_x, center.y)
        new_y_pixel = PixelCoordinate(center.x, center.y + pix_y)

        # Convert to sky coordinate
        new_x_sky = self.to_sky(new_x_pixel)
        new_y_sky = self.to_sky(new_y_pixel)

        #print("center coordinate:", center_coordinate)
        #print("new x coordinate:", new_x_sky)
        #print("new y coordinate:", new_y_sky)
        #print()
        #print("difference with increment in x:")
        difference_x_ra = coordinates.ra_distance(center_coordinate.dec.degree, center_coordinate.ra.degree, new_x_sky.ra.degree)
        difference_x_ra = np.sign(new_x_sky.ra.degree-center_coordinate.ra.degree) * difference_x_ra # CORRECT FOR SIGN !! RA_DISTANCE IS ABSOLUTE VALUE (BECAUSE OF ARCCOS!!!)
        difference_x_dec = new_x_sky.dec.degree - center_coordinate.dec.degree

        difference_x_ra_arcmin = difference_x_ra * 60.
        difference_x_dec_arcmin = difference_x_dec * 60.

        #print("  RA:", difference_x_ra, "deg =", difference_x_ra*60., "arcmin")
        #print("  DEC:", difference_x_dec, "deg =", difference_x_dec*60., "arcmin")
        #print()
        #print("difference with increment in y:")
        difference_y_ra = coordinates.ra_distance(center_coordinate.dec.degree, center_coordinate.ra.degree, new_y_sky.ra.degree)
        difference_y_ra = np.sign(new_y_sky.ra.degree-center_coordinate.ra.degree) * difference_y_ra # CORRECT FOR SIGN !! RA_DISTANCE IS ABSOLUTE VALUE (BECAUSE OF ARCCOS!!!)
        difference_y_dec = new_y_sky.dec.degree - center_coordinate.dec.degree

        difference_y_ra_arcmin = difference_y_ra * 60.
        difference_y_dec_arcmin = difference_y_dec * 60.

        #print("  RA:", difference_y_ra, "deg =", difference_y_ra*60., "arcmin")
        #print("  DEC:", difference_y_dec, "deg=", difference_y_dec*60., "arcmin")

        ra_alignment = None
        dec_alignment = None

        # Get the orientation of the x and y axes with respect to RA and DEC
        ra_alignment_x, dec_alignment_x = get_alignment(difference_x_ra_arcmin, difference_x_dec_arcmin)
        ra_alignment_y, dec_alignment_y = get_alignment(difference_y_ra_arcmin, difference_y_dec_arcmin)

        if ra_alignment_x is not None:

            assert dec_alignment_x is None
            assert ra_alignment_y is None

            alignment = (ra_alignment_x + "x", dec_alignment_y + "y")

        elif ra_alignment_y is not None:

            assert dec_alignment_y is None
            assert ra_alignment_x is None

            alignment = (ra_alignment_y + "y", dec_alignment_x + "x")

        else: raise ValueError("Orientation of the coordinate system is unclear")

        return alignment

    # -----------------------------------------------------------------

    @property
    def orientation_angle(self):

        """
        This function ...
        :return:
        """

        diag_a = self.pixel_scale_matrix[0,1]
        diag_b = self.pixel_scale_matrix[1,0]

        if not np.isclose(diag_a, diag_b, rtol=0.05):
            warnings.warn("The diagonal elements of the pixel scale matrix are not equal: " + repr(diag_a) + " and " + repr(diag_b))

        first = self.pixel_scale_matrix[0,0]

        radians = np.arctan(diag_a / first)

        degrees = radians / math.pi * 180.

        return Angle(degrees, "deg")

    # -----------------------------------------------------------------

    @property
    def standard_orientation_angle(self):

        """
        This function ...
        :return:
        """

        # Angles of 0, +90, +180, -90 ...

        orientation = self.orientation

        # Check if not flipped
        if "x" in orientation[0]:

            sign_x = orientation[0][0]
            sign_y = orientation[1][0]

            if sign_x == sign_y: raise ValueError("The coordinate system is flipped!")

        else:

            sign_x = orientation[1][0]
            sign_y = orientation[0][0]

            if sign_x != sign_y: raise ValueError("The coordinate system is flipped!")

        # There are only four allowed options:
        # (-y, -x), (+x, -y), (+y, +x) and (-x,+y)

        if orientation == ("-x", "+y"): return Angle(0.0, "deg") # This is the common orientation (standard East-West projection and South-North projection)
        elif orientation == ("-y", "-x"): return Angle(-90., "deg")
        elif orientation == ("+x", "-y"): return Angle(180., "deg")
        elif orientation == ("+y", "+x"): return Angle(90, "deg")
        else: raise ValueError("Not a standard orientation angle")

    # -----------------------------------------------------------------

    @property
    def orientation_angle_flipped(self):

        """
        This function ...
        :return:
        """

        # TODO: actually solve the problem of flipped coordinate systems

        orientation = self.orientation

        # There are only four allowed options:
        # (-y, -x), (+x, -y), (+y, +x) and (-x,+y)

        # here I just switched y and x as compared to above
        if orientation == ("-y", "+x"): return Angle(0.0, "deg") # This is the common orientation (standard East-West projection and South-North projection)
        elif orientation == ("-x", "-y"): return Angle(-90., "deg")
        elif orientation == ("+y", "-x"): return Angle(180., "deg")
        elif orientation == ("+x", "+y"): return Angle(90, "deg")
        else: raise ValueError("Not a standard orientation angle")

    # -----------------------------------------------------------------

    def to_header(self, relax=None, key=None):

        """
        This function ...
        :param relax:
        :param key:
        :return:
        """

        header = super(CoordinateSystem, self).to_header(relax, key)

        #header["NAXIS1"] = self._naxis1
        #header["NAXIS2"] = self._naxis2

        # Add the cards to the beginning
        header.insert(0, ("NAXIS2", self._naxis2))
        header.insert(0, ("NAXIS1", self._naxis1))

        return header

    # -----------------------------------------------------------------

    def to_astropy(self):

        """
        This function ...
        :return:
        """

        return wcs.WCS(self.to_header())

    # -----------------------------------------------------------------

    @property
    def csys(self):

        """
        This function ...
        :return:
        """

        ctype = self.wcs.ctype[0]
        if 'RA' in ctype or 'DEC' in ctype:
            if self.wcs.equinox == 2000 or self.wcs.equinox == 2000.:
                return 'fk5'
            elif self.wcs.equinox == 1950 or self.wcs.equinox == 1950.:
                return 'fk4'
            else:
                raise NotImplementedError("Non-fk4/fk5 equinoxes are not allowed")
        elif 'GLON' in ctype or 'GLAT' in ctype:
            return 'galactic'

# -----------------------------------------------------------------

def get_alignment(ra_difference, dec_difference):

    """
    This function ...
    :param ra_difference:
    :param dec_difference:
    :return:
    """

    ra_alignment = None
    dec_alignment = None

    if np.isclose(ra_difference, 1.0, atol=0.0001):  # if RA increment is close to 1 with a positive x increment

        # Test wether the DEC increment is then close to zero
        if np.isclose(dec_difference, 0.0, atol=0.01):

            ra_alignment = "+"

        # Else, something else is going on
        else: raise ValueError("The orientation of the coordinate system is unclear: " + str(ra_difference) + " " + str(dec_difference))

    elif np.isclose(ra_difference, -1.0, atol=0.0001):

        # Test whether the DEC increment is then close to zero
        if np.isclose(dec_difference, 0.0, atol=0.01):

            ra_alignment = "-"

        # Else, something else is going on
        else: raise ValueError("The orientation of the coordinate system is unclear: " + str(ra_difference) + " " + str(dec_difference))

    elif np.isclose(ra_difference, 0.0, atol=0.01):

        # Test whether the DEC increment is then close to one
        if np.isclose(dec_difference, 1.0, atol=0.0001):

            dec_alignment = "+"

        # Or close to -1 ..
        elif np.isclose(dec_difference, -1.0, atol=0.0001):

            dec_alignment = "-"

        # If that is zero, something else is going on
        else: raise ValueError("The orientation of the coordinate system is unclear: " + str(ra_difference) + " " + str(dec_difference))

    else: raise ValueError("The orientation of the coordinate system is unclear: " + str(ra_difference) + " " + str(dec_difference))

    return ra_alignment, dec_alignment

# -----------------------------------------------------------------
