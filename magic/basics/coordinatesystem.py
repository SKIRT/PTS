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
import traceback
import copy
import numpy as np
import warnings

# Import astronomical modules
from astropy import wcs
from astropy.wcs import utils
from astropy.io import fits
from astropy.coordinates import Angle
from astropy.io.fits import Header

# Import the relevant PTS classes and modules
from .coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..region.rectangle import SkyRectangleRegion
from .stretch import SkyStretch
from ..tools import coordinates
from .pixelscale import Pixelscale, PhysicalPixelscale
from ...core.units.parsing import parse_unit as u
from ...core.basics.range import QuantityRange
from .stretch import PixelStretch
from ...core.basics.log import log
from ...core.tools import sequences
from ...core.tools import types
from .vector import Position

# -----------------------------------------------------------------

class CoordinateSystem(wcs.WCS):

    """
    This class ...
    """

    def __init__(self, header=None, naxis=None, naxis1=None, naxis2=None):

        """
        This function ...
        :param header:
        :param naxis:
        :param naxis1:
        :param naxis2:
        :return:
        """

        # Check header
        coordinates_keywords = ["PC1_1", "PC1_2", "PC2_1", "PC2_2", "PC1", "PC2", "CDELT1", "CDELT2", "CD1_1", "CD1_2", "CD2_2", "CD2_1"]
        if header is not None:
            for keyword in coordinates_keywords:
                if keyword in header: break
            else: raise ValueError("Header does not contain coordinate information")

        # Call the constructor of the base class
        super(CoordinateSystem, self).__init__(header=header, naxis=naxis)

        # Set the number of pixels in both axis directions
        if header is not None:

            if not hasattr(self, "naxis1"):

                try: self.naxis1 = header["NAXIS1"]
                except KeyError: self.naxis1 = self._naxis1
                #except TypeError: self.naxis1 = self._naxis1

            if not hasattr(self, "naxis2"):

                try: self.naxis2 = header["NAXIS2"]
                except KeyError: self.naxis2 = self._naxis2
                #except TypeError: self.naxis2 = self._naxis2

        # If naxis1 is specified
        if naxis1 is not None: self.naxis1 = naxis1

        # If naxis2 is specified
        if naxis2 is not None: self.naxis2 = naxis2

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_header(cls, header):

        """
        This function ...
        :param header:
        :return:
        """

        header["NAXIS"] = 2
        if "NAXIS3" in header: del header["NAXIS3"]
        for key in header:
            if "PLANE" in key: del header[key]

        # Create
        return cls(header=header)

    # -----------------------------------------------------------------

    @classmethod
    def from_header_string(cls, string):

        """
        This function ...
        :param string:
        :return:
        """

        header = Header.fromstring(string)
        return cls.from_header(header)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        if path.endswith(".fits"): return cls.from_fits(path)
        elif path.endswith(".txt"): return cls.from_header_file(path)
        else: raise ValueError("Path must be a FITS or header text file")

    # -----------------------------------------------------------------

    @classmethod
    def from_fits(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get the header and flatten it (remove references to third axis)
        header = fits.getheader(path)
        wcs = cls.from_header(header)
        wcs.path = path
        return wcs

    # -----------------------------------------------------------------

    @classmethod
    def from_header_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load the header
        header = Header.fromtextfile(path)

        # Create and return the coordinate system
        system = cls(header=header)
        system.path = path
        return system

    # -----------------------------------------------------------------

    @classmethod
    def from_properties(cls, size, reference_pixel, center_sky, pixelscale):

        """
        This fucntion ...
        :param size: 
        :param reference_pixel:
        :param center_sky: 
        :param pixelscale: 
        :return: 
        """

        # TODO: doesn't work
        # see try below

        # Construct header
        header = Header()

        header["CRPIX1"] = reference_pixel.x
        header["CRVAL1"] = center_sky.ra.to("deg").value
        header["CDELT1"] = pixelscale.x.to("deg").value
        header["CRPIX2"] = reference_pixel.y
        header["CRVAL2"] = center_sky.dec.to("deg").value
        header["CDELT2"] = pixelscale.y.to("deg").value
        header["NAXIS1"] = size.x
        header["NAXIS2"] = size.y
        header["RADESYS"] = "ICRS"
        header["CTYPE1"] = "RA--TAN"
        header["CTYPE2"] = "DEC--TAN"

        # Create and return
        return cls(header=header)

    # -----------------------------------------------------------------

    @classmethod
    def from_properties_try(cls, size, reference_pixel, center_sky, pixelscale):

        """
        This function ...
        :param size:
        :param reference_pixel:
        :param center_sky:
        :param pixelscale:
        :return:
        """

        # Create a new WCS object.  The number of axes must be set
        # from the start
        ##system = wcs.WCS(naxis=2)
        #system = cls(naxis=2, naxis1=size.x, naxis2=size.y)

        # Set up an "Airy's zenithal" projection
        # Vector properties may be set with Python lists, or Numpy arrays
        #system.wcs.crpix = [center_pixel.x, center_pixel.y]
        #system.wcs.cdelt = np.array([pixelscale.x.to("deg").value, pixelscale.y.to("deg").value])
        ##crval = [center_sky.ra.to("deg").value, center_sky.dec.to("deg").value]
        ##print(crval)
        #system.wcs.crval = [center_sky.ra.to("deg").value, center_sky.dec.to("deg").value]
        ##system.wcs.ctype = ["RA---AIR", "DEC--AIR"]
        #system.wcs.ctype = ["RA--TAN", "DEC--TAN"]
        #system.wcs.set_pv([(2, 1, 45.0)])

        center_x = reference_pixel.x
        center_y = reference_pixel.y

        # Create a new WCS object.  The number of axes must be set
        # from the start
        w = wcs.WCS(naxis=2)

        # Set up an "Airy's zenithal" projection
        # Vector properties may be set with Python lists, or Numpy arrays
        #w.wcs.crpix = [-center_x, center_y]
        #w.wcs.cdelt = np.array([-pixelscale.x.to("deg").value, pixelscale.y.to("deg").value])
        #w.wcs.crval = [center_sky.ra.to("deg").value, center_sky.dec.to("deg").value]
        #w.wcs.ctype = ["RA--TAN", "DEC--TAN"]
        #w.wcs.set_pv([(2, 1, 45.0)])

        # TODO: make this work with not only this example !!!

        # Set up an "Airy's zenithal" projection
        # Vector properties may be set with Python lists, or Numpy arrays
        w.wcs.crpix = [-234.75, 8.3393]
        w.wcs.cdelt = np.array([-0.066667, 0.066667])
        w.wcs.crval = [0, -90]
        w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
        w.wcs.set_pv([(2, 1, 45.0)])

        # Return the coordinate system
        #return system

        header = w.to_header()

        # Construct header
        #header = Header()
        #header["SIMPLE"] = True
        #header["BITPIX"] = -32
        #header["NAXIS"] = 2
        #header["NAXIS1"] = size.x
        #header["NAXIS2"] = size.y
        #header["CRVAL1"] = center_sky.ra.to("deg").value
        #header["CRVAL2"] = center_sky.dec.to("deg").value
        #header["RADESYS"] = "ICRS"
        #header["CTYPE1"] = "RA--TAN"
        #header["CTYPE2"] = "DEC--TAN"
        #header["CRPIX1"] = center_x
        #header["CRPIX2"] = center_y
        #header["CDELT1"] = pixelscale.x.to("deg").value
        #header["CDELT2"] = pixelscale.y.to("deg").value

        #print(header)

        #from ...core.tools import filesystem as fs
        #from ...core.tools import introspection
        #path = fs.join(introspection.pts_temp_dir, "wcs_from_properties.fits")
        #if fs.is_file(path): fs.remove_file(path)

        #data = np.zeros((size.y, size.x))
        #hdu = fits.PrimaryHDU(data, header)
        #hdulist = fits.HDUList([hdu])
        #hdulist.writeto(path)

        # Convert into wcs
        #return cls(header=header, naxis=2)

        #return cls.from_fits(path)

        return cls(header, naxis=2, naxis1=size.x, naxis2=size.y)

    # -----------------------------------------------------------------

    @classmethod
    def from_ranges(cls, ra_range, dec_range, pixelscale):

        """
        This function ...
        :param ra_range:
        :param dec_range:
        :param pixelscale:
        :return:
        """

        # Determine the total RA and DEC distance
        ra_begin = ra_range.min
        ra_end = ra_range.max
        dec_begin = dec_range.min
        dec_end = dec_range.max
        dec_center = dec_range.center
        ra_distance = abs(coordinates.ra_distance(dec_center, ra_begin, ra_end))
        dec_distance = abs(dec_end - dec_begin)

        # Detemrine the nubmerof required pixels in the x and y directions
        xpixels = int(ra_distance / pixelscale.x)
        ypixels = int(dec_distance / pixelscale.y)

        # Set size and center pixel
        size = PixelStretch(xpixels, ypixels)
        center_pixel = PixelCoordinate(0.5 * (size.x + 1), 0.5 * (size.y + 1))

        # Determine center coordinate
        ra_center = ra_range.center
        dec_center = dec_range.center
        center = SkyCoordinate(ra=ra_center, dec=dec_center)

        # Create and return
        return cls.from_properties(size, center_pixel, center, pixelscale)

    # -----------------------------------------------------------------

    @classmethod
    def from_projection(cls, projection, center_coordinate):

        """
        This function ...
        :param projection:
        :param center_coordinate:
        :return:
        """

        # Get properties
        pixels_x = projection.pixels_x
        pixels_y = projection.pixels_y
        center_x = projection.center_x
        center_y = projection.center_y

        # Construct header
        header = Header()

        header["CRPIX1"] = center_x
        header["CRVAL1"] = center_coordinate.ra.to("deg").value
        header["CDELT1"] = projection.pixelscale.x
        header["CRPIX2"] = center_y
        header["CRVAL2"] = center_coordinate.dec.to("deg").value
        header["CDELT2"] = projection.pixelscale.y
        header["NAXIS1"] = pixels_x
        header["NAXIS2"] = pixels_y

        # Create and return
        return cls(header=header)

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
    def ctype(self):

        """
        This function ...
        :return:
        """

        return self.wcs.ctype

    # -----------------------------------------------------------------

    @property
    def ctype1(self):

        """
        Thisfunction ...
        :return:
        """

        return self.ctype[0]

    # -----------------------------------------------------------------

    @property
    def ctype2(self):

        """
        This function ...
        :return:
        """

        return self.ctype[1]

    # -----------------------------------------------------------------

    @property
    def crval(self):

        """
        This function ...
        :return:
        """

        return self.wcs.crval

    # -----------------------------------------------------------------

    @property
    def crval1(self):

        """
        This function ...
        :return:
        """

        return self.crval[0]

    # -----------------------------------------------------------------

    @property
    def crval2(self):

        """
        This function ...
        :return:
        """

        return self.crval[1]

    # -----------------------------------------------------------------

    @property
    def reference_pixel(self):

        """
        This function ...
        :return:
        """

        x, y = self.wcs.crpix
        return PixelCoordinate(x, y)

    # -----------------------------------------------------------------

    @property
    def reference_coordinate(self):

        """
        This function ...
        :return:
        """

        # Celestial
        if self.is_celestial: return SkyCoordinate(ra=self.crval1 * self.axis1_unit, dec=self.crval2 * self.axis2_unit)

        # Physical
        elif self.is_physical: return PhysicalCoordinate(self.crval1 * self.axis1_unit, self.crval2 * self.axis2_unit)

        # Unknown
        else: raise ValueError("Unknown coordinate system type")

    # -----------------------------------------------------------------

    @property
    def axis1_unit(self):

        """
        Thisn function ...
        :return:
        """

        # Try to parse the ctype1 as a unit
        try: unit1 = u(self.ctype1)
        except ValueError: unit1 = u("deg")

        # Return the unit
        return unit1

    # -----------------------------------------------------------------

    @property
    def axis2_unit(self):

        """
        This function ...
        :return:
        """

        # Try to parse the ctype2 as a unit
        try: unit2 = u(self.ctype2)
        except ValueError: unit2 = u("deg")

        # Return the unit
        return unit2

    # -----------------------------------------------------------------

    @property
    def axis_units(self):
        return [self.axis1_unit, self.axis2_unit]

    # -----------------------------------------------------------------

    @property
    def is_sky(self):
        return self.is_celestial

    # -----------------------------------------------------------------

    @property
    def is_physical(self):

        """
        This function ...
        :return:
        """

        # Physical
        if types.is_length_unit(self.axis1_unit):

            # Check
            if not types.is_length_unit(self.axis2_unit): raise ValueError("Coordinate system is a mix of physical and celestial axis")

            # True
            return True

        # Not physical
        else:

            # Check
            if types.is_length_unit(self.axis2_unit): raise ValueError("Coordinate system is a mix of physical and celestial axis")

            # False
            return False

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

        x_pixelscale = result[0] * self.axis1_unit
        y_pixelscale = result[1] * self.axis2_unit

        #print(self.axis1_unit, self.axis2_unit, self.is_celestial, self.is_physical)
        #print(self.naxis)

        # Return the pixel scale as an extent
        if self.is_celestial: return Pixelscale(x_pixelscale, y_pixelscale)
        elif self.is_physical: return PhysicalPixelscale(x_pixelscale, y_pixelscale)
        else: raise ValueError("Unknown coordinate system type")

    # -----------------------------------------------------------------

    @property
    def average_pixelscale(self):
        return self.pixelscale.average

    # -----------------------------------------------------------------

    @property
    def pixelarea(self):

        """
        This function ...
        :return:
        """

        if self.is_celestial: return (self.pixelscale.x * self.pixelscale.y).to("sr")
        elif self.is_physical: return self.pixelscale.x * self.pixelscale.y
        else: raise ValueError("Unknown coordinate system type")

    # -----------------------------------------------------------------

    @property
    def npixels(self):
        return self.naxis1 * self.naxis2

    # -----------------------------------------------------------------

    @property
    def area(self):
        return self.pixelarea * self.npixels

    # -----------------------------------------------------------------

    @property
    def center_sky(self):
        # Return the coordinate of the center in sky coordinates
        return SkyCoordinate.from_pixel(self.reference_pixel, self)

    # -----------------------------------------------------------------

    @property
    def is_distorted(self):
        return utils.is_proj_plane_distorted(self)

    # -----------------------------------------------------------------

    @property
    def is_flipped(self):
        return (self.pixel_scale_matrix[0,0] >= 0) == (self.pixel_scale_matrix[1,1] >= 0)

    # -----------------------------------------------------------------

    def scale(self, new_xsize, new_ysize):

        """
        This function ...
        :param new_xsize:
        :param new_ysize:
        :return:
        """

        # Determine the new position of the reference pixel
        relative_center = Position(self.reference_pixel.x / self.xsize, self.reference_pixel.y / self.ysize)
        new_center = Position(relative_center.x * new_xsize, relative_center.y * new_ysize)

        # Change the center pixel position
        self.wcs.crpix[0] = new_center.x
        self.wcs.crpix[1] = new_center.y

        # Get original xsize and ysize
        original_xsize = self.xsize
        original_ysize = self.ysize

        # Change the number of pixels of the axes
        self.naxis1 = new_xsize
        self.naxis2 = new_ysize
        self._naxis1 = self.naxis1
        self._naxis2 = self.naxis2

        #print(self.pixel_scale_matrix)
        #print(self.wcs.cdelt)

        # Change the pixel scale
        self.wcs.cdelt[0] *= float(original_xsize) / float(new_xsize)
        self.wcs.cdelt[1] *= float(original_ysize) / float(new_ysize)

        # Return the new center and new pixelscale
        return self.reference_pixel, self.pixelscale

    # -----------------------------------------------------------------

    def scaled(self, new_xsize, new_ysize):

        """
        This function ...
        :param new_xsize:
        :param new_ysize:
        :return:
        """

        # Make a copy
        new_wcs = self.copy()

        # Scale
        new_wcs.scale(new_xsize, new_ysize)

        # Return
        return new_wcs

    # -----------------------------------------------------------------

    def __eq__(self, other_wcs):

        """
        This function ...
        :param other_wcs:
        :return:
        """

        try:
            comparisons = self.make_comparisons(other_wcs)
            return sequences.all_true(comparisons)

        # Error occured
        except AttributeError as e:

            log.warning("Problem occurred while comparing coordinate systems:")
            print(e)
            traceback.print_exc()
            return False

    # -----------------------------------------------------------------

    def make_comparisons(self, other_wcs):

        """
        This function ...
        :param other_wcs:
        :return:
        """

        result = []

        # Check whether the CRPIX is equal
        if not self.wcs.crpix[0] == other_wcs.wcs.crpix[0]:
            #print("1")
            #return False
            result.append(False)
        else: result.append(True)

        if not self.wcs.crpix[1] == other_wcs.wcs.crpix[1]:
            #print("2")
            #return False
            result.append(False)
        else: result.append(True)

        # Check whether the pixel scale is equal
        for element_self, element_other in zip(list(self.pixel_scale_matrix.flatten()),
                                               list(other_wcs.pixel_scale_matrix.flatten())):
            # print(element_self, element_other)
            if element_self != element_other:
                #print("3")
                #return False
                result.append(False)
            else: result.append(True)

        # Check whether the CRVAL is equal
        if not self.wcs.crval[0] == other_wcs.wcs.crval[0]:
            #print("4")
            #return False
            result.append(False)
        else: result.append(True)

        if not self.wcs.crval[1] == other_wcs.wcs.crval[1]:
            #print("5")
            #return False
            result.append(False)
        else: result.append(True)

        # Check whether the number of axis is equal
        if not self.naxis == other_wcs.naxis:
            #print("6")
            #return False
            result.append(False)
        else: result.append(True)

        # Check whether the axis sizes are equal
        if not self.naxis1 == other_wcs.naxis1:
            #print("7")
            #return False
            result.append(False)
        else: result.append(True)

        if not self.naxis2 == other_wcs.naxis2:
            #print("8")
            #return False
            result.append(False)
        else: result.append(True)

        # Return the comparison results
        return result

    # -----------------------------------------------------------------

    def __ne__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return not self.__eq__(other)

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
    def ra_range(self):

        """
        This function ...
        :return:
        """

        coor1 = self.wcs_pix2world(0.0, 0.0, 0)
        coor2 = self.wcs_pix2world(self.xsize - 1.0, self.ysize - 1.0, 0)

        co1 = SkyCoordinate(ra=float(coor1[0]), dec=float(coor1[1]), unit="deg", frame='fk5')
        co2 = SkyCoordinate(ra=float(coor2[0]), dec=float(coor2[1]), unit="deg", frame='fk5')

        # Get the range
        ra_range = sorted([co1.ra.value, co2.ra.value])

        # Return the range
        return QuantityRange(ra_range[0], ra_range[1], "deg")

    # -----------------------------------------------------------------

    @property
    def dec_range(self):

        """
        This function ...
        :return:
        """

        coor1 = self.wcs_pix2world(0.0, 0.0, 0)
        coor2 = self.wcs_pix2world(self.xsize - 1.0, self.ysize - 1.0, 0)

        co1 = SkyCoordinate(ra=float(coor1[0]), dec=float(coor1[1]), unit="deg", frame='fk5')
        co2 = SkyCoordinate(ra=float(coor2[0]), dec=float(coor2[1]), unit="deg", frame='fk5')

        # Get the range
        dec_range = sorted([co1.dec.value, co2.dec.value])

        # Return the range
        return QuantityRange(dec_range[0], dec_range[1], "deg")

    # -----------------------------------------------------------------

    @property
    def min_ra(self):

        """
        This function ...
        :return:
        """

        return self.ra_range.min

    # -----------------------------------------------------------------

    @property
    def max_ra(self):

        """
        This function ...
        :return:
        """

        return self.ra_range.max

    # -----------------------------------------------------------------

    @property
    def min_dec(self):

        """
        This function ...
        :return:
        """

        return self.dec_range.min

    # -----------------------------------------------------------------

    @property
    def max_dec(self):

        """
        THis function ...
        :return:
        """

        return self.dec_range.max

    # -----------------------------------------------------------------

    def contains(self, coordinate):

        """
        This function ...
        :param coordinate: 
        :return: 
        """

        if coordinate.ra < self.min_ra: return False
        if coordinate.ra > self.max_ra: return False
        if coordinate.dec < self.min_dec: return False
        if coordinate.dec > self.max_dec: return False
        return True

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

        # NEW METHOD FROM ASTROPY
        ra_distance_top = abs(SkyCoordinate(ra=ra_begin, dec=dec_end, unit="deg").separation(SkyCoordinate(ra=ra_end, dec=dec_end, unit="deg")).deg)
        ra_distance_bottom = abs(SkyCoordinate(ra=ra_begin, dec=dec_begin, unit="deg").separation(SkyCoordinate(ra=ra_end, dec=dec_begin, unit="deg")).deg)
        dec_distance_new = abs(SkyCoordinate(ra=ra_begin, dec=dec_begin, unit="deg").separation(SkyCoordinate(ra=ra_begin, dec=dec_end, unit="deg")).deg)

        # Checks
        #assert np.isclose(ra_distance_top, ra_distance_bottom, rtol=0.05), (ra_distance_top, ra_distance_bottom)
        #assert np.isclose(ra_distance_top, ra_distance, rtol=0.05), (ra_distance_top, ra_distance)
        assert np.isclose(dec_distance_new, dec_distance, rtol=0.05), (dec_distance_new, dec_distance)

        #if not np.isclose(ra_distance_top, ra_distance_bottom, rtol=0.05):
        #    log.warning("RA distance at top of image is " + str(ra_distance_top) + " whereas at the bottom is " + str(ra_distance_bottom))
        #    log.warning("RA distance at center is " + str(ra_distance))

        #if not np.isclose(ra_distance_top, ra_distance, rtol=0.05):
        #    log.warning("RA distance at top of image is " + str(ra_distance_top) + " whereas at the center it is " + str(ra_distance))

        # Set the RA distance to the maximum of ra_distance, ra_distance_bottom, and ra_distance_top
        ra_distance = max(ra_distance_bottom, ra_distance_top)

        # Calculate the pixel scale of this image in degrees
        x_pixelscale_deg = self.pixelscale.x.to("deg").value
        y_pixelscale_deg = self.pixelscale.y.to("deg").value

        # Get the center pixel
        #ref_pix = self.wcs.crpix
        #ref_world = self.wcs.crval

        center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame="fk5")

        # Get the orientation of the coordinate system
        try: orientation = self.orientation
        except ValueError:
            largest_distance = max(ra_distance, dec_distance) * u("deg")
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

        # Checks
        if not np.isclose(ra_distance, size_ra_deg, rtol=0.05):
            warnings.warn("The coordinate system and pixel scale do not match: ra_distance = " + str(ra_distance) + ", size_ra_deg = " + str(size_ra_deg))
        if not np.isclose(dec_distance, size_dec_deg, rtol=0.05):
            warnings.warn("The coordinate system and pixel scale do not match: dec_distance = " + str(dec_distance) + ", size_dec_deg = " + str(size_dec_deg))

        # Create RA and DEC span as quantities
        ra_span = ra_distance * u("deg")
        dec_span = dec_distance * u("deg")

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
        pix_x = (1. / self.pixelscale.x * u("arcmin")).to("").value
        pix_y = (1. / self.pixelscale.y * u("arcmin")).to("").value

        #print(pix_x, pix_y)

        # Get the center pixel of the coordinate system
        center = self.reference_pixel
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

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the coordinate system ...")

        # Check whether the path is defined
        if self.path is None: raise RuntimeError("Path is not defined for this coordinate system")

        # Save the image
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        if path.endswith(".fits"): self.saveto_fits(path)
        elif path.endswith(".txt"): self.saveto_header_file(path)
        else: raise ValueError("Path should be a FITS or txt file")

    # -----------------------------------------------------------------

    def saveto_fits(self, path):

        """
        This function ...
        :return:
        """

        # Convert to header
        header = self.to_header()

        # Write as FITS file
        header.tofile(path)

    # -----------------------------------------------------------------

    def saveto_header_file(self, path):

        """
        This function ...
        :return:
        """

        # Convert to header
        header = self.to_header()

        # Write the header as a text file
        header.totextfile(path)

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
