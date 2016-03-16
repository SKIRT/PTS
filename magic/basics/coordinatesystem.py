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
import numpy as np

# Import astronomical modules
from astropy import wcs
from astropy.wcs import utils
from astropy.units import Unit
from astropy.io import fits
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .vector import Extent
from .geometry import Coordinate
from .skygeometry import SkyCoordinate
from ..tools import coordinates
from ...core.tools.logging import log

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

        # Call the constructor of the base class
        super(CoordinateSystem, self).__init__(header)

        # Set the number of pixels in both axis directions
        if header is not None:

            if not hasattr(self, "naxis1"): self.naxis1 = header["NAXIS1"]
            if not hasattr(self, "naxis2"): self.naxis2 = header["NAXIS2"]

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

        return cls(header)

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
        return Extent(x_pixelscale, y_pixelscale)

    # -----------------------------------------------------------------

    @property
    def xy_average_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = self.pixelscale

        x_pixelscale = abs(pixelscale.x.to("arcsec/pix"))
        y_pixelscale = abs(pixelscale.y.to("arcsec/pix"))

        if not np.isclose(x_pixelscale.value, y_pixelscale.value, rtol=0.0005):
            log.warning("Averaging the pixelscale over the x and y direction may not be a good approximation:")
            log.warning("  * x pixelscale (absolute value) = " + str(x_pixelscale))
            log.warning("  * y pixelscale (absolute value) = " + str(y_pixelscale))

        # Return a single value for the pixelscale in arcseconds
        return 0.5 * (x_pixelscale + y_pixelscale)

    # -----------------------------------------------------------------

    @property
    def center_pixel(self):

        """
        This property ...
        :return:
        """

        x, y = self.wcs.crpix
        return Coordinate(x, y)

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
        new_x_pixel = Coordinate(center.x + pix_x, center.y)
        new_y_pixel = Coordinate(center.x, center.y + pix_y)

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

        #return ra_alignment, dec_alignment

    # -----------------------------------------------------------------

    @property
    def orientation_angle(self):

        """
        This function ...
        :return:
        """

        # TODO: Should support any arbritrary orientation angle, not just 0, +90, +180, -90 ...

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
        else: raise ValueError("Unknown orientation angle")

    # -----------------------------------------------------------------

    def to_header(self, relax=None, key=None):

        """
        This function ...
        :param relax:
        :param key:
        :return:
        """

        header = super(CoordinateSystem, self).to_header(relax, key)

        header["NAXIS1"] = self._naxis1
        header["NAXIS2"] = self._naxis2

        return header

    # -----------------------------------------------------------------

    @property
    def csys(self):

        """
        This function ...
        :return:
        """

        ctype = self.wcs.ctype[0]
        if 'RA' in ctype or 'DEC' in ctype:
            if self.wcs.equinox == 2000:
                return 'fk5'
            elif self.wcs.equinox == 1950:
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
