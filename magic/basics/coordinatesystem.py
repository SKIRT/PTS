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
from astropy.coordinates import SkyCoord

# Import the relevant PTS classes and modules
from .vector import Extent, Position
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
            log.warning("  * x pixelscale (absolute value) =", x_pixelscale)
            log.warning("  * y pixelscale (absolute value) =", y_pixelscale)

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
        return Position(x, y)

    # -----------------------------------------------------------------

    @property
    def center_sky(self):

        """
        This property ...
        :return:
        """

        center = self.center_pixel
        return SkyCoord.from_pixel(center.x, center.y, self, mode='wcs')

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

        # Get the center pixel of the coordinate system
        center = self.center_pixel
        center_coordinate = self.center_sky

        # Calculate the new pixel position
        new_position = Position(center.x + pix_x, center.y + pix_y)

        # Convert to sky coordinate
        new_coordinate = self.to_sky(new_position)

        print("center coordinate:", center_coordinate)
        print("new coordinate:", new_coordinate)

        return None

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
