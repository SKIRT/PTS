#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
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
from astropy import units as u
from astropy.io import fits

# Import the relevant AstroMagic classes and modules
from .vector import Extent

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

        x_pixelscale = result[0] * u.Unit("deg/pix")
        y_pixelscale = result[1] * u.Unit("deg/pix")

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
            print("WARNING: averaging the pixelscale over the x and y direction may not be a good approximation:")
            print("          x pixelscale (absolute value) =", x_pixelscale)
            print("          y pixelscale (absolute value) =", y_pixelscale)

        # Return a single value for the pixelscale in arcseconds
        return 0.5 * (x_pixelscale + y_pixelscale)

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

    def to_header(self, relax=None, key=None):

        """
        This function ...
        :return:
        """

        header = super(CoordinateSystem, self).to_header(relax, key)

        header["NAXIS1"] = self._naxis1
        header["NAXIS2"] = self._naxis2

        return header

    # -----------------------------------------------------------------

    def to_sky_coordinate(self, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def to_image_coordinate(self, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def to_sky_extent(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def to_image_extent(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        pass

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
