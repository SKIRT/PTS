#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.kernel Contains the ConvolutionKernel class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import ndimage

# Import the relevant PTS classes and modules
from .frame import Frame
from ...core.tools.logging import log
from ..tools import statistics

# -----------------------------------------------------------------

class ConvolutionKernel(Frame):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, data, *args, **kwargs):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(ConvolutionKernel, self).__init__(data, *args, **kwargs)

        # Check the FWHM
        if self.fwhm is None: raise ValueError("FWHM must be specified if not present in header")

        # Set the WCS to None, but keep the pixelscale
        if self._wcs is not None:

            if self._pixelscale is None: self._pixelscale = self.wcs.pixelscale
            elif not np.isclose(self._pixelscale, self.wcs.pixelscale): raise ValueError("Pixelscale in the header does not correspond to the specified pixelscale")
            self._wcs = None

        elif self._pixelscale is None: raise ValueError("Pixelscale must be specified if not present in header")

        # Normalized
        self._normalized = False

    # -----------------------------------------------------------------

    @property
    def normalized(self):

        """
        This function ...
        :return:
        """

        return self._normalized

    # -----------------------------------------------------------------

    def prepare_for(self, image, sigma_level=10.0):

        """
        This function ...
        :param image:
        :param sigma_level:
        :return:
        """

        # Prepare
        self.prepare(image.pixelscale, sigma_level)

    # -----------------------------------------------------------------

    def prepare(self, pixelscale, sigma_level=10.0):

        """
        This function ...
        :param pixelscale:
        :param sigma_level:
        :return:
        """

        # Truncate
        self.truncate(sigma_level)

        # Adjust pixelscale
        self.adjust_pixelscale(pixelscale)

        # Center
        self.center()

        # Normalize
        self.normalize()

    # -----------------------------------------------------------------

    def truncate(self, sigma_level=10.0):

        """
        This function ...
        :return:
        """

        sigma_pix = statistics.fwhm_to_sigma * self.fwhm_pix
        radius = sigma_level * sigma_pix

        center_x = 0.5 * (self.xsize - 1.)
        center_y = 0.5 * (self.ysize - 1.)

        min_x = int(round(center_x - radius))
        max_x = int(round(center_x + radius))

        min_y = int(round(center_y - radius))
        max_y = int(round(center_y + radius))

        # Crop
        self.crop(min_x, max_x, min_y, max_y)

    # -----------------------------------------------------------------

    def adjust_pixelscale(self, pixelscale):

        """
        This function ...
        :return:
        """

        average_pixelscale = 0.5 * (pixelscale.x + pixelscale.y)

        # Calculate the zooming factor
        factor = (average_pixelscale / self.xy_average_pixelscale).to("").value

        # Rebin to the pixelscale
        new_data = ndimage.interpolation.zoom(self._data, zoom=1.0 / factor)

        # Set the new data and pixelscale
        self._data = new_data
        self._pixelscale = pixelscale

    # -----------------------------------------------------------------

    def center(self):

        """
        This function ...
        :return:
        """

        # FROM CONVOLVE_IMAGE.PRO:

        # size_im      = (size(image))[1]
        # center_pixel = fix((size_im - 1) / 2)

        # get_maximun,image,x_max,y_max

        # ; determine the needed shifts
        # shift_x = center_pixel - x_max
        # shift_y = center_pixel - y_max

        # ; make the shift if nonzero
        # if (shift_x ne 0) or (shift_y ne 0) then begin
        #   if do_we_write eq 1 then print,'Shifting the PSF center by ('+ strtrim(string(shift_x),2) +','+ strtrim(string(shift_y),2) + ') pixels'
        #   image = shift(image,shift_x,shift_y)
        #   image[  0                     :abs(shift_x),*]=0.0
        #   image[  size_im-1-abs(shift_x):size_im-1   ,*]=0.0
        #   image[*,0                     :abs(shift_y)  ]=0.0
        #   image[*,size_im-1-abs(shift_y):size_im-1     ]=0.0
        # endif

        # ; We check that the centering is OK:
        # get_maximun,image,x_max,y_max
        # shift_x = center_pixel - x_max
        # shift_y = center_pixel - y_max
        # if (shift_x ne 0) or (shift_y ne 0) then if do_we_write eq 1 then print,'WARNING: Something went wrong in the image centering routine!!!'

        # if do_we_write eq 1 then print,'The PSF was centered successfully.'
        # if do_we_write eq 1 then print,' '

        pass

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        self.__idiv__(self.sum())
        self._normalized = True

    # -----------------------------------------------------------------

    def check_normalization(self):

        """
        This function ...
        :return:
        """

        assert np.abs(self.sum() - 1) < 1e-8 # same criterion as in astropy.convolution module

# -----------------------------------------------------------------
