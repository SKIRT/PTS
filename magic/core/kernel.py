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

# Import the relevant PTS classes and modules
from .frame import Frame, newFrame
from ...core.tools.logging import log
from ..tools import statistics

# -----------------------------------------------------------------

class ConvolutionKernel(newFrame):
#class ConvolutionKernel(Frame):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    #def __new__(cls, data, **kwargs):
    #def __new__(cls, data, **kwargs):

        #"""
        #This function ...
        #"""

        #obj = np.asarray(data).view(cls)

        #obj = np.ndarray.__new__(cls, data.shape, float, data)

        #for attr in vars(data):
        #    print(attr, getattr(data, attr))
        #    #print(vars(obj))
        #    if not hasattr(obj, attr):
        #        print("here")
        #        setattr(obj, attr, getattr(data, attr))

        #return obj

        #return super(ConvolutionKernel, cls).__new__(data, **kwargs)

    # -----------------------------------------------------------------

    #def __init__(self, data, **kwargs):

        #"""
        #This function ...
        #:param data:
        #:param kwargs:
        #"""

        # Call the constructor of the base class
        #super(ConvolutionKernel, self).__init__(data, **kwargs)
        #Frame.__init__(self, data, **kwargs)

        #for attr in vars(data):

            #setattr(self, attr, getattr(data, attr))

            #value = getattr(data, attr)
            #print(attr, value)
            #self.__dict__[attr] = value

    # -----------------------------------------------------------------

    @classmethod
    def read(cls, path, fwhm=None):

        """
        This function ...
        :param path:
        :param fwhm:
        :return:
        """

        # Open the file
        kernel = super(ConvolutionKernel, cls).read(path, no_filter=True)
        #kernel = cls(Frame.from_file(path, no_filter=True))

        # Check the FWHM
        if kernel.fwhm is None:
            if fwhm is None: raise ValueError("FWHM must be specified if not present in header")
        elif fwhm is not None: assert kernel.fwhm == fwhm

        # Return the kernel
        return kernel

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        self.truncate(sigma_level)

        self.center()

        self.normalize()

    # -----------------------------------------------------------------

    def truncate(self, sigma_level=5.0):

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

        self.center()

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
