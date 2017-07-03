#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.convolution.matching Contains the MatchingKernels class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from photutils import create_matching_kernel
from photutils import TopHatWindow
from astropy.modeling.models import Gaussian2D

# Import the relevant PTS classes and modules
from ..core.kernel import ConvolutionKernel
from .kernels import Kernels, get_fwhm

# -----------------------------------------------------------------

class MatchingKernels(Kernels):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(MatchingKernels, self).__init__()

    # -----------------------------------------------------------------

    def get_kernel(self, from_filter, to_filter, pixelscale, from_fwhm=None, to_fwhm=None, from_model="gaussian",
                   to_model="gaussian", size=51):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param pixelscale:
        :param from_fwhm:
        :param to_fwhm:
        :param from_model:
        :param to_model:
        :param size:
        :return:
        """

        # Make sure from_fwhm and to_fwhm is defined
        if from_fwhm is None: from_fwhm = get_fwhm(from_filter)
        if to_fwhm is None: to_fwhm = get_fwhm(to_filter)

        # Make sure that size is odd number
        if size % 2 == 0: size += 1

        # Determine center based on size
        center = int((size - 1) / 2)

        # Create window
        window = TopHatWindow(0.35)

        y, x = np.mgrid[0:size, 0:size]

        # Determine FWHMs in number of pixels
        from_fwhm_pix = from_fwhm / pixelscale.average
        to_fwhm_pix = to_fwhm / pixelscale.average

        # amplitude, x_mean, y_mean, x_stddev, y_stddev
        gm1 = Gaussian2D(1, center, center, from_fwhm_pix, from_fwhm_pix)
        gm2 = Gaussian2D(1, center, center, to_fwhm_pix, to_fwhm_pix)

        # Generate
        g1 = gm1(x, y)
        g2 = gm2(x, y)

        # Normalize
        g1 /= g1.sum()
        g2 /= g2.sum()

        # Create the kernel, set FWHM, from_filter, to_filter, prepared=True, and pixelscale
        data = create_matching_kernel(g1, g2, window=window)
        kernel = ConvolutionKernel(data, fwhm=to_fwhm, from_filter=from_filter, to_filter=to_filter, prepared=True, pixelscale=pixelscale)

        # Return the kernel
        return kernel

    # -----------------------------------------------------------------

    def get_psf(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        pass

# -----------------------------------------------------------------
