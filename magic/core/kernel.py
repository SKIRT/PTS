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
from .frame import Frame
from ...core.tools.logging import log

# -----------------------------------------------------------------

class ConvolutionKernel(Frame):

    """
    This class ...
    """

    def __new__(cls, data):

        """
        This function ...
        :param cls:
        :param data:
        :return:
        """

        # Create an object of the Box class
        obj = np.asarray(data).view(cls)

        # Return the object
        return obj

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, fwhm=None):

        """
        This function ...
        :param path:
        :param fwhm:
        :return:
        """

        # Open the file
        kernel = super(ConvolutionKernel, cls).from_file(path, no_filter=True)

        # Check the FWHM
        if kernel.fwhm is None:
            if fwhm is None: raise ValueError("FWHM must be specified if not present in header")
        elif fwhm is not None: assert kernel.fwhm == fwhm

        # Return the kernel
        return kernel

    # -----------------------------------------------------------------

    def truncate(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
