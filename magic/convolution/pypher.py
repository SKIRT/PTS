#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.convolution.pypher Contains the PypherKernels class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .kernels import Kernels, kernels_path

# -----------------------------------------------------------------

# Reference: http://pypher.readthedocs.org/en/latest/

# -----------------------------------------------------------------

class PypherKernels(Kernels):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(PypherKernels, self).__init__()

        # The path to the directory where the generated Pypher kernels are saved (to be reused)
        self.path = fs.join(kernels_path, "pypher")

    # -----------------------------------------------------------------

    def get_kernel(self, from_filter, to_filter, high_res=True, from_fwhm=None, to_fwhm=None):

        """
        This function ...
        :param from_filter:
        :param to_filter:
        :param high_res:
        :param from_fwhm:
        :param to_fwhm:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_psf(self, fltr):

        """
        THis function ...
        :param fltr:
        :return:
        """

        pass

# -----------------------------------------------------------------
