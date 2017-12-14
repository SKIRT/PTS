#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.nr Numerical recipes like in SKIRT.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# -----------------------------------------------------------------

def find_nearest(array, value):

    """
    This function ...
    :param array:
    :param value:
    :return:
    """

    idx = (np.abs(array-value)).argmin()
    return array[idx]

# -----------------------------------------------------------------

def locate_clip(array, value):

    """
    This function ...
    :param array: actually a list
    :param value:
    :return:
    """

    #n = array.size
    n = len(array)
    if value < array[0]: return 0
    return locate_basic_impl(array, value, n-1)

# -----------------------------------------------------------------

def locate_basic_impl(xv, x, n):

    """
    This function ...
    :param xv:
    :param x:
    :param n:
    :return:
    """

    jl = -1
    ju = n

    while ju - jl > 1:

        jm = (ju + jl) >> 1
        if x < xv[jm]: ju = jm
        else: jl = jm

    return jl

# -----------------------------------------------------------------

def locate(xv, x):

    """
    This function ...
    :param xv:
    :param x:
    :return:
    """

    n = len(xv)
    if (x == xv[-1]): return n-2
    else: return locate_basic_impl(xv, x, n)

# -----------------------------------------------------------------

# This function resamples the function values \f$y_k\f$ defined on a grid \f$x_k\f$
# (both specified as arrays of the same length) onto an new grid \f$x_l^*\f$. The result is
# returned as an array of function values \f$y_l^*\f$ with the same length as the target
# grid. For new grid points that fall beyond the original grid, the function value is set to
# zero. For new grid points inside the original grid, the function value is interpolated
# using the function specified as template argument. The interpolation functions provided by
# this class can be passed as a template argument for this purpose. */
def resample_log_log(xresv, xoriv, yoriv):

    """
    This function ...
    :param xresv:
    :param xoriv:
    :param yoriv:
    :return:
    """

    Nori = len(xoriv)
    xmin = xoriv[0]
    xmax = xoriv[Nori-1]
    Nres = len(xresv)
    yresv = [None] * Nres

    # Loop over
    for l in range(Nres):

        x = xresv[l]
        if abs(1.0 - x/xmin) < 1e-5: yresv[l] = yoriv[0]
        elif abs(1.0 - x/xmax) < 1e-5: yresv[l] = yoriv[Nori-1]
        elif x < xmin or x > xmax: yresv[l] = 0.0
        else:

            k = locate(xoriv,x)
            yresv[l] = interpolate_log_log(x, xoriv[k], xoriv[k+1], yoriv[k], yoriv[k+1])

    # Return the interpolated values
    return yresv

# -----------------------------------------------------------------

def interpolate_log_log(x, x1, x2, f1, f2):

    """
    This function ...
    :param x:
    :param x1:
    :param x2:
    :param f1:
    :param f2:
    :return:
    """

    # Compute logarithm of coordinate values
    x  = np.log10(x)
    x1 = np.log10(x1)
    x2 = np.log10(x2)

    # Turn off logarithmic interpolation of function value if not all given values are positive
    logf = f1 > 0 and f2 > 0

    # Compute logarithm of function values if required
    if (logf):

        f1 = np.log10(f1)
        f2 = np.log10(f2)

    # Perform the interpolation
    fx = f1 + ((x-x1) / (x2-x1)) * (f2-f1)

    # Compute the inverse logarithm of the resulting function value if required
    if logf: fx = pow(10,fx)

    # Return the interpolated value
    return fx

# -----------------------------------------------------------------
