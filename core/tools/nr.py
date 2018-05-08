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

    idx = find_nearest_index(array, value)
    return array[idx]

# -----------------------------------------------------------------

def find_nearest_index(array, value):

    """
    This function ...
    :param array:
    :param value:
    :return:
    """

    return (np.abs(array-value)).argmin()

# -----------------------------------------------------------------

def find_bounds(array, value):

    """
    This function ...
    :param array:
    :param value:
    :return:
    """

    indices = locate_bounds(array, value)
    if isinstance(indices, tuple):
        lower = indices[0]
        upper = indices[1]
        if lower == -1: return None, array[upper]
        if upper == len(array): return array[lower], None
    else: return array[indices]

# -----------------------------------------------------------------

def locate_bounds(array, value):

    """
    This function ...
    :param array:
    :param value:
    :return:
    """

    nvalues = len(array)
    indices = np.argsort(array)

    previous_index = -1
    for i in range(nvalues):

        index = indices[i]
        value_i = array[index]
        if value_i == value: return index
        elif value_i > value: return (previous_index, index)
        previous_index = index

    # All values are smaller than the target value
    return (previous_index, previous_index+1)

# -----------------------------------------------------------------

def locate_continuous(array, value, scale="linear", out_of_bounds=None):

    """
    This function ...
    :param array:
    :param value:
    :param scale:
    :param out_of_bounds: None, error, or clip
    :return:
    """

    #print(array, value)

    indices = locate_bounds(array, value)

    if isinstance(indices, tuple):

        lower = indices[0]
        upper = indices[1]

        if lower == -1:
            if out_of_bounds is None: return None
            elif out_of_bounds == "error": raise ValueError("Out of bounds")
            elif out_of_bounds == "clip": return float(upper)
            else: raise ValueError("Invalid option for 'out_of_bounds'")
        if upper == len(array):
            if out_of_bounds is None: return None
            elif out_of_bounds == "error": raise ValueError("Out of bounds")
            elif out_of_bounds == "clip": return float(lower)
            else: raise ValueError("Invalid option for 'out_of_bounds'")

        #print(upper, lower)
        if upper == lower + 1:

            x1 = float(lower)
            x2 = float(upper)

            # Compute logarithm of function values if required
            if scale == "linear":

                fx = value
                f1 = array[lower]
                f2 = array[upper]

            elif scale == "logarithmic":

                fx = np.log10(value)
                f1 = np.log10(array[lower])
                f2 = np.log10(array[upper])

            else: raise ValueError("Invalid value for 'scale'")

            # Return
            return interpolate_inverse(fx, x1, x2, f1, f2)

        elif upper == lower - 1: raise NotImplementedError("Not yet implemented")
        else: raise ValueError("Values are not sequential")

    else: return float(indices)

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

def interpolate(x, x1, x2, f1, f2):

    """
    This function ...
    :param x:
    :param x1:
    :param x2:
    :param f1:
    :param f2:
    :return:
    """

    # Perform the interpolation
    fx = f1 + ((x - x1) / (x2 - x1)) * (f2 - f1)

    # Return
    return fx

# -----------------------------------------------------------------

def interpolate_inverse(fx, x1, x2, f1, f2):

    """
    This function ...
    :param fx:
    :param x1:
    :param x2:
    :param f1:
    :param f2:
    :return:
    """

    # Return
    x = x1 + (fx - f1) / (f2 - f1) * (x2 - x1)

    # Return
    return x

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
