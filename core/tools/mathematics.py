#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.mathematics Useful mathematical and statistical functions
#
# The functions in this module are used for ...

# -----------------------------------------------------------------

# Import standard modules
import numpy as np
import itertools
from scipy import optimize

# -----------------------------------------------------------------

# This function performs a two-dimensional polynomial fit.
# This function was based on code provided by Joe Kington, see: http://goo.gl/yXYrmo
def fit_polynomial(x, y, z, order, linear):

    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
        if linear & (i != 0.) & (j != 0.):
            G[:, k] = 0
    m, _, _, _ = np.linalg.lstsq(G, z)

    # Return the parameter list
    return m

# -----------------------------------------------------------------

# Values to two-dimensional polynomial fit.
# This function was based on code provided by Joe Kington, see: http://goo.gl/yXYrmo
def polynomial(x, y, m):

    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j

    # Return the value
    return z

# -----------------------------------------------------------------

# Returns (height, x, y, width_x, width_y) the gaussian parameters of a 2D distribution found by a fit
def fitgaussian(data):

    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)

    # Return the parameter list
    return p

# -----------------------------------------------------------------

# Returns a 2D Gaussian distribution form the given parameters.
def gaussian(height, center_x, center_y, width_x, width_y):

    width_x = float(width_x)
    width_y = float(width_y)

    # Return the distribution
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

# -----------------------------------------------------------------

# This function returns the gaussian parameters of a 2D distribution by calculating its moments
def moments(data):

    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()

    # Return the parameters
    return height, x, y, width_x, width_y

# -----------------------------------------------------------------

