#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# *****************************************************************

## \package pts.magic.tools.general Contains general functions.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# -----------------------------------------------------------------

def split_xyz(data, mask=None, arrays=False):

    """
    This function ...
    :param data:
    :param mask:
    :param arrays:
    :return:
    """

    # Initialize lists to contain the x, y and z values
    x_values = []
    y_values = []
    z_values = []

    # Loop over all x and y values
    for x in range(data.shape[1]):
        for y in range(data.shape[0]):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:

                x_values.append(x)
                y_values.append(y)
                z_values.append(data[y,x])

    if arrays: return np.array(x_values), np.array(y_values), np.array(z_values)
    else: return x_values, y_values, z_values

# -----------------------------------------------------------------
