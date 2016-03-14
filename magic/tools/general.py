#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
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

# From http://stackoverflow.com/questions/398299/looping-in-a-spiral
# Testing with:
# for a,b in spiral(3,3): print (a,b), you get (0, 0) (1, 0) (1, 1) (0, 1) (-1, 1) (-1, 0) (-1, -1) (0, -1) (1, -1)
# for a,b in spiral(5,3): print (a,b), you get (0, 0) (1, 0) (1, 1) (0, 1) (-1, 1) (-1, 0) (-1, -1) (0, -1) (1, -1) (2, -1) (2, 0) (2, 1) (-2, 1) (-2, 0) (-2, -1)
def spiral(N, M):
    x,y = 0,0
    dx, dy = 0, -1

    for dumb in xrange(N*M):
        if abs(x) == abs(y) and [dx,dy] != [1,0] or x>0 and y == 1-x:
            dx, dy = -dy, dx            # corner, change direction

        if abs(x)>N/2 or abs(y)>M/2:    # non-square
            dx, dy = -dy, dx            # change direction
            x, y = -y+dx, x+dy          # jump

        yield x, y
        x, y = x+dx, y+dy

# -----------------------------------------------------------------
