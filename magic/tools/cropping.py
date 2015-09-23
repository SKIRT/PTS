#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

import numpy as np

# *****************************************************************

def crop(data, x_center, y_center, x_radius, y_radius):

    """
    This function ...
    :param data:
    :param x_center:
    :param y_center:
    :param x_radius:
    :param y_radius:
    :return:
    """

    # Determine y and x min and max
    y_min = int(round(y_center - y_radius))
    y_max = int(round(y_center + y_radius))
    x_min = int(round(x_center - x_radius))
    x_max = int(round(x_center + x_radius))

    # Return the cropped data
    return crop_direct(data, x_min, x_max, y_min, y_max)

# *****************************************************************

def crop_direct(data, x_min, x_max, y_min, y_max):

    """
    This function ...
    :param data:
    :param x_min:
    :param x_max:
    :param y_min:
    :param y_max:
    :return:
    """

    if y_min < 0: y_min = 0
    if x_min < 0: x_min = 0
    if y_max >= data.shape[0]: y_max = data.shape[0]
    if x_max >= data.shape[1]: x_max = data.shape[1]

    box = np.copy(data[y_min:y_max, x_min:x_max])

    # Return the new image, and the coordinates of the origin of the new image with respect to the original image
    return box, x_min, x_max, y_min, y_max

# *****************************************************************

def crop_check(data, x_min, x_max, y_min, y_max):

    """
    This function ...
    :param data:
    :param x_min:
    :param x_max:
    :param y_min:
    :param y_max:
    :return:
    """

    if y_min < 0: raise ValueError("y_min < 0")
    if x_min < 0: raise ValueError("x_min < 0")
    if y_max > data.shape[0]: raise ValueError("y_max >= data.ysize")
    if x_max > data.shape[1]: raise ValueError("x_max >= data.xsize")

    return crop_direct(data, x_min, x_max, y_min, y_max)[0]

# *****************************************************************
