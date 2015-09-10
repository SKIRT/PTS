#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import numpy as np

# Import image modules
from astromagic.inpaint import replace_nans

# Import astronomical modules
from photutils.background import Background

# *****************************************************************

def in_paint(data, mask):

    """
    This function ...
    :param data:
    :param mask:
    :return:
    """

    # Fill the data with nans according to the mask
    data_ma = np.ma.array(data.astype(float), mask=mask)
    data_nans = data_ma.filled(np.NaN)

    interpolated = replace_nans(data_nans, 5, 0.5, 2, "localmean")

    return interpolated

# *****************************************************************

def low_res_interpolation(data, downsample_factor, mask=None):

    # Calculate the x and y size of the low-resolution map of the data
    low_res_x_size = data.shape[1] / downsample_factor
    low_res_y_size = data.shape[0] / downsample_factor

    # Create the background by interpolation
    back = Background(data, (low_res_x_size, low_res_y_size), filter_shape=(3, 3), filter_threshold=None, mask=mask,
                      method='sextractor', backfunc=None, interp_order=3, sigclip_sigma=3.0, sigclip_iters=10)

    # Return the background data
    return back.background

# *****************************************************************