#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import Astromagic modules
from .inpaint import replace_nans

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
    #data_ma = np.ma.array(data.astype(float), mask=mask)
    #data_nans = data_ma.filled(np.NaN)

    data_with_nans = np.copy(data)
    data_with_nans[mask] = np.NaN

    interpolated = replace_nans(data_with_nans, 5, 0.5, 2, "localmean")

    # If the interpolated box contains nans, do not fill in the corresponding pixels of the data with these nans,
    # therefore set the pixels that are nan to False in the box_mask (take the difference between the box_mask
    # and the np.isnan(interpolated_box) mask). Then, set the nans to zero in the interpolated_box because
    # False * nan would otherwise still equal to nan.
    #box_mask = masks.subtract(box_mask, np.isnan(interpolated_box))
    #interpolated_box[np.isnan(interpolated_box)] = 0.0

    interpolated[np.isnan(interpolated)] = data[np.isnan(interpolated)]

    # Return the interpolated data
    return interpolated

# *****************************************************************

def low_res_interpolation(frame, downsample_factor, mask=None):

    """
    This function ...
    :param frame:
    :param downsample_factor:
    :param mask:
    :return:
    """

    # Calculate the x and y size of the low-resolution map of the data
    low_res_x_size = int(round(frame.xsize / downsample_factor))
    low_res_y_size = int(round(frame.ysize / downsample_factor))

    # Create the background by interpolation
    back = Background(np.asarray(frame), (low_res_y_size, low_res_x_size), filter_shape=(3, 3), filter_threshold=None, mask=mask,
                      method='sextractor', backfunc=None, interp_order=3, sigclip_sigma=3.0, sigclip_iters=10)

    # Return the background frame
    #return Frame(back.background, frame.wcs, frame.pixelscale, frame.description, frame.selected, frame.unit)

    # Return the data
    return back

# *****************************************************************