#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import copy
import numpy as np

# Import astronomical modules
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy import log

# Import image modules
import tools.general

# *****************************************************************

def sigma_clip_mask_list(data, sigma=3.0, mask=None):

    """
    This function ...
    :param data:
    :param sigma:
    :param mask:
    :return:
    """

    masked_list = sigma_clip(data, sig=sigma, iters=None, copy=False)

    new_mask = copy.deepcopy(mask) if mask is not None else [0]*len(data)

    for i, masked in enumerate(masked_list.mask):

        if masked: new_mask[i] = True

    # Return the new or updated mask
    return new_mask

# *****************************************************************

def sigma_clip_mask(data, sigma=3.0, mask=None):

    """
    This function ...
    :param data:
    :param sigma:
    :param mask:
    :return:
    """

    # Split the x, y and z values of the data, without the masked values
    x_values, y_values, z_values = tools.general.split_xyz(data, mask=mask)

    # Sigma-clip z-values that are outliers
    masked_z_values = sigma_clip(z_values, sig=sigma, iters=None, copy=False)

    # Copy the mask or create a new one if none was provided
    new_mask = np.copy(mask) if mask is not None else np.zeros_like(data, dtype=np.bool)

    for i, masked in enumerate(masked_z_values.mask):

        if masked:

            x = x_values[i]
            y = y_values[i]
            new_mask[y,x] = True

    # Return the new or updated mask
    return new_mask

# *****************************************************************

def sigma_clipped_median(data, sigma=3.0, mask=None):

    """
    This function ...
    :param data:
    :param sigma:
    :param mask:
    :return:
    """

    # Calculate the sigma-clipped mean and median
    _, median, _ = sigma_clipped_stats(data, mask=mask, sigma=sigma)

    # Return the median value
    return median

# *****************************************************************

def sigma_clipped_statistics(data, sigma=3.0, mask=None):

    """
    This function ...
    :param data:
    :param sigma:
    :param mask:
    :return:
    """

    # Calculate the sigma-clipped mean and median
    mean, median, stddev = sigma_clipped_stats(data, mask=mask, sigma=sigma)

    # Return the statistical parameters
    return mean, median, stddev

# *****************************************************************

def sigma_clip_split(input_list, criterium):

    """
    This function ...
    :param input_list:
    :param criterium:
    :return:
    """

    # Initialize an empty list of widths
    determinants = []

    # Inform the user
    log.info("Seperating stars and unidentified objects...")

    # Loop over all the star candidates and calculate their width
    for item in input_list: determinants.append(criterium(item))

    # Use sigma clipping to seperate stars and unidentified objects
    mask = sigma_clip_mask_list(determinants, sigma=3.0)

    # Create a seperate list for the stars and for the ufos
    valid_list = []
    invalid_list = []

    # Loop over all items in the input list, putting them in either the valid or invalid list
    for index, src in enumerate(input_list):

        if mask[index]: invalid_list.append(src)
        else: valid_list.append(src)

    # Return the valid and invalid lists
    return valid_list, invalid_list

# *****************************************************************