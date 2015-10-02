#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import copy
import numpy as np

# Import astronomical modules
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy import log
import pyregion

# Import Astromagic modules
from . import general

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
    x_values, y_values, z_values = general.split_xyz(data, mask=mask)

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

    #data = np.asarray(data).view(np.ndarray)

    #print type(data), type(mask), sigma

    # Calculate the sigma-clipped mean and median
    mean, median, stddev = sigma_clipped_stats(data, mask=mask, sigma=sigma)

    # Return the statistical parameters
    return mean, median, stddev

# *****************************************************************

def sigma_clip_split(input_list, criterium, sigma=3.0, only_high=False, only_low=False, nans="low"):

    """
    This function ...
    :param input_list:
    :param criterium:
    :return:
    """

    # Initialize an empty list of widths
    determinants = []

    # Loop over all the star candidates and calculate their width
    for item in input_list: determinants.append(criterium(item))

    # Use sigma clipping to seperate stars and unidentified objects
    mask = sigma_clip_mask_list(determinants, sigma=sigma)

    # Calculate the mean value of the determinants that are not masked
    mean = np.ma.mean(np.ma.masked_array(determinants, mask=mask))

    # Create a seperate list for the stars and for the ufos
    valid_list = []
    invalid_list = []

    # Loop over all items in the input list, putting them in either the valid or invalid list
    for index, item in enumerate(input_list):

        value = criterium(item)

        if only_high:

            if mask[index] and value > mean: invalid_list.append(item)
            else: valid_list.append(item)

        elif only_low:

            if mask[index] and value < mean: invalid_list.append(item)
            else: valid_list.append(item)

        else:

            if mask[index]: invalid_list.append(item)
            else: valid_list.append(item)

    # Return the valid and invalid lists
    return valid_list, invalid_list

# *****************************************************************

def split_percentage(input_list, criterium, percentage, nans="low"):

    """
    This function...
    """

    if type(input_list).__name__ == "ShapeList":

        new_list = pyregion.ShapeList([])
        nans_list = pyregion.ShapeList([])

    else:

        new_list = []
        nans_list = []

    # Fill in new and nans lists
    for item in input_list:

        if np.isnan(criterium(item)): nans_list.append(item)
        else: new_list.append(item)

    # Do the sorting
    new_list.sort(key=criterium)

    # Determine the splitting point
    split = int(round(len(new_list) - percentage * len(new_list)))

    # Return the two splitted lists
    if type(input_list).__name__ == "ShapeList":

        list_a = pyregion.ShapeList([])
        list_b = pyregion.ShapeList([])

        for i in range(len(new_list)):

            if i < split: list_a.append(new_list[i])
            else: list_b.append(new_list[i])

        if nans == "low":
            for j in range(len(nans_list)): list_a.append(nans_list[j])
        elif nans == "high":
            for j in range(len(nans_list)): list_b.append(nans_list[j])
        elif nans == "discard": pass
        else: raise ValueError("Invalid option for nan")

        return list_a, list_b

    else:

        if nans == "low": return new_list[0:split-1] + nans_list, new_list[split:len(new_list)]
        elif nans == "high": return new_list[0:split-1], new_list[split:len(new_list)] + nans_list
        elif nans == "discard": return new_list[0:split-1], new_list[split:len(new_list)]
        else: raise ValueError("Invalid option for nan")

# *****************************************************************

def cutoff(values, method, limit):

    """
    This function ...
    """

    # Percentage method
    if method == "percentage":

        # Create a sorted list for the input values
        sorted_values = sorted(values)

        # Determine the splitting point
        split = (1.0-limit) * len(sorted_values)
        index = int(round(split))

        # Return the corresponding value in the sorted list
        return sorted_values[index]

    # Sigma-clipping method
    elif method == "sigma_clip":

        # Perform sigma clipping on the input list
        masked_values = sigma_clip(np.array(values), sig=limit, iters=None, copy=False)

        # Calculate the maximum of the masked array
        return np.ma.max(masked_values)

    else: raise ValueError("Invalid cutoff method (must be 'percentage' or 'sigma_clip'")

# *****************************************************************