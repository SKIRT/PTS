#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.statistics Provides statistical functions.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import math
import numpy as np

# Import astronomical modules
from astropy.stats import sigma_clip, sigma_clipped_stats

# Import the relevant PTS classes and modules
from . import general
from ..basics.mask import Mask

# -----------------------------------------------------------------

# Calculate sigma-to-FWHM and FWHM-to-sigma conversion factors
sigma_to_fwhm = (8 * np.log(2))**0.5
fwhm_to_sigma = 1.0 / sigma_to_fwhm

# -----------------------------------------------------------------

def histogram(data, nbins):

    """
    This function ...
    :param data:
    :param nbins:
    :return:
    """

    # Get the bins
    freq, edges = np.histogram(data, nbins)
    nbins = len(edges) - 1

    # Get the bin centers
    centers = []
    for i in range(nbins): centers.append(0.5 * (edges[i] + edges[i + 1]))

    # Get the bin widths
    # widths = []
    # for i in range(len(edges) - 1):
    #    widths.append(edges[i + 1] - edges[i])

    # Get the lower limits of the bins
    lower = []
    for i in range(nbins): lower.append(edges[i])

    # Get the upper limits of the bins
    upper = []
    for i in range(nbins): upper.append(edges[i + 1])

    # Return the bin properties
    return centers, lower, upper

# -----------------------------------------------------------------

def sigma_clip_mask_list(data, sigma=3.0, mask=None):

    """
    This function ...
    :param data:
    :param sigma:
    :param mask:
    :return:
    """

    masked_list = sigma_clip(data, sigma=sigma, iters=None, copy=False)
    new_mask = copy.deepcopy(mask) if mask is not None else [0]*len(data)
    for i, masked in enumerate(masked_list.mask):
        if masked: new_mask[i] = True

    # Return the new or updated mask
    return new_mask

# -----------------------------------------------------------------

def sigma_clip_mask(data, sigma_level=3.0, mask=None, niters=None):

    """
    This function ...
    :param data:
    :param sigma_level:
    :param mask:
    :param niters: None means till convergence is achieved
    :return:
    """

    # Split the x, y and z values of the data, without the masked values
    x_values, y_values, z_values = general.split_xyz(data, mask=mask)

    # Sigma-clip z-values that are outliers
    masked_z_values = sigma_clip(z_values, sigma=sigma_level, iters=niters, copy=False)

    # Copy the mask or create a new one if none was provided
    new_mask = copy.deepcopy(mask) if mask is not None else Mask(np.zeros_like(data))

    for i, masked in enumerate(masked_z_values.mask):

        if masked:

            x = x_values[i]
            y = y_values[i]
            new_mask[y,x] = True

    #if not isinstance(new_mask, Mask): print(new_mask, mask)

    # Assert the mask is of type 'Mask'
    from ..core.mask import Mask as newMask
    assert isinstance(new_mask, Mask) or isinstance(new_mask, newMask)

    # Return the new or updated mask
    return new_mask

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

def sigma_clip_split(input_list, criterion, sigma=3.0, only_high=False, only_low=False, nans="low"):

    """
    This function ...
    :param input_list:
    :param criterion:
    :param sigma:
    :param only_high:
    :param only_low:
    :param nans:
    :return:
    """

    # Initialize an empty list of widths
    determinants = []

    # Loop over all the star candidates and calculate their width
    for item in input_list: determinants.append(criterion(item))

    # Use sigma clipping to seperate stars and unidentified objects
    mask = sigma_clip_mask_list(determinants, sigma=sigma)

    # Calculate the mean value of the determinants that are not masked
    mean = np.ma.mean(np.ma.masked_array(determinants, mask=mask))

    # Create a seperate list for the stars and for the ufos
    valid_list = []
    invalid_list = []

    # Loop over all items in the input list, putting them in either the valid or invalid list
    for index, item in enumerate(input_list):

        value = criterion(item)

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

# -----------------------------------------------------------------

def cutoff(values, method, limit):

    """
    This function ...
    :param values:
    :param method:
    :param limit:
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
        masked_values = sigma_clip(np.array(values), sigma=limit, iters=None, copy=False)

        # Calculate the maximum of the masked array
        return np.ma.max(masked_values)

    else: raise ValueError("Invalid cutoff method (must be 'percentage' or 'sigma_clip'")

# -----------------------------------------------------------------

def inverse_gaussian(center, sigma, amplitude, y):

    """
    This function ...
    :param center:
    :param sigma:
    :param amplitude:
    :param y:
    :return:
    """

    x = np.sqrt( - 2. * sigma**2 * np.log(y * np.sqrt(2. * math.pi * sigma**2) / amplitude)) + center
    return x

# -----------------------------------------------------------------

def gaussian(center, sigma, amplitude, x):

    """
    This function ...
    :param center:
    :param sigma:
    :param amplitude:
    :param x:
    :return:
    """

    normal = 1. / np.sqrt(2. * math.pi * sigma**2 ) * np.exp( - ( x - center)**2 / (2. * sigma**2))
    return normal * amplitude

# -----------------------------------------------------------------

def test_inverse(x, center, sigma, amplitude):

    """
    This function ...
    :param x:
    :param center:
    :param sigma:
    :param amplitude:
    :return:
    """

    y = gaussian(center, sigma, amplitude, x)
    return inverse_gaussian(center, sigma, amplitude, y)

# -----------------------------------------------------------------
