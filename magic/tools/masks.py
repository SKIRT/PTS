#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant AstroMagic classes and modules
from . import regions

# -----------------------------------------------------------------

def annuli_around(region, inner_factor, outer_factor, header, x_size, y_size):

    """
    This function ...
    :param region:
    :param inner_factor:
    :param outer_factor:
    :param header:
    :param x_size:
    :param y_size:
    :return:
    """

    # Create new regions for the background estimation around the stars
    inner_region = regions.expand(region, inner_factor)
    outer_region = regions.expand(region, outer_factor)

    # Create inner and outer masks
    inner_mask = regions.create_mask(inner_region, header, x_size, y_size)
    outer_mask = regions.create_mask(outer_region, header, x_size, y_size)

    # Create the mask
    mask = inner_mask | np.logical_not(outer_mask)

    # Return the mask
    return mask

# -----------------------------------------------------------------

def masked_outside(region, header, x_size, y_size, expand_factor=1.0):

    """
    This function ...
    :param region:
    :param header:
    :param x_size:
    :param y_size:
    :param expand_factor:
    :return:
    """

    # Create a new region ...
    region = regions.expand(region, factor=expand_factor)

    # Create a mask from the region
    mask = np.logical_not(regions.create_mask(region, header, x_size, y_size))

    # Return the mask
    return mask

# -----------------------------------------------------------------

def create_disk_mask(x_size, y_size, x_center, y_center, radius):

    """
    This function ...
    :param x_size:
    :param y_size:
    :param x_center:
    :param y_center:
    :param radius:
    :return:
    """

    # Calculate which pixels should be masked
    y,x = np.ogrid[-y_center:y_size-y_center, -x_center:x_size-x_center]
    mask = x*x + y*y <= radius*radius

    # Return the mask
    return mask

# -----------------------------------------------------------------

def union(mask_a, mask_b):

    """
    This function ...
    :param mask_a:
    :param mask_b:
    :return:
    """

    # Return the unified mask
    return mask_a + mask_b

# -----------------------------------------------------------------

def intersection(mask_a, mask_b):

    """
    This function ...
    :param mask_a:
    :param mask_b:
    :return:
    """

    # Return the intersection
    return mask_a * mask_b

# -----------------------------------------------------------------

def overlap(mask_a, mask_b):

    """
    This function ...
    :param mask_a:
    :param mask_b:
    :return:
    """

    return np.any(intersection(mask_a, mask_b))

# -----------------------------------------------------------------
