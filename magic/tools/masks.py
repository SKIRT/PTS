#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.masks Contains functions for dealing with two-dimensional masks.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..region import tools as regions

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

#def union(*args): # i wanted to do it this way, but didn't succeed ...
def union(mask_a, mask_b):

    """
    This function ...
    :param args:
    :return:
    """

    return mask_a + mask_b

# -----------------------------------------------------------------

#def intersection(*args): i wanted to do it this way, but didn't succeed ...
def intersection(mask_a, mask_b):

    """
    This function ...
    :param args:
    :return:
    """

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

def split_overlap(base_mask, test_mask, return_segments=False):

    """
    This function takes all blobs in the base_mask and checks whether they overlap with the test_mask.
    The function returns two new masks, one mask with all the blobs that overlapped, and another with the blobs
    that did not overlap.
    :param base_mask:
    :param test_mask:
    :return:
    """

    overlapping = np.zeros_like(base_mask, dtype=bool)
    not_overlapping = np.copy(base_mask)

    from photutils import detect_sources
    segments = detect_sources(base_mask.astype('float'), 0.5, 1).data
    overlap = intersection(segments, test_mask)

    # Check which indices are present in the overlap map
    possible = np.array(range(1, np.max(overlap) + 1))
    present = np.in1d(possible, overlap)
    indices = possible[present]

    overlapping_segments = np.zeros_like(base_mask, dtype=int)
    not_overlapping_segments = np.copy(segments)

    # Remove the galaxies from the segmentation map
    for index in indices:
        blob = segments == index
        overlapping[blob] = True
        not_overlapping[blob] = False

        overlapping_segments[blob] = index
        not_overlapping_segments[blob] = 0

    if return_segments: return overlapping, not_overlapping, overlapping_segments, not_overlapping_segments
    else: return overlapping, not_overlapping

# -----------------------------------------------------------------
