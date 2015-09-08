#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import numpy as np

# Import image modules
import regions

# *****************************************************************

class Mask(object):

    """
    This class ...
    """

    # *****************************************************************

    def __init__(self, data):

        """
        The constructor ...
        :param data:
        :return:
        """

        # Set the data array
        self.data = data

        # Set as unactive initially
        self.selected = False

    # *****************************************************************

    def select(self):

        """
        This function ...
        :return:
        """

        self.selected = True

    # *****************************************************************

    def deselect(self):

        """
        This function ...
        :return:
        """

        self.selected = False

# *****************************************************************

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

# *****************************************************************

def annulus_around(shape, inner_factor, outer_factor, x_size, y_size):

    # TODO: fix this for ellipses

    # Create new regions for the background estimation around the stars
    inner_shape = regions.scale_circle(shape, inner_factor)
    outer_shape = regions.scale_circle(shape, outer_factor)

    # Determine ...
    #x_center = outer_shape.coord_list[0]
    #y_center = outer_shape.coord_list[1]
    #outer_radius = outer_shape.coord_list[2]

    #y_min = int(round(y_center - outer_radius))
    #y_max = int(round(y_center + outer_radius))
    #x_min = int(round(x_center - outer_radius))
    #x_max = int(round(x_center + outer_radius))

    #smaller_x_size = x_max - x_min
    #smaller_y_size = y_max - y_min

    # Create inner and outer masks
    inner_mask = create_disk_mask(x_size, y_size, inner_shape.coord_list[0], inner_shape.coord_list[1], inner_shape.coord_list[2])
    outer_mask = create_disk_mask(x_size, y_size, outer_shape.coord_list[0], outer_shape.coord_list[1], outer_shape.coord_list[2])

    # Create the mask
    mask = inner_mask | np.logical_not(outer_mask)

    # Return the mask
    return mask

# *****************************************************************

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

# *****************************************************************

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

# *****************************************************************