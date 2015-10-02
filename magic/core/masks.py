#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np
from scipy import ndimage

# Import Astromagic modules
from . import regions
from .vector import Extent

# *****************************************************************

class Mask(np.ndarray):

    """
    This class ...
    """

    # *****************************************************************

    def __new__(cls, data, selected=False, description=None):

        """
        This function ...
        :param cls:
        :param data:
        :param selected:
        :param description:
        :return:
        """

        obj = np.asarray(data, dtype=bool).view(cls)
        obj.selected = selected
        obj.description = description

        return obj

    # *****************************************************************

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        if obj is None: return
        self.selected = getattr(obj, 'selected', False)
        self.description = getattr(obj, 'description', None)

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

    @property
    def xsize(self): return self.shape[1]

    # *****************************************************************

    @property
    def ysize(self): return self.shape[0]

    # *****************************************************************

    def expand(self, structure_size=2, iterations=100):

        """
        This function ...
        """

        # Define the structure for the expansion
        structure = ndimage.generate_binary_structure(structure_size, structure_size)

        # Make the new mask, made from 100 iterations with the structure array
        self = ndimage.binary_dilation(self, structure, iterations)

        # Check whether this object remains of type Mask
        assert isinstance(self, Mask)

    # *****************************************************************

    def apply(self, frame, fill=0.0):

        """
        This function ...
        """

        # Replace the masked pixel values with the given fill value
        frame[self] = fill

    # *****************************************************************

    def inverse(self):

        """
        This function ...
        """

        # Check whether the new mask will still be an instance of Mask
        assert isinstance(np.logical_not(self), Mask)

        # Return the inverse of this mask
        return np.logical_not(self)

    # *****************************************************************

    def intersection(self, mask):

        """
        This function ...
        :param mask_a:
        :param mask_b:
        :return:
        """

        assert isinstance(self * mask, Mask)

        return self * mask

    # *****************************************************************

    def subtract(self, mask):

        """
        This function ...
        :param mask_a:
        :param mask_b:
        :return:
        """

        assert isinstance(mask.inverse(), Mask)

        # Return ...
        return self.intersection(mask.inverse())

    # *****************************************************************

    def union(self, mask):

        """
        This function ...
        :param mask_a:
        :param mask_b:
        :return:
        """

        assert isinstance(self + mask, Mask)

        # Return ...
        return self + mask

    # *****************************************************************

    def hits_boundary(self, min_pixels=2):

        """
        This function ...
        """

        # Set the number of hits (pixels hitting the boundary) to zero initially
        hits = 0

        for x in range(self.xsize):

            if self[0, x] or self[self.ysize-1, x]:

                # Add one hit
                hits += 1

                # Break the loop, avoid performing unnecessary checks
                if hits >= min_pixels: break

        for y in range(1, self.ysize-1):

            if self[y, 0] or self[y, self.xsize-1]:

                # Add one hit
                hits += 1

                # Break the loop, avoid performing unnecessary checks
                if hits >= min_pixels: break

        # Return whether ...
        return hits >= min_pixels

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

def create_ellipse_mask(x_size, y_size, center, radius, angle):

    """
    This function ...
    :param x_size:
    :param y_size:
    :param x_center:
    :param y_center:
    :param x_radius:
    :param y_radius:
    :param angle:
    :return:
    """

    x_radius = radius.x if isinstance(radius, Extent) else radius
    y_radius = radius.y if isinstance(radius, Extent) else radius

    # Create a region consisting of one ellipse
    region = regions.one_ellipse([center.x, center.y, x_radius, y_radius, angle])

    # Create the mask
    data = region.get_mask(shape=(y_size, x_size))

    # Return the mask
    return Mask(data)

# *****************************************************************

def create_annulus_mask(xsize, ysize, center, inner_radius, outer_radius, angle):

    """
    This function ...
    :param xsize:
    :param ysize:
    :param center:
    :param radius:
    :param angle:
    :return:
    """

    inner_mask = create_ellipse_mask(xsize, ysize, center, inner_radius, angle)
    outer_mask = create_ellipse_mask(xsize, ysize, center, outer_radius, angle)

    # Return the annulus mask
    return inner_mask.union(outer_mask.inverse())

# *****************************************************************

