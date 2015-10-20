#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import math
import numpy as np
from scipy import ndimage
from skimage import morphology

# Import Astromagic modules
from . import regions
from .regions import Region
from .vector import Position, Extent

# Import astronomical modules
import astropy.units as u
from astropy.coordinates import Angle
from photutils import detect_sources

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

    @classmethod
    def from_region(cls, region, shape):

        """
        This function ...
        :param region:
        :return:
        """

        # Return a new Mask object
        return cls(region.get_mask(shape=shape))

    # *****************************************************************

    @classmethod
    def from_ellipse(cls, x_size, y_size, center, radius, angle):

        """
        This function ...
        :param x_size:
        :param y_size:
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        # Create a region consisting of one ellipse
        region = Region.ellipse(center, radius, angle)

        # Create the mask
        data = region.get_mask(shape=(y_size, x_size))

        # Return a new Mask object
        return cls(data)

    # *****************************************************************

    @classmethod
    def from_aperture(cls, x_size, y_size, aperture, expansion_factor=1.0):

        """
        This function ...
        :param aperture:
        :return:
        """

        # TODO: make this work with apertures other than EllipticalAperture

        # Get the parameters of the elliptical aperture
        x_center, y_center = aperture.positions[0]
        center = Position(x=x_center, y=y_center)

        major = aperture.a * expansion_factor
        minor = aperture.b * expansion_factor

        radius = Extent(x=major, y=minor)

        # theta is in radians
        angle = Angle(aperture.theta, u.rad)

        # Return a new Mask object
        return cls.from_ellipse(x_size, y_size, center, radius, angle)

    # *****************************************************************

    def dilated(self, structure=None, connectivity=2, iterations=100):

        """
        This function ...
        :param connectivity:
        :param iterations:
        :return:
        """

        # Define the structure for the expansion
        if structure is None: structure = ndimage.generate_binary_structure(2, connectivity=connectivity)

        # Make the new mask, made from 100 iterations with the structure array
        data = ndimage.binary_dilation(self, structure, iterations)

        # Reassign this object
        return Mask(data, self.selected, self.description)

    # *****************************************************************

    def eroded(self, structure=None, connectivity=2, iterations=100):

        """
        This function ...
        :param connectivity:
        :param iterations:
        :return:
        """

        # Define the structure for the expansion
        if structure is None: structure = ndimage.generate_binary_structure(2, connectivity=connectivity)

        # Make the new mask, made from 100 iterations with the structure array
        data = ndimage.binary_erosion(self, structure, iterations)

        # Reassign this object
        return Mask(data, self.selected, self.description)

    # *****************************************************************

    def opening(self, structure, iterations=1):

        """
        This function ...
        :return:
        """

        data = ndimage.binary_opening(self, structure, iterations)

        # Return the new mask
        return Mask(data)

    # *****************************************************************

    def expanded(self, factor):

        """
        This function ...
        :return:
        """

        # TODO: use radius from aperture of mask (application: saturation)?
        radius = 0.5*self.xsize

        #structure = np.array([[False, True, True, True, False],
        #                      [True, True, True, True, True],
        #                      [True, True, True, True, True],
        #                      [True, True, True, True, True],
        #                      [False, True, True, True, False]])


        extra_pixels = radius * (factor - 1.0)

        iterations = int(round(extra_pixels/4.0))

        print("Factor = " + str(factor))
        print("number of iterations = " + str(iterations))

        disk_structure = morphology.disk(4)

        if iterations < 1: return self
        else: return self.dilated(disk_structure, iterations=iterations)

    # *****************************************************************

    def remove_appendages(self):

        """
        This function ...
        :return:
        """

        structure = np.array([[False, True, True, True, False],
                              [True, True, True, True, True],
                              [True, True, True, True, True],
                              [True, True, True, True, True],
                              [False, True, True, True, False]])

        mask = self.opening(structure)

        segments = detect_sources(mask, 0.5, 1)

        # Get the label of the center segment
        label = segments[int(0.5*segments.shape[0]), int(0.5*segments.shape[1])]

        # Return the new mask with the appendages removed
        return Mask((segments == label))

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

    def masks(self, position):

        """
        This function ...
        """

        # Calculate x and y of the pixel corresponding to the object's position
        x = int(round(position.x))
        y = int(round(position.y))

        # Return the value of the mask in this pixel
        return self[y, x]

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

    # Create a region consisting of one ellipse
    region = Region.ellipse(center, radius, angle)

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

def union(mask_a, mask_b):

    """
    This function ...
    :param mask_a:
    :param mask_b:
    :return:
    """

    # Return the unified mask
    return mask_a + mask_b

# *****************************************************************

def intersection(mask_a, mask_b):

    """
    This function ...
    :param mask_a:
    :param mask_b:
    :return:
    """

    # Return the intersection
    return mask_a * mask_b

# *****************************************************************

def overlap(mask_a, mask_b):

    """
    This function ...
    :param mask_a:
    :param mask_b:
    :return:
    """

    return np.any(intersection(mask_a, mask_b))

# *****************************************************************
