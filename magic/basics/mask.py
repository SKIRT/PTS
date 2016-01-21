#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.mask Contains the Mask class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import ndimage
from skimage import morphology

# Import astronomical modules
from photutils import detect_sources

# Function to create mask from ellipse
from photutils.geometry import elliptical_overlap_grid

# Import the relevant AstroMagic classes and modules
from .vector import Position, Extent

# -----------------------------------------------------------------

class Mask(np.ndarray):

    """
    This class ...
    """

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    @classmethod
    def is_nan(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(np.isnan(frame))

    # -----------------------------------------------------------------

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        if obj is None: return
        self.selected = getattr(obj, 'selected', False)
        self.description = getattr(obj, 'description', None)

    # -----------------------------------------------------------------

    def select(self):

        """
        This function ...
        :return:
        """

        self.selected = True

    # -----------------------------------------------------------------

    def deselect(self):

        """
        This function ...
        :return:
        """

        self.selected = False

    # -----------------------------------------------------------------

    @property
    def xsize(self): return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self): return self.shape[0]

    # -----------------------------------------------------------------

    @classmethod
    def from_shape(cls, shape):

        """
        This function ...
        :param shape:
        :return:
        """

        return cls(np.zeros(shape))

    # -----------------------------------------------------------------

    @classmethod
    def from_region(cls, region, shape):

        """
        This function ...
        :param region:
        :return:
        """

        # Return a new Mask object
        return cls(region.get_mask(shape=shape))

    # -----------------------------------------------------------------

    @classmethod
    def from_rectangle(cls, x_size, y_size, rectangle):

        """
        This function ...
        :param x_size:
        :param y_size:
        :param rectangle:
        :return:
        """

        data = np.zeros((y_size, x_size))

        # Convert into integers
        x_min = int(round(rectangle.x_min))
        x_max = int(round(rectangle.x_max))
        y_min = int(round(rectangle.y_min))
        y_max = int(round(rectangle.y_max))

        data[y_min:y_max, x_min:x_max] = 1

        # Return a new Mask object
        return cls(data)

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, x_size, y_size, ellipse):

        """
        This function ...
        :param ellipse:
        :return:
        """

        rel_center = ellipse.center

        a = ellipse.radius.x if isinstance(ellipse.radius, Extent) else ellipse.radius
        b = ellipse.radius.y if isinstance(ellipse.radius, Extent) else ellipse.radius

        # theta in radians !
        theta = ellipse.angle.radian

        x_min = - rel_center.x
        x_max = x_size - rel_center.x
        y_min = - rel_center.y
        y_max = y_size - rel_center.y

        # Calculate the mask
        fraction = elliptical_overlap_grid(x_min, x_max, y_min, y_max, x_size, y_size, a, b, theta, use_exact=0, subpixels=1)

        #xmin, xmax, ymin, ymax : float
        #    Extent of the grid in the x and y direction.
        #nx, ny : int
        #    Grid dimensions.
        #rx : float
        #    The semimajor axis of the ellipse.
        #ry : float
        #    The semiminor axis of the ellipse.
        #theta : float
        #    The position angle of the semimajor axis in radians (counterclockwise).
        #use_exact : 0 or 1
        #    If set to 1, calculates the exact overlap, while if set to 0, uses a
        #    subpixel sampling method with ``subpixel`` subpixels in each direction.
        #subpixels : int
        #    If ``use_exact`` is 0, each pixel is resampled by this factor in each
        #    dimension. Thus, each pixel is divided into ``subpixels ** 2``
        #    subpixels.

        return cls(fraction)

    # -----------------------------------------------------------------

    def fill_holes(self):

        """
        This function ...
        :return:
        """

        # Create a copy of this mask
        new_mask = self.copy()

        # Perform the segmentation
        segments = detect_sources(self.inverse().astype(float), 0.5, 1).data

        # Find the label of the largest segment (=the background)
        label_counts = np.bincount(segments.flatten())
        background_label = np.argmax(label_counts[1:]) + 1
        # If the source mask is larger than the background (in number of pixels), the above will provide the correct label
        # therefore we do the '[1:]'

        # Create a mask for the holes identified as background
        holes = self.inverse().copy()
        holes[segments == background_label] = False

        # Remove holes from the mask
        new_mask[holes] = True

        # Return the new mask
        return new_mask

    # -----------------------------------------------------------------

    def fill_small_holes(self, radius):

        """
        This function ...
        :param radius:
        :return:
        """

        # This function will come in a later version of sckits-image
        #output_mask = morphology.remove_small_holes(self, radius, connectivity=1)

        #return output_mask

        pass

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    def eroded(self, structure=None, connectivity=2, iterations=100):

        """
        This function ...
        :param connectivity:
        :param iterations:
        :return:
        """

        # Define the structure for the expansion
        if structure is None: structure = ndimage.generate_binary_structure(2, connectivity=connectivity)

        try:
            # Make the new mask, made from 100 iterations with the structure array
            data = ndimage.binary_erosion(self, structure, iterations)
        except:
            print(self)
            print(structure)

        # Reassign this object
        return Mask(data, self.selected, self.description)

    # -----------------------------------------------------------------

    def opening(self, structure, iterations=1):

        """
        This function ...
        :return:
        """

        data = ndimage.binary_opening(self, structure, iterations)

        # Return the new mask
        return Mask(data)

    # -----------------------------------------------------------------

    def closed(self, structure, iterations=1):

        """
        This function ...
        :param structure:
        :return:
        """

        data = ndimage.binary_closing(self, structure, iterations)

        # Return the new mask
        return Mask(data)

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    def remove_appendages(self, super=False):

        """
        This function ...
        :return:
        """

        if super: structure = morphology.disk(5, dtype=bool)
        else:
            structure = np.array([[False, True, True, True, False],
                                  [True, True, True, True, True],
                                  [True, True, True, True, True],
                                  [True, True, True, True, True],
                                  [False, True, True, True, False]])

        mask = self.opening(structure)

        segments = detect_sources(mask, 0.5, 1).data

        # Get the label of the center segment
        label = segments[int(0.5*segments.shape[0]), int(0.5*segments.shape[1])]

        # Return the new mask with the appendages removed
        return Mask((segments == label))

    # -----------------------------------------------------------------

    def apply(self, frame, fill=0.0):

        """
        This function ...
        """

        # Replace the masked pixel values with the given fill value
        frame[self] = fill

    # -----------------------------------------------------------------

    def inverse(self):

        """
        This function ...
        """

        # Return the inverse of this mask
        return np.logical_not(self)

    # -----------------------------------------------------------------

    def intersection(self, mask):

        """
        This function ...
        :param mask_a:
        :param mask_b:
        :return:
        """

        assert isinstance(self * mask, Mask)

        return self * mask

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    def masks(self, position):

        """
        This function ...
        """

        # Calculate x and y of the pixel corresponding to the object's position
        x_pixel = int(round(position.x))
        y_pixel = int(round(position.y))

        # Check whether this box contains the position
        if x_pixel < 0 or y_pixel < 0 or x_pixel >= self.xsize or y_pixel >= self.ysize: return False
        else:

            # Return the value of the mask in this pixel
            return self[y_pixel, x_pixel]

# -----------------------------------------------------------------
