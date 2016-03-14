#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
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
from astropy.io import fits
from photutils import detect_sources

# Import the relevant PTS classes and modules
from .vector import Position
from ...core.tools.logging import log

# -----------------------------------------------------------------

class Mask(np.ndarray):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __new__(cls, data, name=None, description=None):

        """
        This function ...
        :param cls:
        :param data:
        :param description:
        :return:
        """

        obj = np.asarray(data, dtype=bool).view(cls)
        obj.name = name
        obj.description = description

        return obj

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, index=0):

        """
        This function ...
        :param path:
        :param index:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Open the HDU list for the FITS file
        hdulist = fits.open(path)

        # Get the primary HDU
        hdu = hdulist[0]

        # Get the image header
        header = hdu.header

        # Check whether multiple planes are present in the FITS image
        nframes = 1
        if 'NAXIS' in header:
            # If there are 3 axes, get the size of the third
            if header['NAXIS'] == 3: nframes = header['NAXIS3']

        if nframes > 1:

            mask = cls(hdu.data[index])

        else:

            # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
            if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

            mask = cls(hdu.data)

        # Return the new mask
        return mask

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

    @classmethod
    def is_zero(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(frame == 0.0)

    # -----------------------------------------------------------------

    @classmethod
    def is_zero_or_less(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(frame <= 0.0)

    # -----------------------------------------------------------------

    @classmethod
    def empty(cls, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        return cls(np.zeros((y_size, x_size)))

    # -----------------------------------------------------------------

    @classmethod
    def empty_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(np.zeros(frame.shape))

    # -----------------------------------------------------------------

    @classmethod
    def full(cls, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        return cls(np.ones((y_size, x_size)))

    # -----------------------------------------------------------------

    @classmethod
    def full_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(np.ones(frame.shape))

    # -----------------------------------------------------------------

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        if obj is None: return
        self.name = getattr(obj, 'name', None)
        self.description = getattr(obj, 'description', None)

    # -----------------------------------------------------------------

    @property
    def xsize(self): return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self): return self.shape[0]

    # -----------------------------------------------------------------

    @classmethod
    def from_shape(cls, shape, x_size, y_size, invert=False):

        """
        This function ...
        :param shape:
        :param x_size:
        :param y_size:
        :param invert:
        :return:
        """

        # Return a new Mask object
        mask = shape.to_mask(x_size, y_size)

        # Return the mask (inverted if requested)
        if invert: return mask.inverse()
        else: return mask

    # -----------------------------------------------------------------

    @classmethod
    def from_region(cls, region, x_size, y_size):

        """
        This function ...
        :param region:
        :param x_size:
        :param y_size
        :return:
        """

        # Return a new Mask object
        return region.to_mask(x_size, y_size)

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
        if len(label_counts) > 1:

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

        # Return the dilated mask
        #data, name=None, description=None
        return Mask(data, name=self.name, description=self.description)

    # -----------------------------------------------------------------

    def disk_dilation(self, radius=5, iterations=1):

        """
        This function ...
        :param radius:
        :return:
        """

        structure = morphology.disk(radius, dtype=bool)
        data = ndimage.binary_dilation(self, structure, iterations)

        # Return the dilated mask
        #data, name=None, description=None
        return Mask(data, name=self.name, description=self.description)

    # -----------------------------------------------------------------

    def eroded(self, structure=None, connectivity=2, iterations=100):

        """
        This function ...
        :param structure:
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
            #print(self)
            #print(structure)
            data = np.zeros((self.ysize,self.xsize), dtype=bool)

        # Reassign this object
        #data, name=None, description=None
        return Mask(data, name=self.name, description=self.description)

    # -----------------------------------------------------------------

    def opening(self, structure, iterations=1):

        """
        This function ...
        :param structure:
        :param iterations:
        :return:
        """

        data = ndimage.binary_opening(self, structure, iterations)

        # Return the new mask
        #data, name=None, description=None
        return Mask(data, name=self.name, description=self.description)

    # -----------------------------------------------------------------

    def closed(self, structure, iterations=1):

        """
        This function ...
        :param structure:
        :return:
        """

        data = ndimage.binary_closing(self, structure, iterations)

        # Return the new mask
        #data, name=None, description=None
        return Mask(data, name=self.name, description=self.description)

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
        #data, name=None, description=None
        return Mask((segments == label), name=self.name, description=self.description)

    # -----------------------------------------------------------------

    def apply(self, frame, fill_value=0.0):

        """
        This function ...
        """

        # Replace the masked pixel values with the given fill value
        frame[self] = fill_value

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
        :param mask:
        :return:
        """

        assert isinstance(self * mask, Mask)

        return self * mask

    # -----------------------------------------------------------------

    def subtract(self, mask):

        """
        This function ...
        :param mask:
        :return:
        """

        assert isinstance(mask.inverse(), Mask)

        # Return ...
        return self.intersection(mask.inverse())

    # -----------------------------------------------------------------

    def union(self, mask):

        """
        This function ...
        :param mask:
        :return:
        """

        assert isinstance(self + mask, Mask)

        # Return ...
        return self + mask

    # -----------------------------------------------------------------

    def hits_boundary(self, min_pixels=2, where=False):

        """
        This function ...
        """

        # Set the number of hits (pixels hitting the boundary) to zero initially
        hits = 0

        positions = []

        for x in range(self.xsize):
            for y in (0, self.ysize-1):

                # Check if pixel is masked
                if self[y, x]:

                    # Add one hit
                    hits += 1

                    if where: positions.append(Position(x, y))

                    # Break the loop, avoid performing unnecessary checks (if we don't need to output where the hits were)
                    if hits >= min_pixels and not where: break

        for y in range(1, self.ysize-1):
            for x in (0, self.xsize-1):

                # Check if pixel is masked
                if self[y, x]:

                    # Add one hit
                    hits += 1

                    if where: positions.append(Position(x, y))

                    # Break the loop, avoid performing unnecessary checks (if we don't need to output where the hits were)
                    if hits >= min_pixels and not where: break

        # Return whether ...
        if where: return hits >= min_pixels, positions
        else: return hits >= min_pixels

    # -----------------------------------------------------------------

    def masks(self, position):

        """
        This function ...
        :param position:
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

    def covers(self, mask):

        """
        This function ...
        :param mask:
        :return:
        """

        not_covered = mask & self.inverse()
        return not np.any(not_covered)

# -----------------------------------------------------------------
