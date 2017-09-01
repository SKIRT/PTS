#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.alpha Contains the AlphaMask class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import numpy as np

# Import astronomical modules
from astropy.io import fits
from photutils import detect_sources
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ..basics.vector import PixelShape
from ..dist_ellipse import distance_ellipse

# -----------------------------------------------------------------

class AlphaMask(object):
    
    """
    This class ...
    """

    def __init__(self, data, **kwargs):

        """
        The constructor ...
        :param data:
        :param kwargs:
        """

        # Set data
        self._data = data.astype(np.uint8)

        # Set the WCS
        self.wcs = kwargs.pop("wcs", None)

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, ellipse, shape, factor_range, wcs=None):

        """
        Thisf unction ...
        :param ellipse:
        :param shape:
        :param factor_range:
        :param wcs:
        :return:
        """

        center = ellipse.center

        # angle = - region.angle + Angle(-90., "deg")
        angle = ellipse.angle + Angle(90., "deg")

        # Determine the ratio of semimajor and semiminor
        ratio = ellipse.semiminor / ellipse.semimajor
        radius = distance_ellipse(shape, center, ratio, angle) / ellipse.semiminor

        outside_max = radius > factor_range.max
        inside_min = radius < factor_range.min

        test = (factor_range.max - radius) / factor_range.span

        alpha_channel = test
        alpha_channel[inside_min] = 1
        alpha_channel[outside_max] = 0

        # Create alpha mask
        alpha = cls.from_real(alpha_channel, wcs=wcs)

        # Return
        return alpha

    # -----------------------------------------------------------------

    @classmethod
    def from_real(cls, data, **kwargs):

        """
        This function ...
        :param data:
        :param kwargs:
        :return:
        """

        return cls(data * 255, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def between(cls, data, value_range):

        """
        This function ...
        :param data:
        :param value_range:
        :return:
        """

        # Create alpha
        #alpha = np.zeros_like(data)

        if hasattr(data, "data"): data = data.data

        # Create alpha
        test = (data - value_range.min) / value_range.span
        alpha = test

        # Set above max and below min
        above_max = data > value_range.max
        below_min = data < value_range.min
        alpha[above_max] = 1
        alpha[below_min] = 0

        # Create the alpha mask
        return cls.from_real(alpha)

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return PixelShape.from_tuple(self.data.shape)

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.shape[0]

    # -----------------------------------------------------------------

    @property
    def data(self):

        """
        This function ...
        :return:
        """

        return self._data

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def as_real(self):

        """
        This function ...
        :return:
        """

        return self.data.astype(float) / 255

    # -----------------------------------------------------------------

    def invert(self):

        """
        This function ...
        :return:
        """

        self._data = 255 - self._data

    # -----------------------------------------------------------------

    def inverse(self):

        """
        This function ...
        :return:
        """

        # Create new, inverted copy
        new = self.copy()
        new.invert()
        return new

    # -----------------------------------------------------------------

    @property
    def header(self):

        """
        This function ...
        :return:
        """

        # If the WCS for this frame is defined, use it to create a header
        if self.wcs is not None: header = self.wcs.to_header()

        # Else, create a new empty header
        else: header = fits.Header()

        # Add properties to the header
        header['NAXIS'] = 2
        header['NAXIS1'] = self.xsize
        header['NAXIS2'] = self.ysize

        # Return the header
        return header

    # -----------------------------------------------------------------

    def central(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels: number of connected pixels
        :param connectivity:
        :return:
        """

        # Perform the segmentation
        segments = detect_sources(self.data, 1, npixels=npixels, connectivity=connectivity).data

        # To plot the multiple segments that are detected
        # if segments.max() > 1: plotting.plot_box(np.ma.masked_array(box, mask=segments.astype(bool)))

        # Center
        center_x = 0.5 * (self.xsize + 1) - 1
        center_y = 0.5 * (self.ysize + 1) - 1

        # Get the label of the center segment
        #rel_center = self.cutout.rel_position(self.center)
        label = segments[int(round(center_y)), int(round(center_x))]

        # If the center pixel is identified as being part of the background, create an empty mask (the center does not
        # correspond to a segment)
        if label == 0: return self.__class__(np.zeros_like(self.data))

        # Return copy with only largest
        else:
            new = self.copy()
            new._data[segments != label] = 0
            return new

    # -----------------------------------------------------------------

    def largest(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels:
        :param connectivity:
        :return:
        """

        # Perform the segmentation
        segments = detect_sources(self.data, 1, npixels=npixels, connectivity=connectivity).data

        # Get counts for each label
        unique, counts = np.unique(segments, return_counts=True)

        # Get indices of the unique values (sorted on
        sorted_indices = np.argsort(counts)

        # Check last index
        last_index = sorted_indices[-1]
        value = unique[last_index]

        # Get the label
        if value == 0:

            # No other labels: no patches
            if len(sorted_indices) == 1: return self.__class__(np.zeros_like(self.data))

            second_last_index = sorted_indices[-2]
            label = unique[second_last_index]

        else: label = value

        # Return copy with only largest
        new = self.copy()
        #print(segments != label)
        new._data[segments != label] = 0
        return new

    # -----------------------------------------------------------------

    def fill_holes(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels:
        :param connectivity:
        :return:
        """

        # Perform the segmentation
        segments = detect_sources(self.inverse().data, 1, npixels=npixels, connectivity=connectivity).data

        # Find the label of the largest segment (=the background)
        label_counts = np.bincount(segments.flatten())
        if len(label_counts) > 1:

            background_label = np.argmax(label_counts[1:]) + 1
            # If the source mask is larger than the background (in number of pixels), the above will provide the correct label
            # therefore we do the '[1:]'

            # Create a mask for the holes identified as background
            holes = self.inverse().data > 1
            holes[segments == background_label] = False

            # Remove holes from the mask
            self._data[holes] = 255

    # -----------------------------------------------------------------

    def to_rgb(self, colour="black", background_color="white"):

        """
        This function ...
        :param colour:
        :param background_color:
        :return:
        """

        from .rgb import RGBImage
        return RGBImage.from_alpha_mask(self, colour=colour, background_color=background_color)

    # -----------------------------------------------------------------

    def to_rgba(self, colour="black"):

        """
        This function ...
        :param colour:
        :return:
        """

        from .rgba import RGBAImage
        return RGBAImage.from_alpha_mask(self, colour=colour)

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the alpha mask ...")

        # Check whether the path is valid
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path, header=None, update_path=True, colour="black", background_color="white"):

        """
        This function ...
        :param path:
        :param header:
        :param update_path:
        :param colour:
        :param background_color:
        :return:
        """

        # FITS format
        if path.endswith(".fits"): self.saveto_fits(path, header=header, update_path=update_path)

        # PNG format
        elif path.endswith(".png"): self.saveto_png(path, colour=colour, background_color=background_color)

        # Invalid
        else: raise ValueError("Invalid file format")

    # -----------------------------------------------------------------

    def saveto_fits(self, path, header=None, update_path=True):

        """
        This function ...
        :param path:
        :param header:
        :param update_path:
        :return:
        """

        # If a header is not specified, created it from the WCS
        if header is None: header = self.header

        from .fits import write_frame  # Import here because io imports Mask

        # Write to a FITS file
        write_frame(self.data, header, path)

        # Update the path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    def saveto_png(self, path, colour="black", background_color="white", alpha=False):

        """
        This function ...
        :param path:
        :param colour:
        :param background_color:
        :param alpha:
        :return:
        """

        # Get RGB image
        if alpha: image = self.to_rgba(colour=colour)
        else: image = self.to_rgb(colour=colour, background_color=background_color)

        # Save RGB image
        image.saveto(path)

# -----------------------------------------------------------------

def product(*args):

    """
    This function ...
    :param args:
    :return:
    """

    # INTERSECTION = 1 * first * second * ... (1 is neutral element for multiplication)
    # so for one mask, intersection = 1 * mask = mask

    if len(args) == 1: return AlphaMask(args[0].data)

    arrays = [arg.as_real() for arg in args]
    # arrays = []
    # for arg in args:
    #     if isinstance(arg, MaskBase): arrays.append(arg.data)
    #     elif isinstance(arg, oldMask): arrays.append(arg)
    #     else: arrays.append(arg)
    return AlphaMask.from_real(np.product(arrays, axis=0))

# -----------------------------------------------------------------
