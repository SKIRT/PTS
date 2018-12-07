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

min_alpha = 0
max_alpha = 255
half_alpha = 128
nalpha_bins = 256

# -----------------------------------------------------------------

def load_mask_or_alpha_mask(path, index=None, plane=None, hdulist_index=None, no_wcs=False):

    """
    This function ...
    :param path:
    :param index:
    :param plane:
    :param hdulist_index:
    :param no_wcs:
    :return:
    """

    name = None
    description = None
    no_filter = True
    fwhm = None
    add_meta = False

    from . import fits as pts_fits

    def is_binary(data):
        return np.array_equal(data, data.astype(bool))

    from .mask import Mask

    def mask_class_picker(data):
        if is_binary(data): return Mask
        else: return AlphaMask

    # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
    mask = pts_fits.load_frame("mask", path, index, name, description, plane, hdulist_index, no_filter, fwhm,
                               add_meta=add_meta, no_wcs=no_wcs, class_picker=mask_class_picker)

    # Set the path
    mask.path = path

    # Return the mask or alpha mask
    return mask

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
        self._wcs = kwargs.pop("wcs", None)

        # Set the pixelscale
        self._pixelscale = kwargs.pop("pixelscale", None)

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, index=None, plane=None, hdulist_index=None, no_wcs=False):

        """
        This function ...
        :param path:
        :param index:
        :param plane:
        :param hdulist_index:
        :param no_wcs:
        :return:
        """

        name = None
        description = None
        no_filter = True
        fwhm = None
        add_meta = False

        from . import fits as pts_fits

        # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
        mask = pts_fits.load_frame(cls, path, index, name, description, plane, hdulist_index, no_filter, fwhm, add_meta=add_meta, no_wcs=no_wcs)

        # Set the path
        mask.path = path

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, ellipse, shape, factor_range, wcs=None, dynamic_max=False):

        """
        Thisf unction ...
        :param ellipse:
        :param shape:
        :param factor_range:
        :param wcs:
        :param dynamic_max:
        :return:
        """

        center = ellipse.center

        # angle = - region.angle + Angle(-90., "deg")
        angle = ellipse.angle + Angle(90., "deg")

        # Determine the ratio of semimajor and semiminor
        ratio = ellipse.semiminor / ellipse.semimajor
        radius = distance_ellipse(shape, center, ratio, angle) / ellipse.semiminor

        #from ..tools import plotting
        #print("factor:", factor_range)
        #plotting.plot_box(radius)

        max_radius = np.max(radius)
        if dynamic_max and factor_range.max > np.max(radius): max_factor = max_radius
        else: max_factor = factor_range.max
        min_factor = factor_range.min

        outside_max = radius > max_factor
        inside_min = radius < min_factor

        #plotting.plot_mask(outside_max, title="outside maximum factor")
        #plotting.plot_mask(inside_min, title="inside minimum factor")

        test = (factor_range.max - radius) / factor_range.span

        #plotting.plot_box(test, title="test")

        alpha_channel = test
        alpha_channel[inside_min] = 1
        alpha_channel[outside_max] = 0

        #plotting.plot_box(alpha_channel, title="final")

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

        return cls(data * max_alpha, **kwargs)

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

    @classmethod
    def full(cls, shape, wcs=None, pixelscale=None, value=255):

        """
        This function ...
        :param shape:
        :param wcs:
        :param pixelscale:
        :param value:
        :return:
        """

        # Create data
        data = np.full(shape, value, dtype=np.uint8)

        # Create alpha mask
        return cls(data, wcs=wcs, pixelscale=pixelscale)

    # -----------------------------------------------------------------

    @classmethod
    def full_like(cls, other, wcs=None, pixelscale=None, value=255):

        """
        This function ...
        :param other:
        :param wcs:
        :param pixelscale:
        :param value:
        :return:
        """

        # Set WCS and/or pixelscale
        if wcs is None and hasattr(other, 'wcs'): wcs = other.wcs
        if pixelscale is None and hasattr(other, 'pixelscale'): pixelscale = other.pixelscale

        # Create and return
        return cls.full(other.shape, wcs=wcs, pixelscale=pixelscale, value=value)

    # -----------------------------------------------------------------

    @classmethod
    def empty(cls, shape, wcs=None, pixelscale=None):

        """
        This function ...
        :param shape:
        :param wcs:
        :param pixelscale:
        :return:
        """

        # Create data
        data = np.zeros(shape, dtype=np.uint8)

        # Create alpha mask
        return cls(data, wcs=wcs, pixelscale=pixelscale)

    # -----------------------------------------------------------------

    @classmethod
    def empty_like(cls, other, wcs=None, pixelscale=None):

        """
        This function ..
        :param other:
        :param wcs:
        :param pixelscale:
        :return:
        """

        # Set WCS and/or pixelscale
        if wcs is None and hasattr(other, 'wcs'): wcs = other.wcs
        if pixelscale is None and hasattr(other, 'pixelscale'): pixelscale = other.pixelscale

        # Create and return
        return cls.empty(other.shape, wcs=wcs, pixelscale=pixelscale)

    # -----------------------------------------------------------------

    @classmethod
    def from_levels(cls, levels, wcs=None, pixelscale=None):

        """
        This function ...
        :param levels:
        :param wcs:
        :param pixelscale:
        :return:
        """

        if len(levels) == 0: raise ValueError("No level data is given")

        # Initialize
        first = levels[levels.keys()[0]]
        alpha = cls.empty_like(first, wcs=wcs, pixelscale=pixelscale)
        #alpha = cls.full_like(first, wcs=wcs, pixelscale=pixelscale, value=min_level)

        # Set levels
        for level in sorted(levels.keys()):

            mask = levels[level]
            alpha._data[mask] = level

        # Return the alpha mask
        return alpha

    # -----------------------------------------------------------------

    @property
    def shape(self):
        return PixelShape.from_tuple(self.data.shape)

    # -----------------------------------------------------------------

    @property
    def xsize(self):
        return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self):
        return self.shape[0]

    # -----------------------------------------------------------------

    @property
    def x_center(self):
        return 0.5 * (self.xsize - 1)

    # -----------------------------------------------------------------

    @property
    def y_center(self):
        return 0.5 * (self.ysize - 1)

    # -----------------------------------------------------------------

    @property
    def center(self):
        from ..basics.coordinate import PixelCoordinate
        return PixelCoordinate(self.x_center, self.y_center)

    # -----------------------------------------------------------------

    @property
    def data(self):
        return self._data

    # -----------------------------------------------------------------

    @property
    def wcs(self):
        return self._wcs

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, value):
        self._wcs = value

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):
        return self.wcs is not None

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):
        # Return the pixelscale of the WCS is WCS is defined
        if self.wcs is not None: return self.wcs.pixelscale
        else: return self._pixelscale  # return the pixelscale

    # -----------------------------------------------------------------

    @pixelscale.setter
    def pixelscale(self, value):
        from ..basics.pixelscale import Pixelscale
        if not isinstance(value, Pixelscale): value = Pixelscale(value)
        self._pixelscale = value

    # -----------------------------------------------------------------

    @property
    def has_pixelscale(self):
        return self.pixelscale is not None

    # -----------------------------------------------------------------

    @property
    def x_pixelscale(self):
        return self.pixelscale.x

    # -----------------------------------------------------------------

    @property
    def y_pixelscale(self):
        return self.pixelscale.y

    # -----------------------------------------------------------------

    @property
    def has_angular_pixelscale(self):
        from ...core.tools import types
        if not self.has_pixelscale: raise ValueError("No pixelscale")
        return types.is_angle(self.x_pixelscale)

    # -----------------------------------------------------------------

    @property
    def has_physical_pixelscale(self):
        from ...core.tools import types
        if not self.has_pixelscale: raise ValueError("No pixelscale")
        return types.is_length_quantity(self.x_pixelscale)

    # -----------------------------------------------------------------

    @property
    def average_pixelscale(self):
        if self.wcs is not None: return self.wcs.average_pixelscale
        else: return self._pixelscale.average if self._pixelscale is not None else None

    # -----------------------------------------------------------------

    def copy(self):
        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def as_real(self):
        return self.data.astype(float) / max_alpha

    # -----------------------------------------------------------------

    def invert(self):

        """
        This function ...
        :return:
        """

        self._data = max_alpha - self._data

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

    def mask(self, where):

        """
        This function ...
        :param where:
        :return:
        """

        self.data[where] = max_alpha

    # -----------------------------------------------------------------

    def unmask(self, where):

        """
        This function ...
        :param where:
        :return:
        """

        self.data[where] = min_alpha

    # -----------------------------------------------------------------

    def central(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels: number of connected pixels
        :param connectivity:
        :return:
        """

        # Perform the segmentation
        segments = detect_sources(self.data, min_alpha, npixels=npixels, connectivity=connectivity).data

        # To plot the multiple segments that are detected
        # if segments.max() > 1: plotting.plot_box(np.ma.masked_array(box, mask=segments.astype(bool)))

        # Get the label of the center segment
        #rel_center = self.cutout.rel_position(self.center)
        label = segments[int(round(self.y_center)), int(round(self.x_center))]

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
        segments = detect_sources(self.data, min_alpha, npixels=npixels, connectivity=connectivity).data

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
        new._data[segments != label] = 0
        return new

    # -----------------------------------------------------------------

    @property
    def min(self):
        return np.min(self.data)

    # -----------------------------------------------------------------

    @property
    def max(self):
        return np.max(self.data)

    # -----------------------------------------------------------------

    @property
    def transparent(self):
        from .mask import Mask
        return Mask(self.data == min_alpha, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    @property
    def transparent_filled(self):
        transparent = self.transparent
        transparent.fill_holes(connectivity=4)
        return transparent

    # -----------------------------------------------------------------

    @property
    def semitransparent(self):
        from .mask import Mask
        #return Mask(min_alpha < self.data < max_alpha, wcs=self.wcs, pixelscale=self.pixelscale)
        above = self.data > min_alpha
        below = self.data < max_alpha
        return Mask(above * below, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    @property
    def semitransparent_filled(self):

        """
        This function ...
        :return:
        """

        semitransparent = self.semitransparent
        semitransparent.fill_holes(connectivity=4)
        return semitransparent

    # -----------------------------------------------------------------

    @property
    def opaque(self):
        from .mask import Mask
        return Mask(self.data == max_alpha, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    @property
    def opaque_filled(self):
        opaque = self.opaque
        opaque.fill_holes(connectivity=4)
        return opaque

    # -----------------------------------------------------------------

    def equal(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from .mask import Mask
        return Mask(self.data == value, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def equal_filled(self, value):

        """
        Thisf unction ...
        :param value:
        :return:
        """

        equal = self.equal(value)
        equal.fill_holes(connectivity=4)
        return equal

    # -----------------------------------------------------------------

    def above(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from .mask import Mask
        return Mask(self.data > value, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def above_filled(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        above = self.equal(value)
        above.fill_holes(connectivity=4)
        return above

    # -----------------------------------------------------------------

    def above_or_equal(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from .mask import Mask
        return Mask(self.data >= value, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def above_or_equal_filled(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        above_or_equal = self.above_or_equal(value)
        above_or_equal.fill_holes(connectivity=4)
        return above_or_equal

    # -----------------------------------------------------------------

    def below(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from .mask import Mask
        return Mask(self.data < value, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def below_filled(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        below = self.below(value)
        below.fill_holes(connectivity=4)
        return below

    # -----------------------------------------------------------------

    def below_or_equal(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from .mask import Mask
        return Mask(self.data <= value, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def below_or_equal_filled(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        below_or_equal = self.below_or_equal(value)
        below_or_equal.fill_holes(connectivity=4)
        return below_or_equal

    # -----------------------------------------------------------------

    @property
    def above_half(self):

        """
        This function ...
        :return:
        """

        return self.above_or_equal(half_alpha)

    # -----------------------------------------------------------------

    @property
    def above_half_filled(self):

        """
        This function ...
        :return:
        """

        return self.above_or_equal_filled(half_alpha)

    # -----------------------------------------------------------------

    @property
    def below_half(self):

        """
        This function ...
        :return:
        """

        return self.below(half_alpha)

    # -----------------------------------------------------------------

    @property
    def below_half_filled(self):

        """
        This function ...
        :return:
        """

        return self.below_filled(half_alpha)

    # -----------------------------------------------------------------

    def between_levels(self, min_level, max_level, include_min=False, include_max=False):

        """
        This function ...
        :param min_level:
        :param max_level:
        :param include_min:
        :param include_max:
        :return:
        """

        from .mask import Mask

        # Above (or equal) minimum value
        if include_min: above = self.data >= min_level
        else: above = self.data > min_level

        # Below (or equal) maximum value
        if include_max: below = self.data <= max_level
        else: below = self.data < max_level

        # Create mask
        return Mask(above * below, wcs=self.wcs, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def between_levels_filled(self, min_level, max_level, include_min=False, include_max=False):

        """
        This function ...
        :param min_level:
        :param max_level:
        :param include_min:
        :param include_max:
        :return:
        """

        between = self.between_levels(min_level, max_level, include_min=include_min, include_max=include_max)
        between.fill_holes(connectivity=4)
        return between

    # -----------------------------------------------------------------

    @property
    def between_min_max(self):
        return self.between_levels(self.lowest_level, self.highest_level, include_min=False, include_max=True)

    # -----------------------------------------------------------------

    @property
    def between_min_max_filled(self):
        between = self.between_min_max
        between.fill_holes(connectivity=4)
        return between

    # -----------------------------------------------------------------

    @property
    def levels(self):
        return list(np.unique(self.data))

    # -----------------------------------------------------------------

    @property
    def nlevels(self):
        return len(self.levels)

    # -----------------------------------------------------------------

    @property
    def lowest_level(self):
        return min(self.levels)

    # -----------------------------------------------------------------

    @property
    def highest_level(self):
        return max(self.levels)

    # -----------------------------------------------------------------

    @property
    def at_lowest_level(self):
        return self.equal(self.lowest_level)

    # -----------------------------------------------------------------

    @property
    def at_highest_level(self):
        return self.equal(self.highest_level)

    # -----------------------------------------------------------------

    @property
    def levels_between_min_max(self):

        """
        This function ...
        :return:
        """

        levels = self.levels
        if len(levels) < 3: return []
        else: return self.levels[1:-1]

    # -----------------------------------------------------------------

    @property
    def levels_between_transparent_opaque(self):

        """
        Thisf unction ...
        :return:
        """

        levels = self.levels
        if levels[0] == min_alpha: levels = levels[1:]
        if levels[-1] == max_alpha: levels = levels[:-1]
        return levels

    # -----------------------------------------------------------------

    @property
    def lowest_level_not_transparent(self):
        return self.levels_between_transparent_opaque[0]

    # -----------------------------------------------------------------

    @property
    def highest_level_not_opaque(self):
        return self.levels_between_transparent_opaque[-1]

    # -----------------------------------------------------------------

    @property
    def levels_semitransparent(self):
        return self.levels_between_transparent_opaque

    # -----------------------------------------------------------------

    def get_holes(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels:
        :param connectivity:
        :return:
        """

        # Perform the segmentation
        segments = detect_sources(self.inverse().data, self.max-1, npixels=npixels, connectivity=connectivity).data

        #from ..tools import plotting
        #plotting.plot_mask(self.inverse())
        #plotting.plot_box(segments)

        # Find the label of the largest segment (=the background)
        label_counts = np.bincount(segments.flatten())
        if len(label_counts) <= 1: return

        # Determine background label
        background_label = np.argmax(label_counts[1:]) + 1
        # If the source mask is larger than the background (in number of pixels), the above will provide the correct label
        # therefore we do the '[1:]'
        #print(background_label)

        # Create a mask for the holes identified as background
        holes = self.inverse().data > 1
        holes[segments == background_label] = False

        # Get 'edge', unmark as holes
        #core = self.opaque_filled
        core = self.above_half_filled
        edge = core.inverse()
        holes[edge] = False

        #plotting.plot_mask(holes, title="holes")
        #plotting.plot_mask(self.opaque_filled, title="opaque filled")

        # Return
        return holes

    # -----------------------------------------------------------------

    def fill_holes(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels:
        :param connectivity:
        :return:
        """

        # Get the holes
        holes = self.get_holes(npixels=npixels, connectivity=connectivity)

        # Remove holes from the mask
        self._data[holes] = max_alpha

    # -----------------------------------------------------------------

    def filled_holes(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels:
        :param connectivity:
        :return:
        """

        # Create new, inverted copy
        new = self.copy()
        new.fill_holes(npixels=npixels, connectivity=connectivity)
        return new

    # -----------------------------------------------------------------

    def disk_dilate(self, radius=5, niterations=1, nbins=10, max_radius=None):

        """
        This function ...
        :param radius:disk_dilate
        :param niterations:
        :param nbins:
        :param max_radius:
        :return:
        """

        from ..tools import statistics

        # Get all levels
        levels = self.levels_between_min_max

        # Get the minimum and maximum level
        minimum_level = self.lowest_level_not_transparent
        maximum_level = self.highest_level
        level_span = maximum_level - minimum_level

        # Get the bin properties
        centers, lower, upper = statistics.histogram(levels, nbins)

        # Create new data, filled with lowest level
        data = np.full(self.data.shape, self.lowest_level)

        radii = []

        # Loop over the bins, starting with the lower levels
        first = True
        for lower_value, center_value, upper_value in zip(lower, centers, upper):

            # Get integer center level
            center = int(round(center_value))

            # Determine the radius for this level
            if max_radius is not None:
                radius_span = max_radius - radius
                level_radius = radius - radius_span / level_span * (center - maximum_level)
            else: level_radius = radius
            int_level_radius = int(round(level_radius))

            radii.append(int_level_radius)

            #print(center, level_radius)

            # Get mask
            mask = self.between_levels(lower_value, upper_value, include_min=first, include_max=True)
            first = False

            # Dilate the mask
            mask.disk_dilate(radius=int_level_radius, niterations=niterations)

            # Set new data
            data[mask] = center

        # Dilate mask of the maximum value
        highest_mask = self.at_highest_level
        highest_mask.disk_dilate(radius=radius, niterations=niterations)

        # Set
        data[highest_mask] = self.highest_level

        # Set the new data
        self._data = data

        # Return the bins and the radii used for dilation in each bin
        return centers, lower, upper, radii

    # -----------------------------------------------------------------

    def disk_dilated(self, radius=5, niterations=1, nbins=10):

        """
        This function ...
        :param radius:
        :param niterations:
        :param nbins:
        :return:
        """

        new = self.copy()
        new.disk_dilate(radius=radius, niterations=niterations, nbins=nbins)
        return new

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
