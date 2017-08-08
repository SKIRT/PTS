#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.segmentationmap Contains the SegmentationMap class, representing a 2D frame with
#  integer-labeled segments.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from photutils.segmentation import SegmentationImage
from astropy.io import fits

# Import the relevant PTS classes and modules
from ...core.basics.log import log

# -----------------------------------------------------------------

# SEGMENTATIONIMAGE class: (as of version 0.3.dev1880)
# source: http://photutils.readthedocs.io/en/latest/api/photutils.segmentation.SegmentationImage.html#photutils.segmentation.SegmentationImage

# Attributes Summary

# areas	The areas (in pixel**2) of all labeled regions.
# array	The 2D segmentation image.
# data	The 2D segmentation image.
# data_masked	A MaskedArray version of the segmentation image where the background (label = 0) has been masked.
# is_sequential	Determine whether or not the non-zero labels in the segmenation image are sequential (with no missing values).
# labels	The sorted non-zero labels in the segmentation image.
# max	The maximum non-zero label in the segmentation image.
# nlabels	The number of non-zero labels in the segmentation image.
# shape	The shape of the 2D segmentation image.
# slices	The minimal bounding box slices for each labeled region.

# Methods Summary

# area(labels)	The areas (in pixel**2) of the regions for the input labels.
# check_label(label[, allow_zero])	Check for a valid label label number within the segmentation image.
# copy()	Return a deep copy of this class instance.
# keep_labels(labels[, relabel])	Keep only the specified label numbers.
# outline_segments([mask_background])	Outline the labeled segments.
# relabel(labels, new_label)	Relabel one or more label numbers.
# relabel_sequential([start_label])	Relabel the label numbers sequentially, such that there are no missing label numbers (up to the maximum label number).
# remove_border_labels(border_width[, ...])	Remove labeled segments near the image border.
# remove_labels(labels[, relabel])	Remove one or more label numbers.
# remove_masked_labels(mask[, ...])	Remove labeled segments located within a masked region.

# -----------------------------------------------------------------

class SegmentationMap(SegmentationImage):

    """
    This class ...
    """

    def __init__(self, data, **kwargs):

        """
        The constructor ...
        :param data:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SegmentationMap, self).__init__(data)

        # Set the WCS
        self.wcs = kwargs.pop("wcs", None)

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, index=None, plane=None, hdulist_index=None):

        """
        This function ...
        :param path:
        :param index:
        :param plane:
        :param hdulist_index:
        :return:
        """

        name = None
        description = None
        no_filter = True
        fwhm = None
        add_meta = False

        from .fits import load_frame # Import here because .fits imports SegmentationMap

        try:
            # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
            segments = load_frame(cls, path, index, name, description, plane, hdulist_index, no_filter, fwhm, add_meta=add_meta)
        except TypeError: raise IOError("The file is possibly damaged")

        # Set the path
        segments.path = path

        # Return the segmentation map
        return segments

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

    def add_shape(self, shape):

        """
        This function ...
        :return:
        """

        # Create mask
        mask = shape.to_mask(self.xsize, self.ysize)

        # Add the mask as integer type
        self._data += mask.astype(int)

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the segmentation map ...")

        # Check whether the path is valid
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path, header=None):

        """
        This function ...
        :param path:
        :param header:
        :return:
        """

        # If a header is not specified, created it from the WCS
        if header is None: header = self.header

        # FITS format
        if path.endswith(".fits"):

            from .fits import write_frame  # Import here because io imports SegmentationMap

            # Write to a FITS file
            write_frame(self._data, header, path)

        # ASDF format
        elif path.endswith(".asdf"):

            # Import
            from asdf import AsdfFile

            # Create the tree
            tree = dict()

            tree["data"] = self._data
            tree["header"] = header

            # Create the asdf file
            ff = AsdfFile(tree)

            # Write
            ff.write_to(path)

        # Invalid
        else: raise ValueError("Only the FITS or ASDF filetypes are supported")

        # Update the path
        self.path = path

# -----------------------------------------------------------------
