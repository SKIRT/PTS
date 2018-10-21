#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.mask Contains the Mask class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.io import fits
from reproject import reproject_exact, reproject_interp

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ..basics.mask import MaskBase
from ..basics.mask import Mask as oldMask
from ..tools import cropping
from ..basics.pixelscale import Pixelscale
from ...core.tools import types

# -----------------------------------------------------------------

class Mask(MaskBase):
    
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
        super(Mask, self).__init__(data, **kwargs)

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

        from . import fits as pts_fits  # Import here because io imports SegmentationMap

        # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
        mask = pts_fits.load_frame(cls, path, index, name, description, plane, hdulist_index, no_filter, fwhm, add_meta=add_meta, no_wcs=no_wcs)

        # Set the path
        mask.path = path

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @classmethod
    def nans_from_file(cls, path, index=None, plane=None, hdulist_index=None, no_wcs=False):

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

        from . import fits as pts_fits  # Import here because io imports SegmentationMap

        data_converter = lambda data: np.isnan(data)

        # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
        mask = pts_fits.load_frame(cls, path, index, name, description, plane, hdulist_index, no_filter, fwhm,
                                   add_meta=add_meta, no_wcs=no_wcs, data_converter=data_converter)

        # Set the path
        mask.path = path

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @classmethod
    def zeroes_from_file(cls, path, index=None, plane=None, hdulist_index=None, no_wcs=False):

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

        from . import fits as pts_fits  # Import here because io imports SegmentationMap

        data_converter = lambda data: data == 0

        # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
        mask = pts_fits.load_frame(cls, path, index, name, description, plane, hdulist_index, no_filter, fwhm,
                                   add_meta=add_meta, no_wcs=no_wcs, data_converter=data_converter)

        # Set the path
        mask.path = path

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @classmethod
    def empty_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(np.zeros(frame.shape), wcs=frame.wcs, pixelscale=frame.pixelscale)

    # -----------------------------------------------------------------

    @classmethod
    def full_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(np.ones(frame.shape), wcs=frame.wcs, pixelscale=frame.pixelscale)

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        return self._wcs

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self._wcs = value

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):

        """
        This function ...
        :return:
        """

        return self.wcs is not None

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        # Return the pixelscale of the WCS is WCS is defined
        if self.wcs is not None: return self.wcs.pixelscale
        else: return self._pixelscale  # return the pixelscale

    # -----------------------------------------------------------------

    @pixelscale.setter
    def pixelscale(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if not isinstance(value, Pixelscale): value = Pixelscale(value)
        self._pixelscale = value

    # -----------------------------------------------------------------

    @property
    def has_pixelscale(self):

        """
        Thisfunction ...
        :return:
        """

        return self.pixelscale is not None

    # -----------------------------------------------------------------

    @property
    def x_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.pixelscale.x

    # -----------------------------------------------------------------

    @property
    def y_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.pixelscale.y

    # -----------------------------------------------------------------

    @property
    def has_angular_pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: raise ValueError("No pixelscale")
        return types.is_angle(self.x_pixelscale)

    # -----------------------------------------------------------------

    @property
    def has_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: raise ValueError("No pixelscale")
        return types.is_length_quantity(self.x_pixelscale)

    # -----------------------------------------------------------------

    @property
    def average_pixelscale(self):

        """
        This function ...
        :return:
        """

        if self.wcs is not None: return self.wcs.average_pixelscale
        else: return self._pixelscale.average if self._pixelscale is not None else None

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

    def cropped(self, x_min, x_max, y_min, y_max, out_of_bounds="error"):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :param out_of_bounds:
        :return:
        """

        new = self.copy()
        new.crop(x_min, x_max, y_min, y_max, out_of_bounds=out_of_bounds)
        return new

    # -----------------------------------------------------------------

    def crop(self, x_min, x_max, y_min, y_max, out_of_bounds="error"):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :param out_of_bounds:
        :return:
        """

        # Crop the frame
        if out_of_bounds == "error": new_data = cropping.crop_check(self._data, x_min, x_max, y_min, y_max)
        elif out_of_bounds == "adjust": new_data, x_min, x_max, y_min, y_max = cropping.crop_direct(self._data, x_min, x_max, y_min, y_max)
        elif out_of_bounds == "expand": new_data = cropping.crop_absolute(self._data, x_min, x_max, y_min, y_max)
        else: raise ValueError("Invalid option for 'out_of_bounds'")

        # Adapt the WCS
        if self.wcs is not None:

            # Copy the current WCS
            new_wcs = self.wcs.copy()

            # Change the center pixel position
            new_wcs.wcs.crpix[0] -= x_min
            new_wcs.wcs.crpix[1] -= y_min

            # Change the number of pixels
            new_wcs.naxis1 = x_max - x_min
            new_wcs.naxis2 = y_max - y_min

            new_wcs._naxis1 = new_wcs.naxis1
            new_wcs._naxis2 = new_wcs.naxis2

        else: new_wcs = None

        # Check shape of data
        assert new_data.shape[1] == (x_max - x_min) and new_data.shape[0] == (y_max - y_min)

        # Replace the data and WCS
        self._data = new_data
        self._wcs = new_wcs

        # Return the limits
        return x_min, x_max, y_min, y_max

    # -----------------------------------------------------------------

    def cropped_to(self, region, factor=1, out_of_bounds="error"):

        """
        Ths function ...
        :param region:
        :param factor:
        :param out_of_bounds:
        :return:
        """

        new = self.copy()
        new.crop_to(region, factor=factor, out_of_bounds=out_of_bounds)
        return new

    # -----------------------------------------------------------------

    def crop_to(self, region, factor=1., out_of_bounds="error"):

        """
        This function ...
        :param region:
        :param factor:
        :param out_of_bounds:
        :return:
        """

        from ..region.rectangle import PixelRectangleRegion, SkyRectangleRegion

        # Pixel rectangle
        if isinstance(region, PixelRectangleRegion):
            if factor != 1: region = region * factor
            return self.crop(region.x_min_pixel, region.x_max_pixel, region.y_min_pixel, region.y_max_pixel, out_of_bounds=out_of_bounds)

        # Sky rectangle: to pixel rectangle
        elif isinstance(region, SkyRectangleRegion): return self.crop_to(region.to_pixel(self.wcs), factor=factor, out_of_bounds=out_of_bounds)

        # Other kind of shape
        else: return self.crop_to(region.bounding_box, factor=factor, out_of_bounds=out_of_bounds)

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=False, parallel=True, threshold=0.5, dilate=False, rank=2, connectivity=1,
              iterations=2):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :param threshold:
        :param dilate:
        :param rank:
        :param connectivity:
        :param iterations:
        :return:
        """

        from .frame import Frame

        # Check whether the frame has a WCS
        if not self.has_wcs: raise RuntimeError("Cannot rebin a mask without coordinate system")

        # Check whether the WCS is the same
        if self.wcs == reference_wcs: return Frame.ones_like(self)

        #from ..tools import plotting
        #print(self._data.shape)
        #plotting.plot_box(self._data.astype(int), title="before")

        # Calculate rebinned data and footprint of the original image
        if exact: new_data, footprint = reproject_exact((self._data.astype(int), self.wcs), reference_wcs, shape_out=reference_wcs.shape, parallel=parallel)
        else: new_data, footprint = reproject_interp((self._data.astype(int), self.wcs), reference_wcs, shape_out=reference_wcs.shape)

        #from ..tools import plotting
        #print(new_data.shape)
        #plotting.plot_box(new_data, title="after")

        # Get binary mask data
        mask_data = np.logical_or(new_data > threshold, np.isnan(new_data))

        # Replace the data and WCS
        self._data = mask_data
        self._wcs = reference_wcs.copy()

        # Dilate?
        if dilate: self.dilate_rc(rank, connectivity=connectivity, iterations=iterations)  # 1 also worked in test

        # Return the footprint
        return Frame(footprint, wcs=reference_wcs.copy())

    # -----------------------------------------------------------------

    def rebinned(self, reference_wcs, exact=False, parallel=True, threshold=0.5, dilate=False, rank=2, connectivity=1,
                 iterations=2):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :return:
        """

        new = self.copy()
        new.rebin(reference_wcs, exact=exact, parallel=parallel, threshold=threshold, dilate=dilate, rank=rank, connectivity=connectivity, iterations=iterations)
        return new

    # -----------------------------------------------------------------

    def downsample(self, factor, threshold=0.5, dilate=False, rank=2, connectivity=1, iterations=2):

        """
        This function ...
        :param factor:
        :param threshold:
        :param dilate:
        :param rank:
        :param connectivity:
        :param iterations:
        :return:
        """

        from scipy import ndimage
        data = ndimage.interpolation.zoom(self.data.astype(int), zoom=1.0/factor, order=0)
        new_xsize = data.shape[1]
        new_ysize = data.shape[0]

        # Set the data
        new_data = data > threshold

        # Dilate?
        if dilate: self.dilate_rc(rank, connectivity=connectivity, iterations=iterations)

        # Has coordinate system?
        if self.has_wcs: new_wcs = self.wcs.scaled(new_xsize, new_ysize)
        else: new_wcs = None

        # Set the new data and wcs
        self._data = new_data
        self._wcs = new_wcs

        # Return the new coordinate system
        return new_wcs

    # -----------------------------------------------------------------

    def downsampled(self, factor, threshold=0.5, dilate=False, rank=2, connectivity=1, iterations=2):

        """
        This function ...
        :param factor:
        :param threshold:
        :param dilate:
        :param rank:
        :param connectivity:
        :param iterations:
        :return:
        """

        new = self.copy()
        new.downsample(factor, threshold=threshold, dilate=dilate, rank=rank, connectivity=connectivity, iterations=iterations)
        return new

    # -----------------------------------------------------------------

    def get_mask(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        from ..region.region import PixelRegion, SkyRegion

        # Get mask
        if isinstance(region_or_mask, PixelRegion): mask = region_or_mask.to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, SkyRegion): mask = region_or_mask.to_pixel(self.wcs).to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, Mask): mask = region_or_mask
        else: raise ValueError("Argument must be region or mask")

        # Return
        return mask

    # -----------------------------------------------------------------

    def nmasked_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of masked pixels
        return np.sum(self.data[mask])

    # -----------------------------------------------------------------

    def relative_nmasked_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the relative number of masked pixels
        return float(np.sum(self.data[mask])) / np.sum(mask)

    # -----------------------------------------------------------------

    def nunmasked_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of unmasked pixels
        return np.sum(np.logical_not(self.data[mask]))

    # -----------------------------------------------------------------

    def relative_nunmasked_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the relative number of unmasked pixels
        return float(np.sum(np.logical_not(self.data[mask]))) / np.sum(mask)

    # -----------------------------------------------------------------

    def to_rgb(self, colour="black", background_color="white"):

        """
        This function ...
        :param colour:
        :param background_color:
        :return:
        """

        from .rgb import RGBImage
        return RGBImage.from_mask(self, colour=colour, background_color=background_color)

    # -----------------------------------------------------------------

    def to_rgba(self, colour="black"):

        """
        This function ...
        :param colour:
        :return:
        """

        from .rgba import RGBAImage
        return RGBAImage.from_mask(self, colour=colour)

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the mask ...")

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

        # ASDF format
        elif path.endswith(".asdf"): self.saveto_asdf(path, header=header, update_path=update_path)

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

        # Supported types: no binary (unit2, also doesn't exist in numpy) :(
        #from astropy.io.fits.hdu import DTYPE2BITPIX

        # Write to a FITS file
        # write_frame(self._data, header, path)
        #write_frame(self._data.astype(int), header, path)
        write_frame(self._data.astype(np.uint8), header, path)

        # Update the path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    def saveto_asdf(self, path, header=None, update_path=True):

        """
        This function ...
        :param path:
        :param header:
        :param update_path:
        :return:
        """

        # If a header is not specified, created it from the WCS
        if header is None: header = self.header

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

def union(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    rebin = kwargs.pop("rebin", False)

    # UNION = 0 + first + second + ... (0 is neutral element for sum)
    # so for one mask, union = 0 + mask = mask

    # Only one mask
    if len(args) == 1: return Mask(args[0])

    # REBIN?
    if rebin: args = rebin_to_highest_pixelscale(*args)

    #arrays = [arg.data for arg in args]
    arrays = []

    # Loop over the passed
    for arg in args:

        # Check type
        if isinstance(arg, MaskBase): arrays.append(arg.data)
        elif isinstance(arg, oldMask): arrays.append(arg)
        else: arrays.append(arg)

    # Create the union mask
    return Mask(np.sum(arrays, axis=0))

# -----------------------------------------------------------------

def intersection(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    rebin = kwargs.pop("rebin", False)

    # INTERSECTION = 1 * first * second * ... (1 is neutral element for multiplication)
    # so for one mask, intersection = 1 * mask = mask

    # Only one mask
    if len(args) == 1: return Mask(args[0])

    # REBIN?
    if rebin: args = rebin_to_highest_pixelscale(*args)

    #arrays = [arg.data for arg in args]
    arrays = []

    # Loop over the passed masks
    for arg in args:

        if isinstance(arg, MaskBase): arrays.append(arg.data)
        elif isinstance(arg, oldMask): arrays.append(arg)
        else: arrays.append(arg)

    # Create the intersection mask
    return Mask(np.product(arrays, axis=0))

# -----------------------------------------------------------------

def rebin_to_highest_pixelscale(*masks, **kwargs):

    """
    This function ...
    :param masks:
    :param kwargs:
    :return:
    """

    # Get mask names
    names = kwargs.pop("names", None)

    # In place?
    in_place = kwargs.pop("in_place", False)

    # Check
    if len(masks) == 1:

        # Success
        log.success("Only one mask: not rebinning")

        frame = masks[0]
        frame.name = names[0]
        return [frame]

    # Inform the user
    log.info("Rebinning masks to the coordinate system with the highest pixelscale ...")

    highest_pixelscale = None
    highest_pixelscale_wcs = None
    highest_pixelscale_index = None

    # Loop over the frames
    for index, mask in enumerate(masks):

        wcs = mask.wcs
        if wcs is None:

            if names is not None: raise ValueError("Coordinate system of the " + names[index] + " mask is not defined")
            else: raise ValueError("Coordinate system of the mask is not defined")

        if highest_pixelscale is None or wcs.average_pixelscale > highest_pixelscale:

            highest_pixelscale = wcs.average_pixelscale
            highest_pixelscale_wcs = wcs
            highest_pixelscale_index = index

    from ...core.tools.stringify import tostr

    # Debugging
    if names is not None: log.debug("The mask with the highest pixelscale is the '" + names[highest_pixelscale_index] + "' mask ...")
    log.debug("The highest pixelscale is " + tostr(highest_pixelscale))

    # Rebin
    return rebin_to_pixelscale(*masks, names=names, pixelscale=highest_pixelscale, wcs=highest_pixelscale_wcs, in_place=in_place)

# -----------------------------------------------------------------

def rebin_to_pixelscale(*masks, **kwargs):

    """
    THis function ...
    :param masks:
    :param kwargs:
    :return:
    """

    from ...core.tools.stringify import tostr

    # Get input
    names = kwargs.pop("names")
    highest_pixelscale = kwargs.pop("pixelscale")
    highest_pixelscale_wcs = kwargs.pop("wcs")

    # IN PLACE?
    in_place = kwargs.pop("in_place", False)

    # Initialize list for rebinned masks
    if in_place: new_masks = None
    else: new_masks = []

    # Rebin
    index = 0
    for mask in masks:

        # Determine mask name
        name = names[index] if names is not None else ""
        print_name = "'" + names[index] + "' " if names is not None else ""

        # If the current mask is the frame with the highest pixelscale
        if mask.wcs == highest_pixelscale_wcs:

            if names is not None: log.debug("Mask " + print_name + "has highest pixelscale of '" + tostr(highest_pixelscale) + "' and is not rebinned")

            # Not in place, create copy
            if not in_place:

                # Create new and set the name
                new = mask.copy()
                if names is not None: new.name = names[index]

                # Add
                new_masks.append(new)

        # The mask has a lower pixelscale, has to be rebinned
        else:

            # In place?
            if in_place: rebin_mask(name, mask, highest_pixelscale_wcs, in_place=True)

            # New masks
            else:

                # Create rebinned mask
                rebinned = rebin_mask(name, mask, highest_pixelscale_wcs)

                # Set the name
                if names is not None: rebinned.name = names[index]

                # Add the rebinned mask
                new_masks.append(rebinned)

        # Increment the index for the masks
        index += 1

    # Return the rebinned frames
    if not in_place: return new_masks

# -----------------------------------------------------------------

def rebin_mask(name, mask, wcs, in_place=False):

    """
    This function ...
    :param name:
    :param mask:
    :param wcs:
    :param in_place:
    :return:
    """

    # Debugging
    log.debug("Rebinning mask " + name + " ...")

    if in_place:
        mask.rebin(wcs)
        rebinned = None
    else: rebinned = mask.rebinned(wcs)

    # Return rebinned frame
    if not in_place: return rebinned

# -----------------------------------------------------------------
