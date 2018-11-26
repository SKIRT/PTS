#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.frame Contains the Frame class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import math
import numpy as np
import urllib
from scipy import ndimage
import tempfile

# Import astronomical modules
from reproject import reproject_exact, reproject_interp
from astropy.io import fits
from astropy.convolution import convolve, convolve_fft
from astropy.nddata import NDDataArray
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import Gaussian2DKernel
from astropy.units.core import UnitConversionError

# Import the relevant PTS classes and modules
from .cutout import Cutout
from ..basics.vector import Position
from ..region.rectangle import SkyRectangleRegion, PixelRectangleRegion
from ..basics.coordinate import SkyCoordinate, PixelCoordinate
from ..basics.stretch import SkyStretch
from ..tools import cropping
from ...core.basics.log import log
from ..basics.mask import MaskBase
from ...core.tools import filesystem as fs
from ...core.tools import archive
from ...core.units.unit import PhotometricUnit, parse_unit
from ...core.filter.filter import parse_filter
from .mask import Mask
from .alpha import AlphaMask
from ..convolution.kernels import get_fwhm, has_variable_fwhm
from ...core.tools import types
from ...core.units.parsing import parse_unit as u
from ...core.units.unit import get_conversion_factor
from ..basics.vector import PixelShape
from ...core.tools.stringify import tostr
from ...core.units.stringify import represent_unit
from ..basics.pixelscale import Pixelscale, PhysicalPixelscale, angular_or_physical_pixelscale
from ..basics.vector import Pixel
from ..region.region import PixelRegion, SkyRegion
from ...core.units.quantity import add_with_units, subtract_with_units, multiply_with_units, divide_with_units, get_value_and_unit
from ...core.tools.utils import create_lazified_class

# -----------------------------------------------------------------

class AllZeroError(Exception):

    """
    This class ...
    """

    def __init__(self, message):

        """
        Thisf unction ...
        :param message:
        """

        # Call the base class constructor with the parameters it needs
        super(AllZeroError, self).__init__(message)

# -----------------------------------------------------------------

nan_value = float("nan")
inf_value = float("inf")
zero_value = 0.0

# -----------------------------------------------------------------

def get_filter(frame_path):

    """
    Ths function allows getting the filter of a frame without loading the entire frame
    :param frame_path: 
    :return: 
    """

    from ..tools import headers
    header = fits.getheader(frame_path)
    fltr = headers.get_filter(fs.name(frame_path[:-5]), header)
    return fltr

# -----------------------------------------------------------------

def get_filter_name(frame_path):

    """
    This function ...
    :param frame_path:
    :return:
    """

    return str(get_filter(frame_path))

# -----------------------------------------------------------------

class Frame(NDDataArray):

    """
    This class ...
    """

    # Set the default extension
    default_extension = "fits"

    # -----------------------------------------------------------------

    def __init__(self, data, *args, **kwargs):

        """
        This function ...
        :param data:
        :param kwargs:
        """

        # Check the data
        #if data.dtype.type not in types.real_types:
        #    if data.dtype.type in types.integer_types: pass
        #    else: data = np.array(data, dtype=float)

        wcs = kwargs.pop("wcs", None)
        unit = kwargs.pop("unit", None)

        # Set general properties and flags
        self.name = kwargs.pop("name", None)
        self.description = kwargs.pop("description", None)
        self.zero_point = kwargs.pop("zero_point", None)
        self.source_extracted = kwargs.pop("source_extracted", False)
        self.extinction_corrected = kwargs.pop("extinction_corrected", False)
        self.sky_subtracted = kwargs.pop("sky_subtracted", False)

        # Set filter, go through setter
        self._filter = None
        self.filter = kwargs.pop("filter", None)

        # Set FWHM
        self._fwhm = kwargs.pop("fwhm", None)

        # Set pixelscale
        pixelscale = kwargs.pop("pixelscale", None)
        if pixelscale is not None: self._pixelscale = angular_or_physical_pixelscale(pixelscale)
        else: self._pixelscale = None

        # Set wavelength
        self._wavelength = kwargs.pop("wavelength", None)

        # Set meta data
        self.metadata = kwargs.pop("meta", dict())

        # Distance
        self.distance = kwargs.pop("distance", None)

        # PSF FILTER
        self._psf_filter = kwargs.pop("psf_filter", None)

        # The smoothing factor
        self.smoothing_factor = kwargs.pop("smoothing_factor", 1.)

        # The path
        self.path = kwargs.pop("path", None)
        self._from_multiplane = kwargs.pop("from_multiplane", False)

        # Call the constructor of the base class
        super(Frame, self).__init__(data, *args, **kwargs)

        # Set the WCS and unit
        self.wcs = wcs # go through the setter
        self.unit = unit # go through the setter

    # -----------------------------------------------------------------

    def set_meta(self, key, value):

        """
        This function ...
        :param key:
        :param value:
        :return:
        """

        self.metadata[key] = value

    # -----------------------------------------------------------------

    def get_meta(self, key, default=None):

        """
        This function ...
        :param key:
        :param default:
        :return:
        """

        if key in self.metadata: return self.metadata[key]
        elif default is not None: return default
        else: KeyError("No meta entry '" + key + "'")

    # -----------------------------------------------------------------

    @property
    def shape(self):
        return PixelShape.from_tuple(super(Frame, self).shape)

    # -----------------------------------------------------------------

    @property
    def filter_name(self):
        return str(self.filter) if self.filter is not None else None

    # -----------------------------------------------------------------

    @property
    def has_psf_filter(self):
        return self._psf_filter is not None

    # -----------------------------------------------------------------

    @property
    def psf_filter(self):

        """
        This function ...
        :return: 
        """

        # FILTER CAN BE DIFFERENT FROM PSF FILTER !! E.G. FUV IMAGE CONVOLVED TO RESOLUTION OF HERSCHEL BAND

        # BUT: if PSFFLTR WAS NOT FOUND IN HEADER, AND THUS PRESENT AS _PSF_FILTER, ASSUME PSF_FILTER = FILTER (ORIGINAL IMAGE)
        if self._psf_filter is None: return self.filter
        else: return self._psf_filter

    # -----------------------------------------------------------------

    @psf_filter.setter
    def psf_filter(self, value):
        if types.is_string_type(value): value = parse_filter(value)
        self._psf_filter = value

    # -----------------------------------------------------------------

    @property
    def psf_filter_name(self):
        return str(self.psf_filter) if self.psf_filter is not None else None

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return: 
        """

        # Find the FWHM for the filter
        if self._fwhm is not None: return self._fwhm
        elif self.psf_filter is not None and not has_variable_fwhm(self.psf_filter): return get_fwhm(self.psf_filter)
        else: return None

    # -----------------------------------------------------------------

    @fwhm.setter
    def fwhm(self, value):
        self._fwhm = value

    # -----------------------------------------------------------------

    @property
    def has_fwhm(self):
        return self.fwhm is not None

    # -----------------------------------------------------------------

    @property
    def data(self):
        return self._data

    # -----------------------------------------------------------------

    @NDDataArray.wcs.setter
    def wcs(self, wcs):
        self._wcs = wcs

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):
        return self.wcs is not None

    # -----------------------------------------------------------------

    @property
    def has_celestial_wcs(self):
        return self.has_wcs and self.wcs.is_celestial

    # -----------------------------------------------------------------

    @property
    def has_physical_wcs(self):
        return self.has_wcs and self.wcs.is_physical

    # -----------------------------------------------------------------

    @NDDataArray.unit.setter
    def unit(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        # Convert string units to PhotometricUnit or regular unit object
        if types.is_string_type(unit): unit = parse_unit(unit)

        # Set the unit
        self._unit = unit

    # -----------------------------------------------------------------

    @property
    def has_unit(self):
        return self.unit is not None

    # -----------------------------------------------------------------

    @property
    def is_photometric(self):
        return self.has_unit and isinstance(self.unit, PhotometricUnit)

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        return self.distance is not None

    # -----------------------------------------------------------------

    @property
    def filter(self):
        return self._filter

    # -----------------------------------------------------------------

    @filter.setter
    def filter(self, fltr):
        if fltr is None: self._filter = None
        else: self._filter = parse_filter(fltr)

    # -----------------------------------------------------------------

    @property
    def has_filter(self):
        return self.filter is not None

    # -----------------------------------------------------------------

    @property
    def is_broad_band(self):
        from ...core.filter.broad import BroadBandFilter
        return self.filter is not None and isinstance(self.filter, BroadBandFilter)

    # -----------------------------------------------------------------

    @property
    def is_narrow_band(self):
        from ...core.filter.narrow import NarrowBandFilter
        return self.filter is not None and isinstance(self.filter, NarrowBandFilter)

    # -----------------------------------------------------------------

    def __getitem__(self, item):

        """
        This function ...
        :param item:
        :return:
        """

        if isinstance(item, MaskBase): return self._data[item.data]
        elif isinstance(item, Pixel): return self._data[item.y, item.x]
        else: return self._data[item]

    # -----------------------------------------------------------------

    def __setitem__(self, item, value):

        """
        This function ...
        :param item:
        :return:
        """

        if isinstance(item, MaskBase): self._data[item.data] = value
        elif isinstance(item, Pixel): self._data[item.y, item.x] = value
        elif isinstance(item, tuple): self._data[item[0], item[1]] = value
        else: self._data[item] = value

    # -----------------------------------------------------------------

    def is_identical(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        if self.wcs != other.wcs: return False
        return np.all(self.data == other.data)

    # -----------------------------------------------------------------

    def is_close(self, other, rtol=1.e-5, atol=1.e-8, equal_nan=False):

        """
        This function ...
        :param other:
        :param rtol
        :param atol
        :param equal_nan:
        :return:
        """

        return np.all(np.isclose(self.data, other.data, rtol=rtol, atol=atol, equal_nan=equal_nan))

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self._data.__eq__(other)

    # -----------------------------------------------------------------

    def __ne__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self._data.__ne__(other)

    # -----------------------------------------------------------------

    def __gt__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self._data.__gt__(other)

    # -----------------------------------------------------------------

    def __ge__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self._data.__ge__(other)

    # -----------------------------------------------------------------

    def __lt__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self._data.__lt__(other)

    # -----------------------------------------------------------------

    def __le__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self._data.__le__(other)

    # -----------------------------------------------------------------

    @property
    def conversion_info(self):

        """
        Thisfunction ...
        :return:
        """

        info = dict()
        info["distance"] = self.distance
        info["pixelscale"] = self.pixelscale
        info["wavelength"] = self.wavelength
        return info

    # -----------------------------------------------------------------

    def __add__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.copy().__iadd__(value)

    # -----------------------------------------------------------------

    def __radd__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__add__(value)

    # -----------------------------------------------------------------

    def __iadd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_data = add_with_units(self.data, self.unit, other_value, other_unit=other_unit, conversion_info=self.conversion_info)

        # Replace the data
        self._data = new_data

        # Return
        return self

    # -----------------------------------------------------------------

    def __sub__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.copy().__isub__(value)

    # -----------------------------------------------------------------

    def __rsub__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new = self.copy()
        new._data = value - new._data
        return new

    # -----------------------------------------------------------------

    def __isub__(self, other):

        """
        This function ...
        :param value:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)
        if self.unit is not None and other_unit is None:
            log.warning("Subtraction with scalar value while the frame has a unit: assuming same units ...")
            other_unit = self.unit

        # Get new data
        new_data = subtract_with_units(self.data, self.unit, other_value, other_unit=other_unit, conversion_info=self.conversion_info)

        # Replace the data
        self._data = new_data

        # Return
        return self

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.copy().__imul__(value)

    # -----------------------------------------------------------------

    def __rmul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__mul__(value)

    # -----------------------------------------------------------------

    def __imul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_data, new_unit = multiply_with_units(self.data, self.unit, other_value, other_unit=other_unit)

        # Replace the data
        self._data = new_data

        # Replace the unit
        self._unit = new_unit

        # Return
        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.copy().__idiv__(value)

    # -----------------------------------------------------------------

    def __rdiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new = self.copy()
        new._data = value / new._data
        return new

    # -----------------------------------------------------------------

    def __idiv__(self, other):

        """
        This function ...
        :param value:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_data, new_unit = divide_with_units(self.data, self.unit, other_value, other_unit=other_unit)

        # Replace the data
        self._data = new_data

        # Replace the unit
        self._unit = new_unit

        # Return
        return self

    # -----------------------------------------------------------------

    __truediv__ = __div__

    # -----------------------------------------------------------------

    __itruediv__ = __idiv__

    # -----------------------------------------------------------------

    def __pow__(self, power):

        """
        This function ...
        :param power:
        :return:
        """

        return self.copy().__ipow__(power)

    # -----------------------------------------------------------------

    def __ipow__(self, power):

        """
        This function ...
        :param power:
        :return:
        """

        self._data **= power
        return self

    # -----------------------------------------------------------------

    def __abs__(self):

        """
        This function ...
        :return:
        """

        new = self.copy()
        new._data = abs(self._data)
        return new

    # -----------------------------------------------------------------

    def log10(self):

        """
        This function ...
        :return: 
        """

        # Replace the data
        self._data = np.log10(self.data)

        # Set the unit to None
        self.unit = None

    # -----------------------------------------------------------------

    def get_log10(self):

        """
        This function ...
        :return: 
        """

        new = self.copy()
        new.log10()
        return new

    # -----------------------------------------------------------------

    def log(self):

        """
        This function ...
        :return: 
        """

        # Replace the data
        self._data = np.log(self.data)

        # Set the unit to None
        self.unit = None

    # -----------------------------------------------------------------

    def get_log(self):

        """
        THis function ...
        :return: 
        """

        new = self.copy()
        new.log10()
        return new

    # -----------------------------------------------------------------

    @classmethod
    def from_url(cls, url, index=None, name=None, description=None, plane=None, hdulist_index=None, no_filter=False,
                  fwhm=None, add_meta=True, distance=None):

        """
        This function ...
        :param url:
        :param index:
        :param name:
        :param description:
        :param plane:
        :param hdulist_index:
        :param no_filter:
        :param fwhm:
        :param add_meta:
        :param distance:
        :return:
        """

        # Inform the user
        log.info("Downloading file " + url + " ...")

        # Local path
        temp_path = tempfile.gettempdir()
        filename = fs.name(url)
        local_path = fs.join(temp_path, filename)

        # Download
        urllib.urlretrieve(url, local_path)

        if local_path.endswith(".fits"): fits_path = local_path
        else:

            # Inform the user
            log.info("Decompressing kernel file ...")

            # Fits path
            fits_path = fs.join(temp_path, fs.strip_extension(filename))

            # Decompress the kernel FITS file
            archive.decompress_file(local_path, fits_path)

            # Remove the compressed file
            fs.remove_file(local_path)

        # Open the FITS file
        return cls.from_file(fits_path, index, name, description, plane, hdulist_index, no_filter, fwhm, add_meta, distance=distance)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, index=None, name=None, description=None, plane=None, hdulist_index=None, no_filter=False,
                  fwhm=None, add_meta=True, extra_meta=None, silent=False, distance=None, no_wcs=False, density=False,
                  brightness=False, density_strict=False, brightness_strict=False, wcs=None, pixelscale=None):

        """
        This function ...
        :param path:
        :param index:
        :param name:
        :param description:
        :param plane:
        :param hdulist_index: if None, is automatically decided based on where the imageHDU is.
        :param no_filter:
        :param fwhm:
        :param add_meta:
        :param extra_meta:
        :param silent:
        :param distance:
        :param no_wcs:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param wcs:
        :param pixelscale:
        :return:
        """

        # Show which image we are importing
        if not silent: log.info("Reading in file '" + path + "' ...")

        from ..core.fits import load_frame, DamagedFITSFileError
        # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
        try: frame = load_frame(cls, path, index, name, description, plane, hdulist_index, no_filter, fwhm,
                               add_meta=add_meta, extra_meta=extra_meta, distance=distance, no_wcs=no_wcs,
                               density=density, brightness=brightness, density_strict=density_strict,
                               brightness_strict=brightness_strict)
        except DamagedFITSFileError: raise IOError("File is possibly damaged")

        # Set
        if wcs is not None: frame.wcs = wcs
        if pixelscale is not None: frame.pixelscale = pixelscale

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def apply_mask(self, mask, fill=0.0, invert=False, return_complement=False):

        """
        This function ...
        :param mask:
        :param fill:
        :param invert:
        :param return_complement:
        """

        if invert:
            mask = mask.inverse()
            inverse_mask = mask
        elif return_complement: inverse_mask = mask.inverse()
        else: inverse_mask = None

        # Apply the mask
        self[mask] = fill

        # Return the complementary frame
        if return_complement: return self.applied_mask(inverse_mask, fill=fill)

    # -----------------------------------------------------------------

    def apply_mask_nans(self, mask, invert=False, return_complement=False):

        """
        This function ...
        :param mask:
        :param invert:
        :param return_complement:
        :return:
        """

        return self.apply_mask(mask, fill=nan_value, invert=invert, return_complement=return_complement)

    # -----------------------------------------------------------------

    def apply_mask_infs(self, mask, invert=False, return_complement=False):

        """
        This function ...
        :param mask:
        :param invert:
        :param return_complement:
        :return:
        """

        return self.apply_mask(mask, fill=inf_value, invert=invert, return_complement=return_complement)

    # -----------------------------------------------------------------

    def applied_mask(self, mask, fill=0.0, invert=False):

        """
        This function ...
        :param mask:
        :param fill:
        :param invert:
        :return:
        """

        new = self.copy()
        new.apply_mask(mask, fill=fill, invert=invert)
        return new

    # -----------------------------------------------------------------

    def applied_mask_nans(self, mask, invert=False):

        """
        This function ...
        :param mask:
        :param invert:
        :return:
        """

        return self.applied_mask(mask, fill=nan_value, invert=invert)

    # -----------------------------------------------------------------

    def applied_mask_infs(self, mask, invert=False):

        """
        This function ...
        :param mask:
        :param invert:
        :return:
        """

        return self.applied_mask(mask, fill=inf_value, invert=invert)

    # -----------------------------------------------------------------

    def apply_alpha_mask(self, mask, invert=False):

        """
        This function ...
        :param mask:
        :param invert:
        :return:
        """

        data = mask.as_real()
        if invert: data = 1. - data

        # Multiply the data with the alpha mask
        self._data *= data

    # -----------------------------------------------------------------

    def applied_alpha_mask(self, mask, invert=False):

        """
        This function ...
        :param mask:
        :param invert:
        :return:
        """

        new = self.copy()
        new.apply_alpha_mask(mask, invert=invert)
        return new

    # -----------------------------------------------------------------

    @property
    def all_zero(self):

        """
        This function ...
        :return:
        """

        return np.all(np.equal(self._data, 0))

    # -----------------------------------------------------------------

    @property
    def all_nonzero(self):

        """
        This function ...
        :return:
        """

        return not np.any(np.equal(self._data, 0))

    # -----------------------------------------------------------------

    def where(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Mask(np.equal(self._data, value), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def where_not(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Mask(np.not_equal(self._data, value), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def where_smaller_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Mask(np.less(self._data, value), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def where_smaller_than_or_equal(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Mask(np.less_equal(self._data, value), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def where_greater_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Mask(np.greater(self._data, value), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def where_greater_than_or_equal(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return Mask(np.greater_equal(self._data, value), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    @property
    def absolute(self):
        return Frame(np.abs(self.data), wcs=self.wcs, unit=self.unit, pixelscale=self.pixelscale, filter=self.filter)

    # -----------------------------------------------------------------

    def make_absolute(self):

        """
        This function ...
        :return:
        """

        self._data = np.abs(self._data)

    # -----------------------------------------------------------------
    # NANS
    # -----------------------------------------------------------------

    @property
    def nans(self):
        return Mask(np.isnan(self.data), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    @property
    def nans_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.where(self.nans))]

    # -----------------------------------------------------------------

    @property
    def nnans(self):
        return np.sum(self.nans.data)

    # -----------------------------------------------------------------

    @property
    def relative_nnans(self):
        return float(self.nnans) / self.npixels

    # -----------------------------------------------------------------

    @property
    def has_nans(self):
        return np.any(self.nans.data)

    # -----------------------------------------------------------------

    @property
    def all_nans(self):
        return np.all(self.nans.data)

    # -----------------------------------------------------------------
    # INFS
    # -----------------------------------------------------------------

    @property
    def infs(self):
        return Mask(np.isinf(self.data), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    @property
    def infs_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.where(self.infs))]

    # -----------------------------------------------------------------

    @property
    def all_infs(self):
        return np.all(self.infs.data)

    # -----------------------------------------------------------------

    @property
    def ninfs(self):
        return np.sum(self.infs.data)

    # -----------------------------------------------------------------

    @property
    def relative_ninfs(self):
        return float(self.ninfs) / self.npixels

    # -----------------------------------------------------------------

    @property
    def has_infs(self):
        return np.any(self.infs.data)

    # -----------------------------------------------------------------
    # INVALID
    # -----------------------------------------------------------------

    @property
    def invalid(self):
        return self.nans + self.infs

    # -----------------------------------------------------------------

    @property
    def invalid_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.where(self.invalid))]

    # -----------------------------------------------------------------

    @property
    def all_invalid(self):
        return np.all(self.invalid.data)

    # -----------------------------------------------------------------

    @property
    def ninvalid(self):
        return np.sum(self.invalid.data)

    # -----------------------------------------------------------------

    @property
    def relative_ninvalid(self):
        return float(self.ninvalid) / self.npixels

    # -----------------------------------------------------------------

    @property
    def has_invalid(self):
        return np.any(self.invalid.data)

    # -----------------------------------------------------------------
    # VALID
    # -----------------------------------------------------------------

    @property
    def valid(self):
        return self.invalid.inverse()

    # -----------------------------------------------------------------

    @property
    def valid_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.where(self.valid))]

    # -----------------------------------------------------------------

    @property
    def all_valid(self):
        return np.all(self.valid.data)

    # -----------------------------------------------------------------

    @property
    def nvalid(self):
        return np.sum(self.valid.data)

    # -----------------------------------------------------------------

    @property
    def relative_nvalid(self):
        return float(self.nvalid) / self.npixels

    # -----------------------------------------------------------------

    @property
    def has_valid(self):
        return np.any(self.valid.data)

    # -----------------------------------------------------------------
    # ZEROES
    # -----------------------------------------------------------------

    @property
    def zeroes(self):
        return self.where(zero_value)

    # -----------------------------------------------------------------

    @property
    def zeroes_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.where(self.zeroes))]

    # -----------------------------------------------------------------

    @property
    def all_zeroes(self):
        return np.all(self.zeroes.data)

    # -----------------------------------------------------------------

    @property
    def nzeroes(self):
        return np.sum(self.zeroes.data)

    # -----------------------------------------------------------------

    @property
    def relative_nzeroes(self):
        return float(self.nzeroes) / self.npixels

    # -----------------------------------------------------------------

    @property
    def has_zeroes(self):
        return np.any(self.zeroes.data)

    # -----------------------------------------------------------------
    # NON-ZEROES
    # -----------------------------------------------------------------

    @property
    def nonzeroes(self):
        return self.where_not(0.0)

    # -----------------------------------------------------------------

    @property
    def all_nonzeroes(self):
        return np.all(self.nonzeroes.data)

    # -----------------------------------------------------------------

    @property
    def zeroes_x(self):
        return np.where(self.zeroes)[1]

    # -----------------------------------------------------------------

    @property
    def zeroes_y(self):
        return np.where(self.zeroes)[0]

    # -----------------------------------------------------------------

    @property
    def nonzeroes_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.nonzero(self._data))]

    # -----------------------------------------------------------------

    @property
    def nnonzeroes(self):
        return np.sum(self.nonzeroes.data)

    # -----------------------------------------------------------------

    @property
    def relative_nnonzeroes(self):
        return float(self.nnonzeroes) / self.npixels

    # -----------------------------------------------------------------

    @property
    def nonzeroes_x(self):
        return np.nonzero(self._data)[1]

    # -----------------------------------------------------------------

    @property
    def nonzeroes_y(self):
        return np.nonzero(self._data)[0]

    # -----------------------------------------------------------------
    # NEGATIVES
    # -----------------------------------------------------------------

    @property
    def negatives(self):
        return self.where_smaller_than(zero_value)

    # -----------------------------------------------------------------

    @property
    def negatives_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.where(self.negatives))]

    # -----------------------------------------------------------------

    @property
    def nnegatives(self):
        return np.sum(self.negatives.data)

    # -----------------------------------------------------------------

    @property
    def relative_nnegatives(self):
        return float(self.nnegatives) / self.npixels

    # -----------------------------------------------------------------

    @property
    def has_negatives(self):
        return np.any(self.negatives.data)

    # -----------------------------------------------------------------

    @property
    def all_negatives(self):
        return np.all(self.negatives.data)

    # -----------------------------------------------------------------
    # POSITIVES
    # -----------------------------------------------------------------

    @property
    def positives(self):
        return self.where_greater_than(zero_value)

    # -----------------------------------------------------------------

    @property
    def positives_pixels(self):
        return [Pixel(x, y) for y, x in np.transpose(np.where(self.positives))]

    # -----------------------------------------------------------------

    @property
    def npositives(self):
        return np.sum(self.positives.data)

    # -----------------------------------------------------------------

    @property
    def relative_npositives(self):
        return float(self.npositives) / self.npixels

    # -----------------------------------------------------------------

    @property
    def has_positives(self):
        return np.any(self.positives.data)

    # -----------------------------------------------------------------

    @property
    def all_positives(self):
        return np.all(self.positives.data)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def values_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        if isinstance(region_or_mask, PixelRegion): mask = region_or_mask.to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, SkyRegion):
            if not self.has_wcs: raise ValueError("Cannot specify a sky region when frame has no WCS")
            mask = region_or_mask.to_pixel(self.wcs).to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, Mask): mask = region_or_mask
        else: raise ValueError("Argument must be pixel region or mask")

        # Return the values as a Numpy array
        return self.data[mask.data]

    # -----------------------------------------------------------------

    def sum_in(self, region_or_mask, add_unit=False):

        """
        Thisn function ...
        :param region_or_mask:
        :param add_unit:
        :return:
        """

        # Get the values
        values = self.values_in(region_or_mask)

        # Return the nansum
        result = np.nansum(values)

        # Add unit?
        if add_unit and self.has_unit:
            #if self.unit.is_brightness: log.warning("Unit is a surface brightness: adding all pixel values may not be useful before a conversion to a non-brightness unit")
            #if self.is_per_angular_or_intrinsic_area: log.warning("Unit is per angular or intrinsic are")
            if self.is_per_angular_or_intrinsic_area:
                if self.is_per_angular_area: log.warning("Unit is per angular area: the result of adding pixel values may not be useful before a conversion to a non-intensity or brightness unit")
                elif self.is_per_intrinsic_area: log.warning("Unit is per physical area: the result of adding pixel values may not be useful before a conversion to a non-brightness unit")
            return result * self.unit
        else: return result

    # -----------------------------------------------------------------

    def quadratic_sum_in(self, region_or_mask, add_unit=False):

        """
        This function ...
        :param region_or_mask:
        :param add_unit:
        :return:
        """

        # Get the values
        values = self.values_in(region_or_mask)

        # Get non-nans
        nans = np.isnan(values)
        not_nans = np.logical_not(nans)

        # Get the quadratic sum
        result = np.sqrt(np.sum(values[not_nans] ** 2))

        # Add unit?
        if add_unit and self.has_unit:
            #if self.unit.is_brightness: log.warning("Unit is a surface brightness: adding all pixel values may not be useful before a conversion to a non-brightness unit")
            if self.is_per_angular_or_intrinsic_area:
                if self.is_per_angular_area: log.warning("Unit is per angular area: the result of adding pixel values may not be useful before a conversion to a non-intensity or brightness unit")
                elif self.is_per_intrinsic_area: log.warning("Unit is per physical area: the result of adding pixel values may not be useful before a conversion to a non-brightness unit")
            return result * self.unit
        else: return result

    # -----------------------------------------------------------------

    @property
    def values(self):
        #return self.data.flat USE THIS?
        return self.data.flatten()

    # -----------------------------------------------------------------

    @property
    def min(self):
        return np.nanmin(self.data)

    # -----------------------------------------------------------------

    @property
    def max(self):
        return np.nanmax(self.data)

    # -----------------------------------------------------------------

    def get_mask(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        if isinstance(region_or_mask, PixelRegion): mask = region_or_mask.to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, SkyRegion): mask = region_or_mask.to_pixel(self.wcs).to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, Mask): mask = region_or_mask
        else: raise ValueError("Argument must be region or mask")

        # Return
        return mask

    # -----------------------------------------------------------------

    def min_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return minimum value in the masked values
        return np.nanmin(self.data[mask])

    # -----------------------------------------------------------------

    def max_in(self, region_or_mask):

        """
        Thisfunction ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return maximum value in the masked values
        return np.nanmax(self.data[mask])

    # -----------------------------------------------------------------

    def nnans_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of nan pixels
        #return np.sum(np.equal(self.data[mask], nan_value))

        # 2nd attempt
        #array = self.data[mask]
        #from .mask import union
        #return np.sum(union(*[np.equal(array, value) for value in nan_values]))

        return np.sum(np.isnan(self.data[mask]))

    # -----------------------------------------------------------------

    def relative_nnans_in(self, region_or_mask):

        """
        Thisf function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return
        return float(self.nnans_in(mask)) / np.sum(mask)

    # -----------------------------------------------------------------

    def ninfs_in(self, region_or_mask):

        """
        Thisn function ...
        :param region_or_mask:
        :return:
        """

        # get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of nan pixels
        #return np.sum(np.equal(self.data[mask], inf_values[0]) + np.equal(self.data[mask], inf_values[1]))

        # 2nd attempt
        #array = self.data[mask]
        #from .mask import union
        #return np.sum(union(*[np.equal(array, value) for value in inf_values]))

        return np.sum(np.isinf(self.data[mask]))

    # -----------------------------------------------------------------

    def relative_ninfs_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return
        return float(self.ninfs_in(mask)) / np.sum(mask)

    # -----------------------------------------------------------------

    def nzeroes_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of zero pixels
        return np.sum(np.equal(self.data[mask], zero_value))

    # -----------------------------------------------------------------

    def relative_nzeroes_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return
        return float(self.nzeroes_in(mask)) / np.sum(mask)

    # -----------------------------------------------------------------

    def nnegatives_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of negative pixels
        return np.sum(np.less(self.data[mask], zero_value))

    # -----------------------------------------------------------------

    def relative_nnegatives_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return
        return float(self.nnegatives_in(mask)) / np.sum(mask)

    # -----------------------------------------------------------------

    def npositives_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of positive pixels
        return np.sum(np.greater(self.data[mask], zero_value))

    # -----------------------------------------------------------------

    def relative_npositives_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return
        return float(self.npositives_in(region_or_mask)) / np.sum(mask)

    # -----------------------------------------------------------------

    @property
    def npixels(self):

        """
        This function ...
        :return:
        """

        return self.xsize * self.ysize

    # -----------------------------------------------------------------

    def npixels_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Return the number of pixels
        return np.sum(mask)

    # -----------------------------------------------------------------

    @property
    def unique_values(self):

        """
        This function ...
        :return:
        """

        # Loop over the unique values in the gridded data
        values = np.unique(self._data)
        return list(values)

    # -----------------------------------------------------------------

    @property
    def unique_value_masks(self):

        """
        This function ...
        :return:
        """

        masks = []

        for value in self.unique_values:
            where = self.where(value)
            masks.append(where)

        return masks

    # -----------------------------------------------------------------

    @property
    def unique_values_and_masks(self):

        """
        This function ...
        :return:
        """

        returns = []
        for value in self.unique_values:
            where = self.where(value)
            returns.append((value, where))

        return returns

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

        # Convert to pixelscale
        pixelscale = angular_or_physical_pixelscale(value)

        # Set the internal pixelscale
        self._pixelscale = pixelscale

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
    def angular_pixelscale(self):

        """
        Thisn function ...
        :return:
        """

        if self.has_angular_pixelscale: return self.pixelscale
        elif self.has_physical_pixelscale:
            if not self.has_distance: raise ValueError("Distance is not defined to convert physical pixelscale into angular pixelscale")
            return self.pixelscale.to_angular(distance=self.distance)

        # No pixelscale
        else: return None

    # -----------------------------------------------------------------

    @property
    def physical_pixelscale(self):

        """
        Thisf unction ...
        :return:
        """

        if self.has_physical_pixelscale: return self.pixelscale
        elif self.has_angular_pixelscale:
            if not self.has_distance: raise ValueError("Distance is not defined to convert angular pixelscale into physical pixelscale")
            return self.pixelscale.to_physical(distance=self.distance)

        # No pixelscale
        else: return None

    # -----------------------------------------------------------------

    @property
    def average_angular_pixelscale(self):

        """
        This function ...
        :return:
        """

        angular = self.angular_pixelscale
        if angular is None: return None
        else: return angular.average

    # -----------------------------------------------------------------

    @property
    def average_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        physical = self.physical_pixelscale
        if physical is None: return None
        else: return physical.average

    # -----------------------------------------------------------------

    @property
    def pixelarea(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: return None

        try:

            if self.has_angular_pixelscale:
                solid = (self.pixelscale.x * self.pixelscale.y).to("sr")
                return solid
            elif self.has_physical_pixelscale:
                area = (self.pixelscale.x * self.pixelscale.y).to("pc2")
                return area
            else: raise RuntimeError("Something went wrong")

        except UnitConversionError:

            log.warning("Unit conversion error: " + str(self.pixelscale.x) + ", " + str(self.pixelscale.y))

            # NO UNITS? ASSUME IT'S IN DEGREES? IS THIS ANY GOOD? WHY WAS THIS DONE?
            deg = (self.pixelscale.x * self.pixelscale.y).value
            return (deg * u("deg")).to("sr")

    # -----------------------------------------------------------------

    @property
    def angular_pixelarea(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: return None
        return (self.angular_pixelscale.x * self.angular_pixelscale.y).to("sr")

    # -----------------------------------------------------------------

    @property
    def pixel_solid_angle(self):
        return self.angular_pixelarea

    # -----------------------------------------------------------------

    @property
    def physical_pixelarea(self):

        """
        This function ...
        :return:
        """

        if not self.has_pixelscale: return None
        return (self.physical_pixelscale.x * self.physical_pixelscale.y).to("pc2")

    # -----------------------------------------------------------------

    @property
    def pixel_area(self):

        """
        This function ...
        :return:
        """

        return self.physical_pixelarea

    # -----------------------------------------------------------------

    @classmethod
    def example(cls, min_value=0., max_value=1.):

        """
        This function ...
        :return:
        """

        # Generate a shape
        xsize = np.random.randint(250, 1000)
        ysize = np.random.randint(250, 1000)
        shape = (ysize, xsize)

        # Generate random
        frame = cls.random(shape, min_value=min_value, max_value=max_value)

        # Generate random unit
        unit = np.random.choice(["Jy", "W/micron", "W/m2", "count/pix", "Lsun/pc2", "erg/s/micron/arcsec2"])

        # Generate a random wavelength
        wavelength_micron = np.random.uniform(0.1,1000)
        wavelength = wavelength_micron * u("micron")

        # Generate random pixelscale
        pixelscale_arcsec = np.random.uniform(0.1,20.)
        pixelscale = pixelscale_arcsec * u("arcsec")

        # Generate a random FWHM
        fwhm_arcsec = np.random.uniform(0.5*pixelscale_arcsec, 10.*pixelscale_arcsec)
        fwhm = fwhm_arcsec * u("arcsec")

        # Set the properties
        frame.unit = unit
        frame.wavelength = wavelength
        frame.pixelscale = pixelscale
        frame.fwhm = fwhm

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @classmethod
    def random(cls, shape, min_value=0., max_value=1., **kwargs):

        """
        This function ...
        :param shape:
        :param min_value:
        :param max_value:
        :param kwargs:
        :return: 
        """

        # Create a new frame
        data = (max_value - min_value) * np.random.random(shape) + min_value
        new = cls(data, **kwargs)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def random_normal(cls, shape, mean=0.0, sigma=1., **kwargs):

        """
        This function ...
        :param shape:
        :param mean:
        :param sigma:
        :param kwargs:
        :return:
        """

        data = sigma * np.random.randn(shape[0], shape[1]) + mean
        new = cls(data, **kwargs)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def random_integers(cls, shape, min_value, max_value, **kwargs):

        """
        This function ...
        :param shape:
        :param min_value:
        :param max_value:
        :param kwargs:
        :return:
        """

        data = np.random.random_integers(min_value, max_value, size=shape)
        new = cls(data, **kwargs)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def ones(cls, shape, wcs=None, filter=None, unit=None):

        """
        This function ...
        :param shape: 
        :param wcs: 
        :param filter: 
        :return: 
        """

        # Create a new frame
        new = cls(np.ones(shape), wcs=wcs, filter=filter, unit=unit)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def initialize_ones(cls, *args, **kwargs):
        return cls.ones(*args, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def zeros(cls, shape, wcs=None, filter=None, unit=None):

        """
        This function ...
        :param shape:
        :param wcs:
        :param filter:
        :return:
        """

        # Create a new frame
        new = cls(np.zeros(shape), wcs=wcs, filter=filter, unit=unit)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def initialize_zeroes(cls, *args, **kwargs):
        return cls.zeros(*args, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def ones_like(cls, frame):

        """
        This function ...
        :param frame: 
        :return: 
        """

        if hasattr(frame, "wcs"): wcs = frame.wcs
        else: wcs = None

        # Create a new frame
        if hasattr(frame, "_data"): new = cls(np.ones_like(frame._data), wcs=wcs)
        else: new = cls(np.zeros_like(frame), wcs=wcs)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def zeros_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        if hasattr(frame, "wcs"): wcs = frame.wcs
        else: wcs = None

        # Create a new frame
        if hasattr(frame, "_data"): new = cls(np.zeros_like(frame._data), wcs=wcs)
        else: new = cls(np.zeros_like(frame), wcs=wcs)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def initialize_nans(cls, shape, wcs=None, filter=None, unit=None):

        """
        This function ...
        :param shape:
        :param wcs:
        :param filter:
        :param unit:
        :return:
        """

        # Create a new frame
        new = cls(np.full(shape, nan_value), wcs=wcs, filter=filter, unit=unit)
        return new

    # -----------------------------------------------------------------

    @classmethod
    def nans_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Return a NaN-filled copy of the frame
        nans = cls.zeros_like(frame)
        nans.fill(np.nan)
        return nans

    # -----------------------------------------------------------------

    @classmethod
    def filled_like(cls, frame, value):

        """
        This function ...
        :param frame:
        :param value:
        :return:
        """

        frame = cls.zeros_like(frame)
        frame.fill(value)
        return frame

    # -----------------------------------------------------------------

    @property
    def is_constant(self):

        """
        This function ...
        :return:
        """

        return np.nanmax(self._data) == np.nanmin(self._data)

    # -----------------------------------------------------------------

    def is_constant_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        # Get mask
        if isinstance(region_or_mask, PixelRegion): mask = region_or_mask.to_mask(self.xsize, self.ysize)
        elif isinstance(region_or_mask, Mask): mask = region_or_mask
        else: raise ValueError("Argument must be pixel region or mask")

        # Return whether constant
        return np.nanmax(self.data[mask]) == np.nanmin(self.data[mask])

    # -----------------------------------------------------------------

    @property
    def xsize(self): return self.shape.x

    # -----------------------------------------------------------------

    @property
    def ysize(self): return self.shape.y

    # -----------------------------------------------------------------

    @property
    def fwhm_pix(self):

        """
        This function ...
        :return:
        """

        return (self.fwhm / self.average_pixelscale).to("").value if self.fwhm is not None else None

    # -----------------------------------------------------------------

    @property
    def sigma(self):

        """
        THis function ...
        :return:
        """

        from ..tools import statistics
        return self.fwhm * statistics.fwhm_to_sigma if self.fwhm is not None else None

    # -----------------------------------------------------------------

    @property
    def sigma_pix(self):

        """
        Thisfunction ...
        :return:
        """

        from ..tools import statistics
        return self.fwhm_pix * statistics.fwhm_to_sigma if self.fwhm is not None else None

    # -----------------------------------------------------------------

    @property
    def header(self):

        """
        This function ...
        """

        # If the WCS for this frame is defined, use it to create a header
        if self.wcs is not None: header = self.wcs.to_header()

        # Else, create a new empty header
        else: header = fits.Header()

        # Add properties to the header
        header['NAXIS'] = 2
        header['NAXIS1'] = self.xsize
        header['NAXIS2'] = self.ysize

        # ISSUE: see bug #4592 on Astropy GitHub (WCS.to_header issue)
        # temporary fix !!
        # I don't know whether this is a good fix.. but it seems to fix it for a particular situation
        #if "PC1_1" in header:

            #if "NAXIS1" in header: header.remove("NAXIS1")
            #if "NAXIS2" in header: header.remove("NAXIS2")
            #if "CDELT1" in header: header.remove("CDELT1")
            #if "CDELT2" in header: header.remove("CDELT2")
            #header.rename_keyword("PC1_1", "CD1_1")
            #header.rename_keyword("PC2_2", "CD2_2")

        # Return the header
        return header

    # -----------------------------------------------------------------

    @property
    def wavelength(self):

        """
        This function ...
        :return:
        """

        # First return the wavelength, if set
        if self._wavelength is not None: return self._wavelength
        elif self.has_filter: return self.filter.wavelength
        else: return None

    # -----------------------------------------------------------------

    @wavelength.setter
    def wavelength(self, value):

        """
        This function ...
        :return:
        """

        self._wavelength = value

    # -----------------------------------------------------------------

    @property
    def wavelength_micron(self):

        """
        This function ...
        :return:
        """

        wavelength = self.wavelength
        if wavelength is None: return None
        else: return wavelength.to("micron").value

    # -----------------------------------------------------------------

    @property
    def pivot_wavelength_or_wavelength(self):

        """
        This fucntion ...
        :return: 
        """

        if self.filter is not None: return self.filter.pivot
        else: return self.wavelength
        
    # -----------------------------------------------------------------

    def cutout_around(self, position, radius, as_frame=False):

        """
        This function ...
        :param position:
        :param radius:
        :param as_frame:
        :return:
        """

        # Create cutout
        cutout = Cutout.cutout(self, position, radius)

        # Return as frame
        if as_frame:

            frame = Frame(np.asarray(cutout), name=self.name, pixelscale=self.pixelscale, fwhm=self.fwhm, distance=self.distance)
            frame.set_meta("xmin", cutout.x_min)
            frame.set_meta("xmax", cutout.x_max)
            frame.set_meta("ymin", cutout.y_min)
            frame.set_meta("ymax", cutout.y_max)
            return frame

        # Return cutout
        else: return cutout

    # -----------------------------------------------------------------

    def flip_horizontally(self):

        """
        This function ...
        :return:
        """

        self._data = np.fliplr(self._data)

    # -----------------------------------------------------------------

    def flip_vertically(self):

        """
        This function ...
        :return:
        """

        self._data = np.flipud(self._data)

    # -----------------------------------------------------------------

    def convert_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False,
                   brightness_strict=False, wavelength=None, silent=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param wavelength:
        :param silent:
        :return:
        """

        # Check the unit of this frame is defined
        if not self.has_unit: raise ValueError("The unit of the frame is not defined")

        # Parse "to unit": VERY IMPORTANT, BECAUSE DOING SELF.UNIT = TO_UNIT WILL OTHERWISE REPARSE AND WILL BE OFTEN INCORRECT!! (NO DENSITY OR BRIGHTNESS INFO)
        #to_unit = PhotometricUnit(to_unit, density=density, brightness=brightness, brightness_strict=brightness_strict, density_strict=density_strict)
        to_unit = parse_unit(to_unit, density=density, brightness=brightness, brightness_strict=brightness_strict, density_strict=density_strict)

        # Already in the correct unit
        if to_unit == self.unit:
            if not silent: log.debug("Frame is already in the desired unit")
            return 1.

        # Debugging
        if not silent: log.debug("Converting the frame from unit " + tostr(self.unit, add_physical_type=True) + " to unit " + tostr(to_unit, add_physical_type=True) + " ...")

        # Get the conversion factor
        factor = self._get_conversion_factor(to_unit, distance=distance, wavelength=wavelength, silent=silent)

        # Debugging
        if not silent: log.debug("Conversion factor: " + str(factor))

        # Convert
        self.convert_by_factor(factor, to_unit)

        # Return the conversion factor
        return factor

    # -----------------------------------------------------------------

    def _get_conversion_factor(self, to_unit, distance=None, wavelength=None, silent=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param wavelength:
        :param silent:
        :return:
        """

        # New: one central place to implement this
        return get_conversion_factor(self.unit, to_unit, parse=False, silent=silent, distance=distance,
                                     wavelength=wavelength, conversion_info=self.conversion_info)

    # -----------------------------------------------------------------

    def converted_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False,
                     brightness_strict=False, wavelength=None, silent=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param wavelength:
        :param silent:
        :return:
        """

        new = self.copy()
        new.convert_to(to_unit, distance=distance, density=density, brightness=brightness, density_strict=density_strict,
                       brightness_strict=brightness_strict, wavelength=wavelength, silent=silent)
        return new

    # -----------------------------------------------------------------

    def convert_by_factor(self, factor, new_unit):

        """
        This function SHOULD BE USED WITH CARE! UNIT CONVERSION IS LEFT UP TO THE USER
        :param factor:
        :param new_unit:
        :return:
        """

        # Multiply the frame with the conversion factor
        self.__imul__(factor)

        # Set the new unit
        self.unit = new_unit

    # -----------------------------------------------------------------

    def converted_by_factor(self, factor, new_unit):

        """
        This function ...
        :param factor:
        :param new_unit:
        :return:
        """

        # Create multiplicated frame
        new = self * factor

        # Set the new unit
        new.unit = new_unit

        # Return the new frame
        return new

    # -----------------------------------------------------------------

    @property
    def physical_type(self):

        """
        This function ...
        :return:
        """

        return self.unit.physical_type

    # -----------------------------------------------------------------

    @property
    def physical_base_type(self):

        """
        This function ...
        :return:
        """

        return self.unit.physical_base_type

    # -----------------------------------------------------------------

    @property
    def is_bolometric(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_bolometric

    # -----------------------------------------------------------------

    @property
    def corresponding_bolometric_unit(self):

        """
        This function ...
        :return:
        """

        return self.unit.corresponding_bolometric_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_bolometric_unit(self):

        """
        This function ...
        :return:
        """

        # Convert, return the factor
        return self.convert_to(self.corresponding_bolometric_unit)

    # -----------------------------------------------------------------

    def converted_to_corresponding_bolometric_unit(self):

        """
        This function ...
        :return:
        """

        return self.converted_to(self.corresponding_bolometric_unit)

    # -----------------------------------------------------------------

    @property
    def is_spectral_density(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_spectral_density

    # -----------------------------------------------------------------

    @property
    def is_wavelength_density(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_wavelength_density

    # -----------------------------------------------------------------

    @property
    def is_frequency_density(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_frequency_density

    # -----------------------------------------------------------------

    @property
    def is_neutral_density(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_neutral_density

    # -----------------------------------------------------------------

    @property
    def corresponding_wavelength_density_unit(self):
        return self.unit.corresponding_wavelength_density_unit

    # -----------------------------------------------------------------

    def get_corresponding_wavelength_density_unit(self, wavelength_unit=None):
        return self.unit.get_corresponding_wavelength_density_unit(wavelength_unit=wavelength_unit)

    # -----------------------------------------------------------------

    def convert_to_corresponding_wavelength_density_unit(self, distance=None, wavelength=None, wavelength_unit=None):

        """
        This function ...
        :param distance:
        :param wavelength:
        :param wavelength_unit:
        :return:
        """

        # Get unit
        if wavelength_unit is not None: unit = self.get_corresponding_wavelength_density_unit(wavelength_unit=wavelength_unit)
        else: unit = self.corresponding_wavelength_density_unit

        # Convert, return the factor
        return self.convert_to(unit, distance=distance, wavelength=wavelength)

    # -----------------------------------------------------------------

    def converted_to_corresponding_wavelength_density_unit(self, wavelength_unit=None):

        """
        This function ...
        :param wavelength_unit:
        :return:
        """

        return self.converted_to(self.corresponding_wavelength_density_unit)

    # -----------------------------------------------------------------

    @property
    def corresponding_frequency_density_unit(self):
        return self.unit.corresponding_frequency_density_unit

    # -----------------------------------------------------------------

    def get_corresponding_frequency_density_unit(self, frequency_unit=None):
        return self.unit.get_corresponding_frequency_density_unit(frequency_unit=frequency_unit)

    # -----------------------------------------------------------------

    def convert_to_corresponding_frequency_density_unit(self, wavelength=None, frequency_unit=None, distance=None):

        """
        This function ...
        :param wavelength:
        :param frequency_unit:
        :param distance:
        :return:
        """

        # Get unit
        if frequency_unit is not None: unit = self.get_corresponding_frequency_density_unit(frequency_unit=frequency_unit)
        else: unit = self.corresponding_frequency_density_unit

        # Convert, return the factor
        return self.convert_to(unit, wavelength=wavelength, distance=distance)

    # -----------------------------------------------------------------

    def converted_to_corresponding_frequency_density_unit(self):

        """
        This function ...
        :return:
        """

        return self.converted_to(self.corresponding_frequency_density_unit)

    # -----------------------------------------------------------------

    @property
    def corresponding_neutral_density_unit(self):
        return self.unit.corresponding_neutral_density_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_neutral_density_unit(self):

        """
        This function ...
        :return:
        """

        # Convert, return the factor
        return self.convert_to(self.corresponding_neutral_density_unit)

    # -----------------------------------------------------------------

    def converted_to_corresponding_neutral_density_unit(self):

        """
        This function ...
        :return:
        """

        return self.converted_to(self.corresponding_neutral_density_unit)

    # -----------------------------------------------------------------

    @property
    def is_brightness(self):
        return self.unit.is_brightness

    # -----------------------------------------------------------------

    @property
    def corresponding_brightness_unit(self):
        return self.unit.corresponding_brightness_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_brightness_unit(self):

        """
        This function ...
        :return:
        """

        return self.convert_to(self.corresponding_brightness_unit)

    # -----------------------------------------------------------------

    def converted_to_corresponding_brightness_unit(self):

        """
        This function ...
        :return:
        """

        return self.converted_to(self.corresponding_brightness_unit)

    # -----------------------------------------------------------------

    @property
    def corresponding_non_brightness_unit(self):
        return self.unit.corresponding_non_brightness_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_non_brightness_unit(self):

        """
        This function ...
        :return:
        """

        # Convert, return the factor
        return self.convert_to(self.corresponding_non_brightness_unit)

    # -----------------------------------------------------------------

    def converted_to_corresponding_non_brightness_unit(self):

        """
        This function ...
        :return:
        """

        return self.converted_to(self.corresponding_non_brightness_unit)

    # -----------------------------------------------------------------

    @property
    def is_surface_brightness(self):
        return self.unit.is_surface_brightness

    # -----------------------------------------------------------------

    @property
    def is_intrinsic_brightness(self):
        return self.unit.is_intrinsic_brightness

    # -----------------------------------------------------------------

    @property
    def is_per_angular_or_intrinsic_area(self):
        return self.unit.is_per_angular_or_intrinsic_area

    # -----------------------------------------------------------------

    @property
    def corresponding_angular_or_intrinsic_area_unit(self):
        return self.unit.corresponding_angular_or_intrinsic_area_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_angular_or_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        return self.convert_to(self.corresponding_angular_or_intrinsic_area_unit)

    # -----------------------------------------------------------------

    def converted_to_corresponding_angular_or_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        return self.converted_to(self.corresponding_angular_or_intrinsic_area_unit)

    # -----------------------------------------------------------------

    @property
    def corresponding_non_angular_or_intrinsic_area_unit(self):
        return self.unit.corresponding_non_angular_or_intrinsic_area_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_non_angular_or_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        return self.convert_to(self.corresponding_non_angular_or_intrinsic_area_unit)

    # -----------------------------------------------------------------

    def converted_to_corresponding_non_angular_or_intrinsic_area_unit(self):

        """
        This function ...
        :return:
        """

        return self.converted_to(self.corresponding_non_angular_or_intrinsic_area_unit)

    # -----------------------------------------------------------------

    @property
    def is_per_angular_area(self):
        return self.unit.is_per_angular_area

    # -----------------------------------------------------------------

    @property
    def is_per_intrinsic_area(self):
        return self.unit.is_per_intrinsic_area

    # -----------------------------------------------------------------

    @property
    def corresponding_angular_area_unit(self):
        return self.unit.corresponding_angular_area_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_angular_area_unit(self, distance=None):

        """
        This function ...
        :param distance: 
        :return: 
        """

        if distance is None: distance = self.distance
        # Convert, return the factor
        return self.convert_to(self.corresponding_angular_area_unit, distance=distance)

    # -----------------------------------------------------------------

    def converted_to_corresponding_angular_area_unit(self, distance=None):

        """
        This function ...
        :param distance:  
        :return: 
        """

        if distance is None: distance = self.distance
        return self.converted_to(self.corresponding_angular_area_unit, distance=distance)

    # -----------------------------------------------------------------

    @property
    def corresponding_intrinsic_area_unit(self):
        return self.unit.corresponding_intrinsic_area_unit

    # -----------------------------------------------------------------

    def convert_to_corresponding_intrinsic_area_unit(self, distance=None):

        """
        This function ...
        :param distance:
        :return:
        """

        if distance is None: distance = self.distance

        # Convert, return the factor
        return self.convert_to(self.corresponding_intrinsic_area_unit, distance=distance)

    # -----------------------------------------------------------------

    def converted_to_corresponding_intrinsic_area_unit(self, distance=None):

        """
        This function ...
        :param distance:
        :return:
        """

        if distance is None: distance = self.distance
        return self.converted_to(self.corresponding_intrinsic_area_unit, distance=distance)

    # -----------------------------------------------------------------

    def sum(self, add_unit=False, per_area="warning"):

        """
        This function ...
        :param add_unit:
        :param per_area:
        :return:
        """

        # Check
        unit = conversion_factor = None
        if self.has_unit and self.is_per_angular_or_intrinsic_area:

            # Set message
            if self.is_per_angular_area: message = "Unit is per angular area: the result of adding pixel values may not be useful before a conversion to a non-intensity or brightness unit"
            else: message = "Unit is per physical area: the result of adding pixel values may not be useful before a conversion to a non-brightness unit"

            # Warning or error
            if per_area == "warning": log.warning(message)
            elif per_area == "error": raise ValueError(message)
            elif per_area == "convert":
                unit = self.corresponding_non_angular_or_intrinsic_area_unit
                conversion_factor = self._get_conversion_factor(unit)
            else: raise ValueError("Invalid option for 'per_area'")

        # Calculate
        result = np.nansum(self.data)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def quadratic_sum(self, add_unit=False, per_area="warning"):

        """
        This function ...
        :param add_unit:
        :param per_area:
        :return:
        """

        # Check
        unit = conversion_factor = None
        if self.has_unit and self.is_per_angular_or_intrinsic_area:

            # Set message
            if self.is_per_angular_area: message = "Unit is per angular area: the result of adding pixel values may not be useful before a conversion to a non-intensity or brightness unit"
            else: message = "Unit is per physical area: the result of adding pixel values may not be useful before a conversion to a non-brightness unit"

            # Warning or error
            if per_area == "warning": log.warning(message)
            elif per_area == "error": raise ValueError(message)
            elif per_area == "convert":
                unit = self.corresponding_non_angular_or_intrinsic_area_unit
                conversion_factor = self._get_conversion_factor(unit)
            else: raise ValueError("Invalid option for 'per_area'")

        # Calculate
        result = np.sqrt(np.sum(self._data[self.nans.inverse()]**2))

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def mean(self, add_unit=False, not_per_area="warning"):

        """
        This function ...
        :param add_unit:
        :param not_per_area:
        :return:
        """

        # Check?
        unit = conversion_factor = None
        if self.is_photometric and not self.is_per_angular_or_intrinsic_area:

            # Set message
            message = "Unit is not per angular or intrinsic area: the result of taking the mean value may not be useful before a conversion to a pixelscale-independent unit (e.g. surface brightness)"

            # Warning or error
            if not_per_area == "warning": log.warning(message)
            elif not_per_area == "error": raise ValueError(message)
            elif not_per_area == "convert":
                unit = self.corresponding_angular_or_intrinsic_area_unit
                conversion_factor = self._get_conversion_factor(unit)
            else: raise ValueError("Invalid option for 'per_area'")

        # Calculate
        result = np.nanmean(self.data)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def average(self, add_unit=False, not_per_area="warning"):

        """
        This function ...
        :param add_unit:
        :param not_per_area:
        :return:
        """

        return self.mean(add_unit=add_unit, not_per_area=not_per_area)

    # -----------------------------------------------------------------

    def median(self, add_unit=False, not_per_area="warning"):

        """
        This function ...
        :param add_unit:
        :param not_per_area:
        :return:
        """

        # Check?
        unit = conversion_factor = None
        if self.is_photometric and not self.is_per_angular_or_intrinsic_area:

            # Set message
            message = "Unit is not per angular or intrinsic area: the result of taking the median value may not be useful before a conversion to a pixelscale-independent unit (e.g. surface brightness)"

            # Warning or error
            if not_per_area == "warning": log.warning(message)
            elif not_per_area == "error": raise ValueError(message)
            elif not_per_area == "convert":
                unit = self.corresponding_angular_or_intrinsic_area_unit
                conversion_factor = self._get_conversion_factor(unit)
            else: raise ValueError("Invalid option for 'per_area'")

        # Calculate
        result = np.nanmedian(self.data)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def stddev(self, add_unit=False, not_per_area="warning"):

        """
        This function ...
        :param add_unit:
        :param not_per_area:
        :return:
        """

        # Check?
        unit = conversion_factor = None
        if self.is_photometric and not self.is_per_angular_or_intrinsic_area:

            # Set message
            message = "Unit is not per angular or intrinsic area: the result of taking the standard deviation may not be useful before a conversion to a pixelscale-independent unit (e.g. surface brightness)"

            # Warning or error
            if not_per_area == "warning": log.warning(message)
            elif not_per_area == "error": raise ValueError(message)
            elif not_per_area == "convert":
                unit = self.corresponding_angular_or_intrinsic_area_unit
                conversion_factor = self._get_conversion_factor(unit)
            else: raise ValueError("Invalid option for 'per_area'")

        # Calculate
        result = np.nanstd(self.data)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def normalize(self, to=1.0, silent=False):

        """
        This function ...
        :param to:
        :param silent:
        :return:
        """

        # Calculate the sum of all the pixels
        sum = np.nansum(self)

        # Check whether the sum is nonnegative
        if sum < 0: raise RuntimeError("The sum of the frame is negative")

        # Check if the sum is not zero
        if sum == 0: raise AllZeroError("The frame cannot be normalized")

        # Calculate the conversion factor
        if hasattr(to, "unit"): # quantity
            factor = to.value / sum
            unit = to.unit
        else:
            factor = to / sum
            unit = None

        # Debugging
        if not silent: log.debug("Multiplying the frame with a factor of " + tostr(factor) + " to normalize to " + tostr(to) + " ...")

        # Multiply the frame with the conversion factor
        try: self.__imul__(factor)
        except TypeError:
            print(np.nanmax(self.data))
            from ..tools import plotting
            plotting.plot_frame(self)
            exit()

        # Set the unit
        self.unit = unit

        # Return the sum
        return sum

    # -----------------------------------------------------------------

    def normalized(self, to=1., return_sum=False, silent=False):
        
        """
        This function ...
        :param to:
        :param return_sum:
        :param silent:
        :return: 
        """

        new = self.copy()
        normalization_sum = new.normalize(to=to, silent=silent)
        if return_sum: return new, normalization_sum
        else: return new

    # -----------------------------------------------------------------

    def is_normalized(self, to=1., rtol=1.e-5, atol=1.e-8):

        """
        This function ...
        :param to:
        :param rtol
        :parma atol:
        :return:
        """

        # In some unit
        if hasattr(to, "unit"): # quantity

            total = self.sum() * self.unit
            return np.isclose(total.to(to.unit).value, to.value, atol=atol, rtol=rtol)

        # Numerically
        else: return np.isclose(self.sum(), to, atol=atol, rtol=rtol)

    # -----------------------------------------------------------------

    def convolved(self, *args, **kwargs):

        """
        This function ...
        :return:
        """

        new = self.copy()
        new.convolve(*args, **kwargs)
        return new

    # -----------------------------------------------------------------

    def convolve(self, kernel, allow_huge=True, fft=True, preserve_nans=True):

        """
        This function ...
        :param kernel:
        :param allow_huge:
        :param fft:
        :param preserve_nans:
        :return:
        """

        # Get the kernel FWHM and PSF filter
        kernel_fwhm = kernel.fwhm
        kernel_psf_filter = kernel.psf_filter

        # Skip the calculation for a constant frame
        if self.is_constant:
            self._fwhm = kernel_fwhm
            self._psf_filter = kernel_psf_filter
            return

        # Check whether the kernel is prepared
        if not kernel.prepared:
            log.warning("The convolution kernel is not prepared, creating a prepared copy ...")
            kernel = kernel.copy()
            kernel.prepare(self.pixelscale)

        # Check where the NaNs are at
        nans_mask = np.isnan(self._data)

        # Assert that the kernel is normalized
        if not kernel.normalized: raise RuntimeError("The kernel is not properly normalized: sum is " + repr(kernel.sum()) + " , difference from unity is " + repr(kernel.sum() - 1.0))

        # Get the current minimum and maximum of the frame
        min_value = self.min
        max_value = self.max

        # Debugging
        log.debug("The minimum and maximum value of the frame before convolution is " + tostr(min_value) + " and " + tostr(max_value))

        # Do the convolution on this frame
        log.debug("Convolving ...")
        if fft: new_data = convolve_fft(self._data, kernel.data, boundary='fill', nan_treatment="interpolate", normalize_kernel=False, allow_huge=allow_huge)
        else: new_data = convolve(self._data, kernel.data, boundary='fill', nan_treatment="interpolate", normalize_kernel=False)

        # Determine new min and max value
        new_min_value = np.nanmin(self.min)
        new_max_value = np.nanmax(self.max)

        # Debugging
        log.debug("The minimum and maximum value of the frame after convolution is " + tostr(new_min_value) + " and " + tostr(new_max_value))

        # Put back NaNs
        if preserve_nans: new_data[nans_mask] = nan_value

        # Don't mess up the scale
        # Set values lower than min to min value
        # and values above max to max value
        new_data[new_data < min_value] = min_value
        new_data[new_data > max_value] = max_value

        # Replace the data and FWHM
        self._data = new_data
        self._fwhm = kernel_fwhm
        self._psf_filter = kernel_psf_filter

    # -----------------------------------------------------------------

    def smoothed(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        new = self.copy()
        new.smooth(*args, **kwargs)
        return new

    # -----------------------------------------------------------------

    def smooth(self, factor, allow_huge=True, fft=True, preserve_nans=True, kernel_sigma_level=5.0):

        """
        This function ...
        :param factor:
        :param allow_huge:
        :param fft:
        :param preserve_nans:
        :param kernel_sigma_level:
        :return:
        """

        # Check whether the FWHM of the frame is defined
        if self.fwhm is None: raise ValueError("Cannot smooth if the FWHM of the frame is not defined")

        # Check whether the pixelscale is defined
        if self.pixelscale is None: raise ValueError("Cannot smooth if the pixelscale of the frame is not defined")

        # Determine the new FWHM
        new_fwhm = self.fwhm * factor

        # Debugging
        log.debug("The FWHM after smoothing will be " + tostr(new_fwhm))

        # Create convolution kernel
        from .kernel import ConvolutionKernel
        kernel = ConvolutionKernel.gaussian(new_fwhm, self.pixelscale, sigma_level=kernel_sigma_level)

        # Get the original PSF filter and FWHM
        original_fwhm = self.fwhm
        original_psf_filter = self.psf_filter

        # Convolve
        self.convolve(kernel, allow_huge=allow_huge, fft=fft, preserve_nans=preserve_nans)

        # Set the FWHM and PSF filter back to the original
        self.fwhm = original_fwhm
        self.psf_filter = original_psf_filter

        # BUT SET THE SMOOTHING FACTOR
        self.smoothing_factor *= factor

    # -----------------------------------------------------------------

    def interpolated(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        new = self.copy()
        new.interpolate(*args, **kwargs)
        return new

    # -----------------------------------------------------------------

    def interpolate(self, region_or_mask, sigma=None, max_iterations=10, plot=False, not_converge="keep",
                    min_max_in=None, smoothing_factor=None, replace_nans=None):

        """
        Thisfunction ...
        :param region_or_mask:
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :param smoothing_factor:
        :param interpolate_nans:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Get a mask of the original NaN pixels
        original_nans = self.nans

        # Set originally NaN pixels to something else? zero?
        if replace_nans: self[original_nans] = replace_nans

        # Set nans at masked pixels
        original_values = self[mask]
        self[mask] = nan_value

        # Interpolate the nans
        try: self.interpolate_nans(sigma=sigma, max_iterations=max_iterations, plot=plot, not_converge=not_converge,
                                   min_max_in=min_max_in, smoothing_factor=smoothing_factor)
        except RuntimeError as e:

            # Reset the original values (e.g. infs)
            self[mask] = original_values

            # Set original Nans Back to Nan
            self[original_nans] = nan_value

            # Reraise the error
            raise e

        # Set original NaNs back to NaN
        self[original_nans] = nan_value

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def interpolated_nans(self, **kwargs):

        """
        Thisfunction ...
        :param kwargs:
        :return:
        """

        new = self.copy()
        new.interpolate_nans(**kwargs)
        return new

    # -----------------------------------------------------------------

    def interpolate_nans_if_below(self, threshold=0.7, sigma=None, max_iterations=None, min_max_in=None, smoothing_factor=None):

        """
        This function ...
        :param threshold:
        :param sigma:
        :param max_iterations:
        :param min_max_in:
        :param smoothing_factor:
        :return:
        """

        # Get the number of relative nnans
        relnans = self.relative_nnans

        # UNDER THRESHOLD
        if relnans < threshold:

            # Debugging
            log.debug("The relative number of NaN values in the frame is " + str(relnans * 100) + "%")

            # Interpolate, returning the nans
            nans = self.interpolate_nans(sigma=sigma, max_iterations=max_iterations, min_max_in=min_max_in, smoothing_factor=smoothing_factor)

        # ABOVE THRESHOLD
        else:

            # Give warning
            log.warning("The number of NaN values in the frame is very high (" + str(relnans * 100) + "%)")

            # Set nans to None
            nans = None

        # Return the nans, None if nothing is done
        return nans

    # -----------------------------------------------------------------

    def interpolate_negatives_if_below(self, threshold=0.7, sigma=None, max_iterations=None, min_max_in=None):

        """
        This function ...
        :param threshold:
        :param sigma:
        :param max_iterations:
        :param min_max_in:
        :return:
        """

        # Get the number of relative negatives
        relnegatives = self.relative_nnegatives

        # UNDER THRESHOLD
        if relnegatives < threshold:

            # Debugging
            log.debug("The relative number of negative values in the frame is " + str(relnegatives * 100) + "%")

            # Interpolate, returning the negatives
            negatives = self.interpolate_negatives(sigma=sigma, max_iterations=max_iterations, min_max_in=min_max_in)

        # ABOVE THRESHOLD
        else:

            # Give warning
            log.warning("The number of negative values in the frame is very high (" + str(relnegatives * 100) + "%)")

            # Set negatives to None
            negatives = None

        # Return the negatives, None if nothing is done
        return negatives

    # -----------------------------------------------------------------

    def interpolate_nans(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None,
                         smoothing_factor=None, error_on_max=True):

        """
        This function ...
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge: what to do with NaN values when the number of NaN pixels does not converge to zero?
        #                     -> "error', or "keep"
        :param min_max_in:
        :param smoothing_factor:
        :param error_on_max:
        :return:
        """

        # Determine sigma
        if sigma is None:

            # Check whether we have the necessary information
            if self.fwhm is None: raise ValueError("FWHM of the frame should be defined or sigma should be passed")
            if self.pixelscale is None: raise ValueError("Pixelscale of the frame is not defined")

            # Get the sigma in pixels
            sigma = self.sigma_pix

            # Smoothing factor
            if smoothing_factor is not None:
                if smoothing_factor < 1: raise ValueError("Smoothing factor cannot be smaller than one")
                log.debug("Original sigma of the frame resolution is " + tostr(sigma) + " pixels")
                log.debug("Interpolated regions will be smoother by a factor of " + str(smoothing_factor))
                sigma = sigma * smoothing_factor

        # Smoothing factor is passed but also sigma
        elif smoothing_factor is not None: raise ValueError("Smoothing factor cannot be passed when sigma is passed: multiply the specified sigma with the desired smoothing factor")

        # Debugging
        log.debug("Creating a kernel with a sigma of " + tostr(sigma) + " pixels ...")

        # We smooth with a Gaussian kernel with stddev passed by the user
        # Create the kernel
        kernel = Gaussian2DKernel(stddev=sigma)

        # Interpolate
        return self.interpolate_nans_with_kernel(kernel, plot=plot, max_iterations=max_iterations,
                                                 not_converge=not_converge, min_max_in=min_max_in, error_on_max=error_on_max)

    # -----------------------------------------------------------------

    def interpolate_nans_with_kernel(self, kernel, plot=False, max_iterations=10, not_converge="keep", min_max_in=None,
                                     error_on_max=True):

        """
        This function ...
        :param kernel:
        :param plot:
        :param max_iterations:
        :param not_converge:
        :param min_max_in:
        :param error_on_max:
        :return:
        """

        from ..tools import plotting

        # Get the current minimum and maximum of the frame
        if min_max_in is not None:
            min_value = self.min_in(min_max_in)
            max_value = self.max_in(min_max_in)
        else:
            min_value = self.min
            max_value = self.max

        # Debugging
        log.debug("The minimum and maximum value of the frame before interpolation is " + tostr(min_value) + " and " + tostr(max_value))

        # Debugging
        log.debug("Interpolation iteration 1 ...")

        # Generate the interpolated result
        result = interpolate_replace_nans(self.data, kernel)
        niterations = 1

        # Plot
        if plot: plotting.plot_box(result, title="Result after iteration " + str(niterations))
        if plot: plotting.plot_mask(np.isnan(result), title="NaNs after iteration " + str(niterations))

        # Get the current number of nans
        previous_nnans = None # don't get it for performance
        nnans = np.sum(np.isnan(result))

        # Debugging
        log.debug("The number of NaN values after iteration 1 is " + str(nnans))

        # Are there still NaNs?
        while nnans > 0:

            # Check number of iterations
            if max_iterations is not None and niterations == max_iterations:
                if error_on_max: raise RuntimeError("The maximum number of iterations has been reached without success")
                else:
                    log.warning("The maximum number of iterations has been reached without success (" + str(nnans) + " remaining NaNs)")
                    break # break the loop

            # Debugging
            log.debug("Interpolation iteration " + str(niterations+1) + " ...")

            # Perform next interpolation
            result = interpolate_replace_nans(result, kernel)

            # Increment the niterations counter
            niterations += 1

            # Get the current number of nans
            previous_nnans = nnans
            nnans = np.sum(np.isnan(result))

            # Not converging
            if nnans >= previous_nnans:
                if not_converge == "keep":
                    log.warning("The number of NaNs could not converge to zero: " + str(nnans) + " NaN values will remain")
                    break # break the loop
                elif not_converge == "error": raise RuntimeError("The number of NaNs is not converging to zero (nnans = " + str(nnans) + ", previous nnans = " + str(previous_nnans) + ")")
                else: raise ValueError("Invalid option for 'not_converge'")

            # Debugging
            log.debug("The number of NaN values after iteration " + str(niterations) + " is " + str(nnans))

            # Plot
            if plot: plotting.plot_box(result, title="Result after iteration " + str(niterations))
            if plot: plotting.plot_mask(np.isnan(result), title="NaNs after iteration " + str(niterations))

        # Determine new min and max value
        if min_max_in is not None:
            mask = self.get_mask(min_max_in)
            new_min_value = np.nanmin(result[mask])
            new_max_value = np.nanmax(result[mask])
        else:
            new_min_value = np.nanmin(result)
            new_max_value = np.nanmax(result)

        # Debugging
        log.debug("The minimum and maximum value of the frame after interpolation is " + tostr(new_min_value) + " and " + tostr(new_max_value))

        # Don't mess up the scale
        # Set values lower than min to min value
        # and values above max to max value
        result[result < min_value] = min_value
        result[result > max_value] = max_value

        # Get the mask of NaNs
        original_nans = self.nans

        # Replace the data
        self._data = result

        # Return the original nans
        return original_nans

    # -----------------------------------------------------------------

    def interpolated_infs(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        return self.interpolated(self.infs, **kwargs)

    # -----------------------------------------------------------------

    def interpolate_infs(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None):

        """
        This function ...
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :return:
        """

        return self.interpolate(self.infs, sigma=sigma, max_iterations=max_iterations, plot=plot, not_converge=not_converge, min_max_in=min_max_in)

    # -----------------------------------------------------------------

    def interpolated_negatives(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        return self.interpolated(self.negatives, **kwargs)

    # -----------------------------------------------------------------

    def interpolate_negatives(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None):

        """
        Thisf unction ...
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :return:
        """

        return self.interpolate(self.negatives, sigma=sigma, max_iterations=max_iterations, plot=plot, not_converge=not_converge, min_max_in=min_max_in)

    # -----------------------------------------------------------------

    def interpolated_zeroes(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        return self.interpolated(self.zeroes, **kwargs)

    # -----------------------------------------------------------------

    def interpolate_zeroes(self, sigma=None, max_iterations=10, plot=False, not_converge="keep", min_max_in=None):

        """
        This function ...
        :param sigma:
        :param max_iterations:
        :param plot:
        :param not_converge:
        :param min_max_in:
        :return:
        """

        return self.interpolate(self.zeroes, sigma=sigma, max_iterations=max_iterations, plot=plot, not_converge=not_converge, min_max_in=min_max_in)

    # -----------------------------------------------------------------

    def rebinned(self, reference_wcs, exact=False, parallel=True, convert=None):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :param convert:
        :return:
        """

        new = self.copy()
        new.rebin(reference_wcs, exact=exact, parallel=parallel, convert=convert)
        return new

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=False, parallel=True, convert=None):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :param convert:
        :return:
        """

        # Check if the reference WCS is not None
        if reference_wcs is None: raise ValueError("Reference coordinate system has to be specified")
        #if not isinstance(reference_wcs, CoordinateSystem): ... # extra check?

        # Check whether the frame has a WCS
        if not self.has_wcs: raise RuntimeError("Cannot rebin a frame without coordinate system")

        # Check whether the WCS is the same
        if self.wcs == reference_wcs: return Frame.ones_like(self)

        # Check the unit
        original_unit = self._test_angular_or_intrinsic_area("rebinning", convert=convert)

        # Calculate rebinned data and footprint of the original image
        if exact: new_data, footprint = reproject_exact((self._data, self.wcs), reference_wcs, shape_out=reference_wcs.shape, parallel=parallel)
        else: new_data, footprint = reproject_interp((self._data, self.wcs), reference_wcs, shape_out=reference_wcs.shape)

        # Replace the data and WCS
        self._data = new_data
        self._wcs = reference_wcs.copy()

        # Convert back?
        if original_unit is not None: self.convert_to(original_unit)

        # Return the footprint
        return Frame(footprint, wcs=reference_wcs.copy())

    # -----------------------------------------------------------------

    def _test_angular_or_intrinsic_area(self, task, convert=None):

        """
        This function ...
        :param task:
        :param convert:
        :return:
        """

        original_unit = None

        # The unit of the frame is not defined
        if self.unit is None: log.warning("The unit of this frame is not defined. Be aware of the fact that " + task + " a frame not in brightness units gives an incorrect result")
        elif not self.is_photometric: log.warning("The frame does not have a photometric unit. Cannot test whether the unit is per angular or intrinsic area.")
        elif not self.is_per_angular_or_intrinsic_area:

            # Pixelscale is not defined
            if not self.has_pixelscale: raise RuntimeError("The pixelscale of the frame is not defined")

            # Whether or not to convert is not specified
            if convert is None: raise RuntimeError("The frame is not defined in units per angular or physical area, but in units of " + str(self.unit) + " [" + str(self.unit.physical_type) + "]. Convert to a related intensity or brightness unit prior to " + task)

            # Convert to angular or intrinsic area unit
            elif convert is True:

                original_unit = self.unit
                self.convert_to_corresponding_angular_or_intrinsic_area_unit()

            # Don't convert:
            elif convert is False:

                # Give warning
                log.warning("Not converting the unit prior to " + task + " will provide incorrect pixel values")

                # Return None: no conversion needed afterwards
                return None

            # Invalid
            else: raise ValueError("Invalid value for 'convert'")

        # Return the original unit
        return original_unit

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

        # Pixel rectangle
        if isinstance(region, PixelRectangleRegion):
            if factor != 1: region = region * factor
            return self.crop(region.x_min_pixel, region.x_max_pixel, region.y_min_pixel, region.y_max_pixel, out_of_bounds=out_of_bounds)

        # Sky rectangle: to pixel rectangle
        elif isinstance(region, SkyRectangleRegion): return self.crop_to(region.to_pixel(self.wcs), factor=factor, out_of_bounds=out_of_bounds)

        # Other kind of shape
        else: return self.crop_to(region.bounding_box, factor=factor, out_of_bounds=out_of_bounds)

    # -----------------------------------------------------------------

    def soften_edges(self, region, factor_range, invert=False):

        """
        This function ...
        :param region:
        :param factor_range:
        :param invert:
        :return:
        """

        # Create alpha mask
        alpha = AlphaMask.from_ellipse(region, self.shape, factor_range, wcs=self.wcs)

        # Apply the alpha mask
        self.apply_alpha_mask(alpha, invert=invert)

        # Return the alpha mask
        return alpha

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
        return PixelCoordinate(self.x_center, self.y_center)

    # -----------------------------------------------------------------

    @property
    def center_sky(self):
        if not self.has_wcs: raise ValueError("Cannot determine the center as sky coordinate: coordinate system not defined")
        else: return self.center.to_sky(self.wcs)

    # -----------------------------------------------------------------

    @property
    def pixel_center(self):
        return Pixel.for_coordinate(self.center)

    # -----------------------------------------------------------------

    @property
    def reference_pixel(self):
        return self.wcs.reference_pixel if self.has_wcs else None

    # -----------------------------------------------------------------

    @property
    def reference_coordinate(self):
        return self.wcs.reference_coordinate if self.has_wcs else None

    # -----------------------------------------------------------------

    @property
    def center_value(self):
        return self.data[self.pixel_center.y, self.pixel_center.x]

    # -----------------------------------------------------------------

    def padded(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        new = self.copy()
        new.pad(nx, ny)
        return new

    # -----------------------------------------------------------------

    def pad(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        if nx == 0 and ny == 0: return

        new_data = np.pad(self._data, ((ny,0), (nx,0)), 'constant')

        if self.wcs is not None:

            new_wcs = copy.deepcopy(self.wcs)

            new_wcs.wcs.crpix[0] += nx
            new_wcs.wcs.crpix[1] += ny

            # Change the number of pixels
            new_wcs.naxis1 = new_data.shape[1]
            new_wcs.naxis2 = new_data.shape[0]
            new_wcs._naxis1 = new_wcs.naxis1
            new_wcs._naxis2 = new_wcs.naxis2

        else: new_wcs = None

        # Set the new data and the new WCS
        self._data = new_data
        self._wcs = new_wcs

    # -----------------------------------------------------------------

    def unpad(self, nx=0, ny=0):

        """
        This function ...
        :param nx:
        :param ny:
        :return:
        """

        if nx == 0 and ny == 0: return

        # Slice
        new_data = self._data[ny:, nx:]

        if self.wcs is not None:

            new_wcs = self.wcs.copy()

            new_wcs.wcs.crpix[0] -= nx
            new_wcs.wcs.crpix[1] -= ny

            # Change the number of pixels
            new_wcs.naxis1 = new_data.shape[1]
            new_wcs.naxis2 = new_data.shape[0]
            new_wcs._naxis1 = new_wcs.naxis1
            new_wcs._naxis2 = new_wcs.naxis2

        else: new_wcs = None

        # Set the new data and the new WCS
        self._data = new_data
        self._wcs = new_wcs

    # -----------------------------------------------------------------

    def dilate_nans(self, radius=5, niterations=1):

        """
        This function ...
        :param radius:
        :param niterations:
        :return:
        """

        # Get original nans
        nans = self.nans
        original_nans = nans.copy()

        # Dilate
        #nans.dilate_rc(rank, connectivity=connectivity, iterations=iterations)
        nans.disk_dilate(radius=radius, niterations=niterations)

        # Apply
        self.apply_mask_nans(nans)

        # Return the original nans
        return original_nans

    # -----------------------------------------------------------------

    def downsampled(self, factor, order=3, dilate_nans=True, dilate_infs=True, convert=None):

        """
        This function ...
        :param factor:
        :param order:
        :param dilate_nans:
        :param dilate_infs:
        :param convert:
        :return:
        """

        new = self.copy()
        new.downsample(factor, order=order, dilate_nans=dilate_nans, dilate_infs=dilate_infs, convert=convert)
        return new

    # -----------------------------------------------------------------

    def downsample(self, factor, order=3, dilate_nans=True, dilate_infs=True, convert=None):

        """
        This function ...
        :param factor:
        :param order:
        :param dilate_nans
        :param dilate_infs:
        :param convert:
        :return:
        """

        # Check if factor is 1: nothing needs to be done
        if factor == 1: return

        # Get mask of nans and infs
        if self.has_nans: nans = self.nans
        else: nans = None
        if self.has_infs: infs = self.infs
        else: infs = None

        # Put zeroes where nans and infs were
        if nans is not None: self.apply_mask(nans)
        if infs is not None: self.apply_mask(infs)

        # Check the unit
        original_unit = self._test_angular_or_intrinsic_area("downsampling", convert=convert)

        # Set mode
        #mode = "constant" # default
        #mode = "nearest"
        #mode = "reflect"
        mode = "wrap"

        # Calculate the downsampled array
        zoom_factor = 1./factor
        new_data = ndimage.interpolation.zoom(self._data, zoom=zoom_factor, order=order, mode=mode)
        new_xsize = new_data.shape[1]
        new_ysize = new_data.shape[0]

        # Does the frame have a WCS?
        new_wcs = None
        if self.wcs is not None: new_wcs, nans, infs = self._scale_wcs_and_masks(new_xsize, new_ysize, nans, infs, dilate_nans=dilate_nans, dilate_infs=dilate_infs)

        # No WCS
        else: nans, infs = self._downsample_masks(factor, new_xsize, new_ysize, nans, infs, dilate_nans=dilate_nans, dilate_infs=dilate_infs)

        # Apply masks
        if nans is not None: new_data[nans.data] = nan_value
        if infs is not None: new_data[infs.data] = inf_value

        # Set the new data and wcs
        self._data = new_data
        self._wcs = new_wcs

        # Convert back?
        if original_unit is not None: self.convert_to(original_unit)

        # Return the new coordinate system
        return new_wcs

    # -----------------------------------------------------------------

    def _scale_wcs_and_masks(self, new_xsize, new_ysize, nans, infs, dilate_nans=True, dilate_infs=True):

        """
        This function ...
        :param new_xsize:
        :param new_ysize:
        :param nans:
        :param infs:
        :param dilate_nans:
        :param dilate_infs:
        :return:
        """

        # Create the new WCS
        new_wcs = self.wcs.scaled(new_xsize, new_ysize)

        # Rebin nans
        if nans is not None: nans.rebin(new_wcs, dilate=dilate_nans)

        # Rebin infs
        if infs is not None: infs.rebin(new_wcs, dilate=dilate_infs)

        # Return the new WCS and masks
        return new_wcs, nans, infs

    # -----------------------------------------------------------------

    def _downsample_masks(self, factor, new_xsize, new_ysize, nans, infs, dilate_nans=True, dilate_infs=True):

        """
        This function ...
        :param factor:
        :param new_xsize:
        :param new_ysize:
        :param dilate_nans:
        :param dilate_infs:
        :return:
        """

        # Downsample nans
        if nans is not None:

            # Downsample with dilation
            new_nans = nans.downsampled(factor, dilate=dilate_nans)

            # Check the shapes
            if new_nans.shape[1] != new_xsize or new_nans.shape[0] != new_ysize:

                # Give warning
                log.warning("Could not downsample the mask of NaN values: all NaNs have been replaced by zero")

                # Set nans to None
                new_nans = None

        # No nans mask
        else: new_nans = None

        # Downsample infs
        if infs is not None:

            # Downsample with dilation
            new_infs = infs.downsampled(factor, dilate=dilate_infs)

            # Check the shapes
            if new_infs.shape[1] != new_xsize or new_infs.shape[0] != new_ysize:

                # Give warning
                log.warning("Could not downsample the mask of infinite values: all infinities have been replaced by zero")

                # Set infs to None
                new_infs = None

        # No infs mask
        else: new_infs = None

        # Return the masks
        return new_nans, new_infs

    # -----------------------------------------------------------------

    def downsample_to_ratio(self, ratio, order=3, integer=True, above_or_below=None, even_or_odd=None, dilate_nans=True,
                            dilate_infs=True):

        """
        This function ...
        :param ratio:
        :param order:
        :param integer:
        :param above_or_below:
        :param even_or_odd:
        :param dilate_nans:
        :param dilate_infs:
        :return:
        """

        from ...core.tools import numbers

        # Integer downsampling factors
        if integer:

            # Not necessarily above or below
            if above_or_below is None:

                if even_or_odd is None: downsample_factor = numbers.nearest_integer(ratio)
                elif even_or_odd == "even": downsample_factor = numbers.nearest_even_integer(ratio)
                elif even_or_odd == "odd": downsample_factor = numbers.nearest_odd_integer(ratio)
                else: raise ValueError("Invalid value for 'even_or_odd'")

            # Always above
            elif above_or_below == "above":

                if even_or_odd is None: downsample_factor = numbers.nearest_integer_below(ratio, below=ratio)
                elif even_or_odd == "even": downsample_factor = numbers.nearest_even_integer_below(ratio, below=ratio)
                elif even_or_odd == "odd": downsample_factor = numbers.nearest_odd_integer_below(ratio, below=ratio)
                else: raise ValueError("Invalid value for 'even_or_odd'")

            # Always below
            elif above_or_below == "below":

                if even_or_odd is None: downsample_factor = numbers.nearest_integer_above(ratio, above=ratio)
                elif even_or_odd == "even": downsample_factor = numbers.nearest_even_integer_above(ratio, above=ratio)
                elif even_or_odd == "odd": downsample_factor = numbers.nearest_odd_integer_above(ratio, above=ratio)
                else: raise ValueError("Invalid value for 'even_or_odd'")

            # Invalid
            else: raise ValueError("Invalid value for 'above_or_below'")

        # Not necessarily integer
        else: downsample_factor = ratio

        # Debugging
        log.debug("Downsampling with a factor of " + str(downsample_factor) + " ...")

        # Downsample
        self.downsample(downsample_factor, order=order, dilate_nans=dilate_nans, dilate_infs=dilate_infs)

        # Return the downsample factor
        return downsample_factor

    # -----------------------------------------------------------------

    def downsample_to_npixels(self, npixels, order=3, integer=True, above_or_below=None, even_or_odd=None,
                              dilate_nans=True, dilate_infs=True):

        """
        This function ...
        :param npixels:
        :param order:
        :param integer:
        :param above_or_below:
        :param even_or_odd:
        :param dilate_nans:
        :param dilate_infs:
        :return:
        """

        # Check
        if self.xsize <= npixels and self.ysize <= npixels:
            log.warning("Downsampling is not necessary: xsize = " + str(self.xsize) + " and ysize = " + str(self.ysize) + " below " + str(npixels))
            return 1

        # Determine the ratio
        ratio = max(self.xsize, self.ysize) / float(npixels) # scale so that the longest axis falls below the target number of pixels

        # Debugging
        log.debug("The ratio between the number of pixels along longest axis and the maximum number of pixels is " + str(ratio))

        # Downsample
        return self.downsample_to_ratio(ratio, order=order, integer=integer, above_or_below=above_or_below,
                                        even_or_odd=even_or_odd, dilate_nans=dilate_nans, dilate_infs=dilate_infs)

    # -----------------------------------------------------------------

    def downsample_to_pixelscale(self, pixelscale, order=3, integer=True, above_or_below=None, even_or_odd=None,
                                 dilate_nans=True, dilate_infs=True):

        """
        This function ...
        :param pixelscale:
        :param order:
        :param integer:
        :param above_or_below: None, 'above', or 'below'
        :param even_or_odd:
        :param dilate_nans:
        :param dilate_infs:
        :return:
        """

        # Check
        if not self.has_pixelscale: raise ValueError("Pixelscale of the frame is undefined")

        # Check
        if pixelscale < self.pixelscale.average: raise ValueError("Pixelscale of the frame is greater than the target pixelscale")

        # Determine the downsample factor
        ratio = (pixelscale / self.pixelscale.average).to("").value

        # Debugging
        log.debug("The ratio between the target pixelscale and the original pixelscale is " + str(ratio))

        # Downsample
        return self.downsample_to_ratio(ratio, order=order, integer=integer, above_or_below=above_or_below,
                                        even_or_odd=even_or_odd, dilate_nans=dilate_nans, dilate_infs=dilate_infs)

    # -----------------------------------------------------------------

    def upsampled(self, factor, order=3):

        """
        This function ...
        :param factor:
        :param order:
        :return:
        """

        new = self.copy()
        new.upsample(factor, order=order)
        return new

    # -----------------------------------------------------------------

    def upsample_integers(self, factor): # For segmentation map class??

        """
        This function ...
        :param factor:
        :return:
        """

        # Check whether the upsampling factor is an integer or not
        if int(factor) == factor:

            new_data = ndimage.zoom(self._data, factor, order=0)

            new_xsize = new_data.shape[1]
            new_ysize = new_data.shape[0]

            relative_center = Position(self.reference_pixel.x / self.xsize, self.reference_pixel.y / self.ysize)

            new_center = Position(relative_center.x * new_xsize, relative_center.y * new_ysize)

            new_wcs = copy.deepcopy(self.wcs)
            # Change the center pixel position
            new_wcs.wcs.crpix[0] = new_center.x
            new_wcs.wcs.crpix[1] = new_center.y

            # Change the number of pixels
            new_wcs.naxis1 = new_xsize
            new_wcs.naxis2 = new_ysize
            new_wcs._naxis1 = new_wcs.naxis1
            new_wcs._naxis2 = new_wcs.naxis2

            # Change the pixel scale
            new_wcs.wcs.cdelt[0] *= float(self.xsize) / float(new_xsize)
            new_wcs.wcs.cdelt[1] *= float(self.ysize) / float(new_ysize)

            # return Frame(data, wcs=new_wcs, name=self.name, description=self.description, unit=self.unit, zero_point=self.zero_point, filter=self.filter, sky_subtracted=self.sky_subtracted, fwhm=self.fwhm)

            # Set the new data and wcs
            self._data = new_data
            self._wcs = new_wcs

        # Upsampling factor is not an integer
        else:

            old = self.copy()

            self.downsample(1. / factor)

            # print("Checking indices ...")
            indices = np.unique(old._data)

            # print("indices:", indices)

            # Loop over the indices
            for index in list(indices):
                # print(index)

                index = int(index)

                where = old._data == index

                # Calculate the downsampled array
                data = ndimage.interpolation.zoom(where.astype(float), zoom=factor)
                upsampled_where = data > 0.5

                self[upsampled_where] = index

    # -----------------------------------------------------------------

    def upsample(self, factor, order=3):

        """
        This function ...
        :param factor:
        :param order:
        :return:
        """

        # Just do inverse of downsample
        self.downsample(1./factor, order=order)

    # -----------------------------------------------------------------

    def fill(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self._data.fill(value)

    # -----------------------------------------------------------------

    def fill_nans(self):

        """
        This function ...
        :return:
        """

        self.fill(nan_value)

    # -----------------------------------------------------------------

    def fill_infs(self):

        """
        This function ...
        :return:
        """

        self.fill(inf_value)

    # -----------------------------------------------------------------

    def fill_zeroes(self):

        """
        This function ...
        :return:
        """

        self.fill(zero_value)

    # -----------------------------------------------------------------

    def fill_in(self, region_or_mask, value):

        """
        This function ...
        :param region_or_mask:
        :param value:
        :return:
        """

        # Get mask
        mask = self.get_mask(region_or_mask)

        # Fill
        self._data[mask] = value

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def fill_nans_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        return self.fill_in(region_or_mask, nan_value)

    # -----------------------------------------------------------------

    def fill_infs_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        return self.fill_in(region_or_mask, inf_value)

    # -----------------------------------------------------------------

    def fill_zeroes_in(self, region_or_mask):

        """
        This function ...
        :param region_or_mask:
        :return:
        """

        return self.fill_in(region_or_mask, zero_value)

    # -----------------------------------------------------------------

    def rotate(self, angle):

        """
        This function ...
        :param angle:
        :return:
        """

        data = ndimage.interpolation.rotate(self.data, angle.to("deg").value, reshape=False, order=1, mode='constant', cval=float('nan'))

        # Replace data
        self._data = data

        # Rotate wcs
        if self.wcs is not None:

            rotated_wcs = self.wcs.deepcopy()
            #rotated_wcs.rotate_wcs(angle)  # STILL UNTESTED
            rotate_wcs(rotated_wcs, angle)

            # Replace wcs
            self.wcs = rotated_wcs

        # Return mask of padded pixels
        return self.nans

    # -----------------------------------------------------------------

    def rotation_mask(self, angle):

        """
        This function ...
        :param angle:
        :return:
        """

        data = ndimage.interpolation.rotate(self.data, angle.to("deg").value, reshape=False, order=1, mode='constant', cval=float('nan'))
        return Mask(np.isnan(data), wcs=self.wcs.copy() if self.wcs is not None else None, pixelscale=self.pixelscale)

    # -----------------------------------------------------------------

    def rotated(self, angle):

        """
        This function ...
        :param angle:
        :return:
        """

        # Calculate the rotated array
        #frame[np.isnan(frame)] = 0.0
        data = ndimage.interpolation.rotate(self.data, angle.to("deg").value, reshape=False, order=1, mode='constant', cval=float('nan'))
        #new_frame = misc.imrotate(frame, angle, interp="bilinear")

        # Convert the wcs to header
        #header = self.wcs.to_header()

        # Rotate the header (Sebastien's script)
        #from ..tools import rotation
        #rotated_header = rotation.rotate_header(header, angle)

        # Create the new WCS
        #rotated_wcs = CoordinateSystem(rotated_header)

        rotated_wcs = self.wcs.deepcopy()
        #rotated_wcs.rotate_wcs(angle) # STILL UNTESTED
        rotate_wcs(rotated_wcs, angle)

        # Return the rotated frame
        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        return Frame(data,
                     wcs=rotated_wcs,
                     name=self.name,
                     description=self.description,
                     unit=self.unit,
                     zero_point=self.zero_point,
                     filter=self.filter,
                     sky_subtracted=self.sky_subtracted,
                     fwhm=self.fwhm)

    # -----------------------------------------------------------------

    def shifted(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        # TODO: change the WCS !!!

        # Transform the data
        data = ndimage.interpolation.shift(self.data, (extent.y, extent.x))

        # Return the shifted frame
        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        return Frame(data,
                     wcs=None,
                     name=self.name,
                     description=self.description,
                     unit=self.unit,
                     zero_point=self.zero_point,
                     filter=self.filter,
                     sky_subtracted=self.sky_subtracted,
                     fwhm=self.fwhm)

    # -----------------------------------------------------------------

    def centered_around(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        center = Position(x=0.5*self.xsize, y=0.5*self.ysize)
        shift = position - center

        # Return the shifted frame
        return self.shifted(shift)

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):
        return self.wcs.bounding_box

    # -----------------------------------------------------------------

    @property
    def corners(self):

        """
        This function ...
        :return:
        """

        # Get coordinate values
        coordinate_values = self.wcs.calc_footprint(undistort=True)

        # Initialize a list to contain the coordinates of the corners
        corners = []

        for ra_deg, dec_deg in coordinate_values:

            # Create sky coordinate
            coordinate = SkyCoordinate(ra=ra_deg, dec=dec_deg, unit="deg", frame="fk5")

            # Add the coordinate of this corner to the list
            corners.append(coordinate)

        # Return the list of coordinates
        return corners

    # -----------------------------------------------------------------

    @property
    def coordinate_range(self):
        return self.wcs.coordinate_range

    # -----------------------------------------------------------------

    @property
    def coordinate_box(self):

        """
        This function ...
        :return:
        """

        # Get coordinate range
        center, ra_span, dec_span = self.coordinate_range

        #ra = center.ra.to(unit).value
        #dec = center.dec.to(unit).value

        #ra_span = ra_span.to(unit).value
        #dec_span = dec_span.to(unit).value

        # Create rectangle
        center = SkyCoordinate(center.ra, center.dec)
        radius = SkyStretch(0.5 * ra_span, 0.5 * dec_span)
        box = SkyRectangleRegion(center, radius)

        # Return the box
        return box

    # -----------------------------------------------------------------

    def contains(self, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        pixel = coordinate.to_pixel(self.wcs)
        return 0.0 <= pixel.x < self.xsize and 0.0 <= pixel.y < self.ysize

    # -----------------------------------------------------------------

    def replace_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        nans = self.nans
        self._data[nans] = value
        return nans

    # -----------------------------------------------------------------

    def replaced_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new = self.copy()
        new.replace_nans(value)
        return new

    # -----------------------------------------------------------------

    def replace_nans_by_zeroes(self):

        """
        This function ...
        :return:
        """

        return self.replace_nans(zero_value)

    # -----------------------------------------------------------------

    def replace_nans_by_infs(self):

        """
        This function ...
        :return:
        """

        return self.replace_nans(inf_value)

    # -----------------------------------------------------------------

    def replace_infs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        infs = self.infs
        self._data[infs] = value
        return infs

    # -----------------------------------------------------------------

    def replace_infs_by_zeroes(self):

        """
        This function ...
        :return:
        """

        return self.replace_infs(zero_value)

    # -----------------------------------------------------------------

    def replace_infs_by_nans(self):

        """
        This function ...
        :return:
        """

        return self.replace_infs(nan_value)

    # -----------------------------------------------------------------

    def replace_zeroes(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        zeroes = self.zeroes
        self._data[zeroes] = value
        return zeroes

    # -----------------------------------------------------------------

    def replace_zeroes_by_infs(self):

        """
        Thisj function ...
        :return:
        """

        return self.replace_zeroes(inf_value)

    # -----------------------------------------------------------------

    def replace_zeroes_by_nans(self):

        """
        This function ...
        :return:
        """

        return self.replace_zeroes(nan_value)

    # -----------------------------------------------------------------

    def replace_negatives(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        negatives = self.negatives
        self._data[negatives] = value
        return negatives

    # -----------------------------------------------------------------

    def replace_positives(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        positives = self.positives
        self._data[positives] = value
        return positives

    # -----------------------------------------------------------------

    def replace_negatives_by_zeroes(self):

        """
        Thisn function ...
        :return:
        """

        return self.replace_negatives(zero_value)

    # -----------------------------------------------------------------

    def replace_positives_by_zeroes(self):

        """
        Thisn function ...
        :return:
        """

        return self.replace_positives(zero_value)

    # -----------------------------------------------------------------

    def replace_negatives_by_infs(self):

        """
        This function ...
        :return:
        """

        return self.replace_negatives(inf_value)

    # -----------------------------------------------------------------

    def replace_negatives_by_nans(self):

        """
        This function ...
        :return:
        """

        return self.replace_negatives(nan_value)

    # -----------------------------------------------------------------

    def replace(self, mask, value):

        """
        This function ...
        :param mask:
        :param value:
        :return:
        """

        # Replace
        self.data[mask] = value

    # -----------------------------------------------------------------

    def replace_where_equal_to(self, value, replacement):

        """
        This function ...
        :param value:
        :param replacement:
        :return:
        """

        # Get the mask
        mask = self.where(value)

        # Replace
        self.replace(mask, replacement)

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def replace_where_not_equal_to(self, value, replacement):

        """
        This function ...
        :param value:
        :param replacement:
        :return:
        """

        # Get the mask
        mask = self.where_not(value)

        # Replace
        self.replace(mask, replacement)

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def replace_where_greater_than(self, value, replacement):

        """
        This function ...
        :param value:
        :param replacement
        :return:
        """

        # Get the mask
        mask = self.where_greater_than(value)

        # Replace
        self.replace(mask, replacement)

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def cutoff_greater(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.replace_where_greater_than(value, value)

    # -----------------------------------------------------------------

    def cutoff_above_zero(self):

        """
        This function ...
        :return:
        """

        return self.cutoff_greater(zero_value)

    # -----------------------------------------------------------------

    def replace_where_smaller_than(self, value, replacement):

        """
        Thisnf unction ...
        :param value:
        :param replacement:
        :return:
        """

        # Get the mask
        mask = self.where_smaller_than(value)

        # Replace
        self.data[mask] = replacement

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def cutoff_smaller(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.replace_where_smaller_than(value, value)

    # -----------------------------------------------------------------

    def cutoff_below_zero(self):

        """
        Thisn function ...
        :return:
        """

        return self.cutoff_smaller(zero_value)

    # -----------------------------------------------------------------

    def replace_by_zeroes_where_greater_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.replace_where_greater_than(value, zero_value)

    # -----------------------------------------------------------------

    def replace_by_nans_where_equal_to(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.replace_where_equal_to(value, nan_value)

    # -----------------------------------------------------------------

    def replace_by_infs_where_equal_to(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.replace_where_equal_to(value, inf_value)

    # -----------------------------------------------------------------

    def replace_by_zeroes_where_equal_to(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.replace_where_equal_to(value, zero_value)

    # -----------------------------------------------------------------

    def replace_by_nans_where_greater_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.replace_where_greater_than(value, nan_value)

    # -----------------------------------------------------------------

    def replace_by_infs_where_greater_than(self, value):

        """
        Thisnf unction ...
        :param value:
        :return:
        """

        return self.replace_where_greater_than(value, inf_value)

    # -----------------------------------------------------------------

    def replace_by_zeroes_where_smaller_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.replace_where_smaller_than(value, zero_value)

    # -----------------------------------------------------------------

    def replace_by_nans_where_smaller_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.replace_where_smaller_than(value, nan_value)

    # -----------------------------------------------------------------

    def replace_by_infs_where_smaller_than(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.replace_where_smaller_than(value, inf_value)

    # -----------------------------------------------------------------

    def replace_by_zeroes(self, mask):

        """
        This function ...
        :param mask:
        :return:
        """

        self.replace(mask, zero_value)

    # -----------------------------------------------------------------

    def replace_by_nans(self, mask):

        """
        This function ...
        :param mask:
        :return:
        """

        self.replace(mask, nan_value)

    # -----------------------------------------------------------------

    def replace_by_infs(self, mask):

        """
        This function ...
        :param mask:
        :return:
        """

        self.replace(mask, inf_value)

    # -----------------------------------------------------------------

    def box_like(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        data = self._data[box.y_min:box.y_max, box.x_min:box.x_max]

        # Create the new box and return it
        return Cutout(data, box.x_min, box.x_max, box.y_min, box.y_max)

    # -----------------------------------------------------------------

    def to_rgba(self, interval="pts", scale="log", alpha="absolute", peak_alpha=1., colours="red", normalize_in=None,
                return_minmax=False, around_zero=False, symmetric=False):

        """
        This function ...
        :param interval:
        :param scale:
        :param alpha:
        :param peak_alpha:
        :param colours:
        :param normalize_in:
        :param return_minmax:
        :param around_zero:
        :param symmetric:
        :param interval:
        :return:
        """

        from .rgba import RGBAImage
        return RGBAImage.from_frame(self, interval=interval, scale=scale, alpha=alpha, peak_alpha=peak_alpha,
                                    colours=colours, normalize_in=normalize_in, return_minmax=return_minmax,
                                    around_zero=around_zero, symmetric=symmetric)

    # -----------------------------------------------------------------

    @property
    def file_size(self):

        """
        This function ...
        :return:
        """

        # Give warning message
        log.warning("The size of the frame in itself is perhaps not accurately represented by the size of the file it originates from, because it possibly contains muliple frames")

        # Check if path is defined
        if self.path is None: raise ValueError("Path is not defined")

        # Return the file size, if the frame has a path
        return fs.file_size(self.path)

    # -----------------------------------------------------------------

    @property
    def data_size(self):
        return (self._data.nbytes * u("byte")).to("GB")

    # -----------------------------------------------------------------

    @property
    def hdu(self):
        from astropy.io.fits import PrimaryHDU
        return PrimaryHDU(self.data, self.header)

    # -----------------------------------------------------------------

    def save(self, header=None, origin=None, extra_header_info=None, add_meta=True):

        """
        This function ...
        :param header:
        :param origin:
        :param extra_header_info:
        :param add_meta:
        :return:
        """

        # Inform the user
        log.info("Saving the frame ...")

        # Check if path is defined
        if self.path is None: raise RuntimeError("Path for the frame is not defined")

        # Check if original file was not multiplane
        if self._from_multiplane: raise RuntimeError("Cannot save frame into a multiplane image")

        # Save
        self.saveto(self.path, header=header, origin=origin, extra_header_info=extra_header_info, add_meta=add_meta)

    # -----------------------------------------------------------------

    def saveto(self, path, **kwargs):

        """
        This function ...
        :param path:
        :param kwargs:
        """

        if path.endswith("fits"): self.saveto_fits(path, **kwargs)
        elif path.endswith("asdf"): self.saveto_asdf(path, **kwargs)
        elif path.endswith("png"): self.saveto_png(path, **kwargs)
        else: raise ValueError("Unknown file format: " + fs.get_extension(path))

    # -----------------------------------------------------------------

    def saveto_fits(self, path, header=None, origin=None, extra_header_info=None, add_meta=True, update_path=True):

        """
        This function ...
        :param path:
        :param header:
        :param origin:
        :param extra_header_info:
        :param add_meta:
        :param update_path:
        :return:
        """

        if header is None: header = self.header

        # Set unit, FWHM and filter description
        if self.unit is not None:
            header.set("SIGUNIT", represent_unit(self.unit), "Unit of the map")
            header.set("PHYSTYPE", self.unit.physical_type, "Physical type of the unit")
        if self.fwhm is not None: header.set("FWHM", self.fwhm.to("arcsec").value, "[arcsec] FWHM of the PSF")

        # Set filter
        if self.filter is not None: header.set("FILTER", str(self.filter), "Filter used for this observation")
        else: header.set("FILTER", "n/a", "This image does not correspond to a certain observational filter")

        # Set pixelscale
        if not self.has_wcs and self.has_pixelscale:

            # Angular
            if self.has_angular_pixelscale:
                header.set("XPIXSIZE", repr(self.pixelscale.x.to("arcsec").value), "[arcsec] Pixelscale for x axis")
                header.set("YPIXSIZE", repr(self.pixelscale.y.to("arcsec").value), "[arcsec] Pixelscale for y axis")

            # Physical
            elif self.has_physical_pixelscale:
                header.set("XPIXSIZE", repr(self.pixelscale.x.to("pc").value), "[pc] Pixelscale for x axis")
                header.set("YPIXSIZE", repr(self.pixelscale.y.to("pc").value), "[pc] Pixelscale for y axis")

            # Invalid pixelscale
            else: raise RuntimeError("We shouldn't get here")

        # Set distance
        if self.distance is not None: header.set("DISTANCE", repr(self.distance.to("Mpc").value), "[Mpc] Distance to the object")

        # Set PSF FILTER
        if self.psf_filter is not None: header.set("PSFFLTR", str(self.psf_filter), "Filter to which the PSF of the frame corresponds")

        # Set smoothing factor
        if self.smoothing_factor != 1: header.set("SMOOTHF", repr(self.smoothing_factor), "Applied smoothing factor")

        # Add origin description
        if origin is not None: header["ORIGIN"] = origin
        else: header["ORIGIN"] = "Frame class of PTS package"

        # Add meta information
        if add_meta:
            for key in self.metadata:
                #print(key, self.metadata[key])
                header[key] = self.metadata[key]

        # Add extra info
        if extra_header_info is not None:
            for key in extra_header_info: header[key] = extra_header_info[key]

        # Write
        from .fits import write_frame

        #try:
        write_frame(self._data, header, path)
        #except IOError("Something went wrong during write")

        # Replace the path
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

        if header is None: header = self.header

        # Import
        from asdf import AsdfFile

        # Create the tree
        tree = dict()

        # Add data and header
        tree["data"] = self._data
        tree["header"] = header

        # Create the asdf file
        ff = AsdfFile(tree)

        # Write
        ff.write_to(path)

        # Update the path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    def saveto_png(self, path, interval="pts", scale="log", alpha="absolute", peak_alpha=1., colours="red", normalize_in=None):

        """
        This function ...
        :param path:
        :param interval:
        :param scale:
        :param alpha:
        :param peak_alpha:
        :param colours:
        :param normalize_in:
        :return:
        """

        # Get image values
        image, vmin, vmax = self.to_rgba(interval=interval, scale=scale, alpha=alpha, peak_alpha=peak_alpha,
                                         colours=colours, normalize_in=normalize_in, return_minmax=True)

        # Save
        image.saveto(path)

        # Return vmin and vmax
        return vmin, vmax

    # -----------------------------------------------------------------

    def to_hdu(self):

        """
        This function ...
        :return:
        """

        # Create and return the HDU
        return fits.ImageHDU(self.data, header=self.header)

# -----------------------------------------------------------------

def sum_frames(*args, **kwargs):

    """
    This function ...
    :param args:
    :return:
    """

    #arrays = [np.array(arg) for arg in args] # how it used to be
    arrays = [arg.data for arg in args]
    return Frame(np.sum(arrays, axis=0), **kwargs)

# -----------------------------------------------------------------

def sum_frames_quadratically(*args):

    """
    This function ...
    :param args:
    :return:
    """

    #arrays = [np.array(arg)**2 for arg in args]

    if len(args) == 0: raise ValueError("No frames are passed")

    arrays = []

    for arg in args:
        if arg is None: continue
        else: arrays.append(arg.data**2)

    # Return
    return Frame(np.sqrt(np.sum(arrays, axis=0)))

# -----------------------------------------------------------------

def linear_combination(frames, coefficients, checks=True):

    """
    This function ...
    :param frames:
    :param coefficients:
    :param checks:
    :return:
    """

    if checks:
        from .list import check_uniformity
        unit, wcs, pixelscale, psf_filter, fwhm, distance = check_uniformity(*frames)
    else: unit = wcs = pixelscale = psf_filter = fwhm = distance = None

    #print("unit", unit)
    #print("wcs", wcs)
    #print("pixelscale", pixelscale)
    #print("psf_filter", psf_filter)
    #print("fwhm", fwhm)
    #print("distance", distance)

    #print(coefficients)
    return sum_frames(*[frame*coefficient for coefficient, frame in zip(coefficients, frames)], unit=unit, wcs=wcs, pixelscale=pixelscale, psf_filter=psf_filter, fwhm=fwhm, distance=distance)

# -----------------------------------------------------------------

def log10(frame):

    """
    This function ...
    :param frame:
    :return:
    """

    return frame.get_log10()

# -----------------------------------------------------------------

def create_rotate_matrix(degs):

    """
    Return a rotation matrix for counterclockwise rotation by ``deg`` degrees.
    """

    rads = math.radians(degs)
    s = math.sin(rads)
    c = math.cos(rads)

    # Return the matrix
    return np.array([[c, -s], [s, c]])

# -----------------------------------------------------------------

def rotate_wcs(wcs, angle):

    """
    This function ...
    :param wcs:
    :param angle:
    :return:
    """

    if hasattr(wcs.wcs, 'cd'): wcs.wcs.cd = np.dot(create_rotate_matrix(angle.to("deg").value), wcs.wcs.cd)
    else: wcs.wcs.pc = np.dot(create_rotate_matrix(angle.to("deg").value), wcs.wcs.pc)

# -----------------------------------------------------------------

def regularize_frame(frame, absolute=False, dilate_nans=False, dilation_radius=5, dilation_niterations=1,
                     no_nans=False, no_infs=False, no_negatives=False, no_positives=False, cutoff_above=None,
                     cutoff_below=None, replace_nans=zero_value, replace_infs=zero_value):

    """
    This function ...
    :param frame:
    :param absolute:
    :param dilate_nans:
    :param dilation_radius:
    :param dilation_niterations:
    :param no_nans:
    :param no_infs:
    :param no_negatives:
    :param no_positives:
    :param cutoff_above:
    :param cutoff_below:
    :param replace_nans:
    :param replace_infs:
    :return:
    """

    # Absolute?
    if absolute: frame.make_absolute()

    # Dilate
    if dilate_nans:
        nans = frame.dilate_nans(radius=dilation_radius, niterations=dilation_niterations)
        # plotting.plot_mask(nans, title="previous nans")
        # plotting.plot_mask(frame.nans, title="new nans")

    # Remove nans or infs
    if no_infs: frame.replace_infs(replace_infs) #frame.replace_infs_by_zeroes()
    if no_nans: frame.replace_nans(replace_nans) #frame.replace_nans_by_zeroes()
    if no_negatives: frame.replace_negatives_by_zeroes()
    if no_positives: frame.replace_positives_by_zeroes()
    if cutoff_above is not None: frame.cutoff_greater(cutoff_above)
    if cutoff_below is not None: frame.cutoff_smaller(cutoff_below)

# -----------------------------------------------------------------

StaticFrame = create_lazified_class(Frame, "StaticFrame")

# -----------------------------------------------------------------
