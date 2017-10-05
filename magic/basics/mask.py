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
import copy
import numpy as np
from scipy import ndimage

# Import astronomical modules
from astropy.io import fits
from photutils import detect_sources

# Import the relevant PTS classes and modules
from .vector import Position
from ...core.basics.log import log
from .vector import Pixel

# -----------------------------------------------------------------

def get_mask_names(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..tools import headers

    # Load the header
    header = fits.getheader(path)

    # Get the number of planes
    nplanes = headers.get_number_of_frames(header)

    # Initialize a dictionary to contain the mask names and corresponding descriptions
    masks = dict()

    # Look at the properties of each plane
    for i in range(nplanes):

        # Get name and description of plane
        name, description, plane_type = headers.get_frame_name_and_description(header, i, always_call_first_primary=False)
        if plane_type == "mask": masks[name] = description

    # Return the masks with their name and description
    return masks

# -----------------------------------------------------------------

class MaskBase(object):

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
        self._data = data.astype(bool)

    # -----------------------------------------------------------------

    @classmethod
    def where(cls, data, value, **kwargs):

        """
        This function ...
        :param data:
        :param value:
        :param kwargs:
        :return:
        """

        return cls(data == value, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def above(cls, data, threshold, **kwargs):

        """
        This function ...
        :param data:
        :param threshold:
        :param kwargs:
        :return:
        """

        return cls(data > threshold, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def below(cls, data, threshold, **kwargs):

        """
        This function ...
        :param data:
        :param threshold:
        :param kwargs:
        :return:
        """

        return cls(data < threshold, **kwargs)

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
        data = region.to_mask(x_size, y_size)
        return cls(data)

    # -----------------------------------------------------------------

    @classmethod
    def intersection(cls, *args):

        """
        This function ...
        :param args:
        :return:
        """

        # Initialize data
        data = np.ones(args[0].shape)

        # Combine the masks
        for mask in args: data *= mask

        # Return the mask
        return cls(data)

    # -----------------------------------------------------------------

    @classmethod
    def union(cls, *args):

        """
        This function ...
        :param args:
        :return:
        """

        # Initialize data
        data = np.zeros(args[0].shape)

        # Combine the masks
        for mask in args: data += mask

        # Return the mask
        return cls(data)

    # -----------------------------------------------------------------

    @property
    def data(self):

        """
        This function ...
        :return:
        """

        return self._data

    # -----------------------------------------------------------------

    def astype(self, dtype):

        """
        This function ...
        :param dtype: 
        :return: 
        """

        return self.data.astype(dtype)

    # -----------------------------------------------------------------

    #@property
    #def array(self):

        #"""
        #The underlying array
        #"""

        #return self._data

    # -----------------------------------------------------------------

    def __array__(self):

        """
        Array representation of the mask (e.g., for matplotlib).
        """

        return np.array(self.data)

    # -----------------------------------------------------------------

    def __array_prepare__(self, array, context=None):

        """
        I don't know if this function should be here
        """

        return array

    # -----------------------------------------------------------------

    @property
    def dtype(self):

        """
        `numpy.dtype` of this object's data.
        """

        return self.data.dtype

    # -----------------------------------------------------------------

    def __getitem__(self, item):

        """
        This function ...
        :param item:
        :return:
        """

        if isinstance(item, MaskBase): return self._data[item.data]
        elif isinstance(item, Pixel): return self._data[item.y, item.x]
        elif isinstance(item, tuple):
            #return self._data[item] #  ## WHY IS THIS NOT WORKING?????!!!
            #print(item)
            #print(self.shape)
            if isinstance(item[0], slice) and isinstance(item[1], slice):
                slice_y = item[0]
                slice_x = item[1]
                #indices_y = slice_y.indices(self.ysize)
                #indices_x = slice_x.indices(self.xsize)
                #print(indices_x)
                #print(indices_y)
                return self._data[slice_y, slice_x]
            elif isinstance(item[0], int) and isinstance(item[1], int):
                #print(self.shape)
                #print(item)
                #try:
                return self._data.__getitem__(item)
                #except IndexError:
            else: raise NotImplementedError("Not implemented")
            #return self._data[item[0], item[1]]
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
        elif isinstance(item, tuple):
            #print(self._data)
            #print(type(self._data))
            #print(item)
            #self._data[0, 0] = value
            #self._data[item[0], item[1]] = value
            ## WHY IS THIS NOT WORKING?????!!!
            # Error: TypeError: __array__() takes exactly 1 argument (2 given)
            try: self._data.__setitem__(item, value)
            except TypeError:
                #print(item)
                #print(value)
                #data = np.array(self._data)
                slice_y = item[0]
                slice_x = item[1]
                indices_y = slice_y.indices(self.ysize)
                indices_x = slice_x.indices(self.xsize)
                for index_y, y in enumerate(indices_y):
                    for index_x, x in enumerate(indices_x):
                        self._data[y, x] = value[index_y, index_x]
                #for index in range(indices[0], indices[1], indices[2]):
        else: self._data[item] = value

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return self._data.shape

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
    def x_center(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.xsize - 1)

    # -----------------------------------------------------------------

    @property
    def y_center(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.ysize - 1)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        Thisf unction ....
        :return:
        """

        from .coordinate import PixelCoordinate
        return PixelCoordinate(self.x_center, self.y_center)

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

    def dilate_rc(self, rank, connectivity, iterations=1):

        """
        This function ...
        :param rank: 2?
        :param connectivity:
        :param iterations:
        :return:
        """

        # Define the structure for the expansion
        structure = ndimage.generate_binary_structure(rank, connectivity=connectivity)

        # Dilate
        self.dilate(structure, iterations)

    # -----------------------------------------------------------------

    def dilate(self, structure, iterations=1):

        """
        This function ...
        :return:
        """

        # Dilate
        self._data = ndimage.binary_dilation(self._data, structure, iterations)

    # -----------------------------------------------------------------

    def erode_rc(self, rank, connectivity, iterations=1):

        """
        This function ...
        :param rank: 2?
        :param connectivity:
        :param iterations:
        :return:
        """

        # Create structure
        structure = ndimage.generate_binary_structure(rank, connectivity=connectivity)

        # Erode
        self.erode(structure, iterations)

    # -----------------------------------------------------------------

    def erode(self, structure, iterations):

        """
        This function ...
        :return:
        """

        # Erode
        self._data = ndimage.binary_erosion(self._data, structure, iterations)

    # -----------------------------------------------------------------

    def eroded_rc(self, rank=2, connectivity=2, iterations=1):

        """
        This function ...
        :param rank: 
        :param connectivity: 
        :param iterations: 
        :return: 
        """

        new = self.copy()
        new.erode_rc(rank, connectivity, iterations)
        return new

    # -----------------------------------------------------------------

    def eroded(self, structure, iterations):

        """
        This function ...
        :param structure: 
        :param iterations: 
        :return: 
        """

        new = self.copy()
        new.erode(structure, iterations)
        return new

    # -----------------------------------------------------------------

    def open(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def close(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def __add__(self, other):

        """
        This function ...
        :return:
        """

        new = self.copy()
        if isinstance(other, MaskBase): new._data += other.data
        else: new._data += other
        return new

    # -----------------------------------------------------------------

    def __iadd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        if isinstance(other, MaskBase): self._data += other.data
        else: self._data += other
        return self

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        new = self.copy()
        if isinstance(other, MaskBase): new._data *= other.data
        else: new._data *= other
        return new

    # -----------------------------------------------------------------

    def __imul__(self, other):

        """
        This function ...
        :return:
        """

        if isinstance(other, MaskBase): self._data *= other.data
        else: self._data *= other
        return self

    # -----------------------------------------------------------------

    def invert(self):

        """
        This function ...
        :return:
        """

        self._data = np.logical_not(self._data)

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

    def central(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels: number of connected pixels
        :param connectivity:
        :return:
        """

        # Perform the segmentation
        segments = detect_sources(self.data.astype(float), 0.5, npixels=npixels, connectivity=connectivity).data

        # To plot the multiple segments that are detected
        # if segments.max() > 1: plotting.plot_box(np.ma.masked_array(box, mask=segments.astype(bool)))

        # Get the label of the center segment
        #rel_center = self.cutout.rel_position(self.center)
        label = segments[int(round(self.y_center)), int(round(self.x_center))]

        # If the center pixel is identified as being part of the background, create an empty mask (the center does not
        # correspond to a segment)
        if label == 0: return self.__class__(np.zeros_like(self.data, dtype=bool))

        # Create a mask of the center segment
        else: return self.__class__(segments == label)

    # -----------------------------------------------------------------

    def largest(self, npixels=1, connectivity=8):

        """
        This function ...
        :param npixels:
        :param connectivity:
        :return:
        """

        # Perform the segmentation
        segments = detect_sources(self.data.astype(float), 0.5, npixels=npixels, connectivity=connectivity).data

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
            if len(sorted_indices) == 1: return self.__class__(np.zeros_like(self.data, dtype=bool))

            second_last_index = sorted_indices[-2]
            label = unique[second_last_index]

        else: label = value

        # Return the patch as a mask
        return self.__class__(segments == label)

    # -----------------------------------------------------------------

    def fill_holes(self):

        """
        This function ...
        :return:
        """

        # Create a copy of this mask
        #new_mask = self.copy()

        # Perform the segmentation
        segments = detect_sources(self.inverse().data.astype(float), 0.5, 1).data

        # Find the label of the largest segment (=the background)
        label_counts = np.bincount(segments.flatten())
        if len(label_counts) > 1:

            background_label = np.argmax(label_counts[1:]) + 1
            # If the source mask is larger than the background (in number of pixels), the above will provide the correct label
            # therefore we do the '[1:]'

            # Create a mask for the holes identified as background
            holes = self.inverse().data
            holes[segments == background_label] = False

            # Remove holes from the mask
            self._data[holes] = True

        # Return the new mask
        #return new_mask

    # -----------------------------------------------------------------

    @property
    def nmasked(self):

        """
        This function ...
        :return:
        """

        return np.sum(self.data)

    # -----------------------------------------------------------------

    @property
    def nunmasked(self):

        """
        This function ...
        :return:
        """

        return np.sum(np.logical_not(self.data))

    # -----------------------------------------------------------------

    def masks(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Calculate x and y of the pixel corresponding to the object's position
        pixel = Pixel.for_coordinate(position, round_first=True)

        # Check whether pixel exists in this mask
        if not pixel.exists_in(self): return False

        # Return whether the corresponding pixel is masked
        return self.data[pixel.y, pixel.x]  # Return the value of the mask in this pixel

    # -----------------------------------------------------------------

    def covers(self, other_mask):

        """
        This function ...
        :param other_mask:
        :return:
        """

        if isinstance(other_mask, MaskBase): other = other_mask.data
        else: other = other_mask

        not_covered = other & self.inverse().data
        return not np.any(not_covered)

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
        mask = cls(shape.to_mask(x_size, y_size))

        # Return the mask (inverted if requested)
        if invert: return mask.inverse()
        else: return mask

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
    def from_file(cls, path, index=0, plane=None):

        """
        This function ...
        :param path:
        :param index:
        :param plane:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file '" + path + "' ...")

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

            if plane is not None:

                from ..tools import headers

                for i in range(nframes):
                    # Get name and description of frame
                    name, description, plane_type = headers.get_frame_name_and_description(header, i, always_call_first_primary=False)
                    if plane == name and plane_type == "mask":
                        index = i
                        break

                # If a break is not encountered, a matching plane name is not found
                else: raise ValueError("Plane with name '" + plane + "' not found")

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
    def is_inf(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        return cls(np.isinf(frame))

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
    def xsize(self):

        """
        This property ...
        :return:
        """

        return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This property ...
        :return:
        """

        return self.shape[0]

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

        from skimage import morphology  ###### TODO: this can cause a very weird error: (on Nancy (Ubuntu 14.04.4 LTS) with NUMPY VERSION 1.9.0)
        # *** libmkl_mc3.so *** failed with error : /raid6/home/sjversto/Enthought/Canopy_64bit/User/bin/../lib/libmkl_mc3.so: undefined symbol: i_free
        # *** libmkl_def.so *** failed with error : /raid6/home/sjversto/Enthought/Canopy_64bit/User/bin/../lib/libmkl_def.so: undefined symbol: i_free
        # MKL FATAL ERROR: Cannot load neither libmkl_mc3.so nor libmkl_def.so
        # POTENTIAL FIX HERE: http://stackoverflow.com/questions/14495334/python-matplotlib-mkl-fatal-error-on-ubuntu-12-04
        # IT WORKS WITH NUMPY VERSION 1.8.1 !!!

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

        from skimage import morphology  ###### TODO: this can cause a very weird error: (on Nancy (Ubuntu 14.04.4 LTS) with NUMPY VERSION 1.9.0)
        # *** libmkl_mc3.so *** failed with error : /raid6/home/sjversto/Enthought/Canopy_64bit/User/bin/../lib/libmkl_mc3.so: undefined symbol: i_free
        # *** libmkl_def.so *** failed with error : /raid6/home/sjversto/Enthought/Canopy_64bit/User/bin/../lib/libmkl_def.so: undefined symbol: i_free
        # MKL FATAL ERROR: Cannot load neither libmkl_mc3.so nor libmkl_def.so
        # POTENTIAL FIX HERE: http://stackoverflow.com/questions/14495334/python-matplotlib-mkl-fatal-error-on-ubuntu-12-04
        # IT WORKS WITH NUMPY VERSION 1.8.1 !!!

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
        else: return self[y_pixel, x_pixel] # Return the value of the mask in this pixel

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
