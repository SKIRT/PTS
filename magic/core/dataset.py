#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.dataset Contains the DataSet class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import numpy as np
from collections import OrderedDict

# Import astronomical modules
from astropy.io import fits
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .frame import Frame
from .fits import get_frame_names, get_mask_names, get_plane_names
from ...core.basics.log import log
from .datacube import DataCube
from .image import Image
from .mask import Mask
from ..basics.coordinatesystem import CoordinateSystem
from ...core.basics.configurable import Configurable
from ..region.list import SkyRegionList
from ..tools import headers
from .list import NamedImageList, NamedFrameList, FrameList, ImageList
from ...core.tools import types
from ...core.filter.filter import parse_filter
from .mask import intersection, union
from ..region.rectangle import SkyRectangleRegion
from ..basics.coordinate import SkyCoordinate
from ..basics.stretch import SkyStretch
from ...core.basics.map import Map
from ..tools import coordinates
from ...core.units.parsing import parse_unit as u
from ...core.basics.table import SmartTable
from ...core.basics.containers import FileList, NamedFileList
from ...core.basics.range import QuantityRange
from ...core.tools.utils import create_lazified_class
from ...core.tools import sequences
from ..basics.pixelscale import Pixelscale
from ...core.tools import strings

# -----------------------------------------------------------------

class DataSetTable(SmartTable):

    """
    This function ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Name"] = (str, None, "Name for the fitting run")
    _column_info["Path"] = (str, None, "Name of the model used")
    _column_info["Error path"] = (str, None, "Path of the error map")
    _column_info["Mask path"] = (str, None, "Path of the mask")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DataSetTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_entry(self, name, path, error_path, mask_path):

        """
        THis function ...
        :param name: 
        :param path: 
        :param error_path: 
        :param mask_path: 
        :return: 
        """

        values = [name, path, error_path, mask_path]
        self.add_row(values)

# -----------------------------------------------------------------

class DataSet(object):

    """
    This class...
    """

    default_extension = "dat"

    # -----------------------------------------------------------------

    def __init__(self):

        """
        The constructor ...
        """

        # The path to the dataset file
        self.path = None

        # The paths to the images
        self.paths = OrderedDict()

        # The paths to the error maps
        self.error_paths = dict()

        # The paths to the masks
        self.mask_paths = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_cwd(cls, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        return cls.from_directory(fs.cwd(), **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, **kwargs):

        """
        This function ...
        :param path:
        :param kwargs:
        :return:
        """

        # Get function that obtains name for the dataset from the filename
        if "get_name" in kwargs: get_name = kwargs.pop("get_name")
        else: get_name = None

        # Add kwargs
        kwargs["extension"] = "fits"
        kwargs["returns"] = ["path", "name"]

        # Set the paths
        #paths, names = fs.files_in_path(path, **kwargs)

        # Create a new dataset instance
        dataset = cls()

        # Add the paths
        #for path, name in zip(paths, names):
        for path, name in fs.files_in_path(path, **kwargs):

            # Get actual name from filename
            if get_name is not None: name = get_name(name)

            # Add the image path
            dataset.add_path(name, path)

        # Return the dataset instance
        return dataset

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, check=True):

        """
        This function ...
        :param path:
        :param check:
        :return:
        """

        # Construct a new dataset
        dataset = cls()

        # Load the table
        #table = tables.from_file(path)
        table = DataSetTable.from_file(path)

        # Loop over the entries in the table
        for i in range(len(table)):

            # Get the paths
            name = table["Name"][i]
            path = table["Path"][i]
            error_path = table["Error path"][i] if not table["Error path"].mask[i] else None
            mask_path = table["Mask path"][i] if not table["Mask path"].mask[i] else None

            # Add the paths to the dataset
            dataset.add_path(name, path, check=check)
            if error_path is not None: dataset.add_error_path(name, error_path, check=check)
            if mask_path is not None: dataset.add_mask_path(name, mask_path, check=check)

        # Return the dataset
        return dataset

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        for name in self.paths: yield name

    # -----------------------------------------------------------------

    def __len__(self):
        return len(self.paths)

    # -----------------------------------------------------------------

    def get_depending_paths(self):

        """
        This function ...
        :return:
        """

        paths = dict()

        # Loop over the image paths
        for name in self.paths:

            # Generate a unique name
            unique_name = "image_" + name

            # Set the path
            paths[unique_name] = self.paths[name]

        # Loop over the error map paths
        for name in self.error_paths:

            # Generate a unique name
            unique_name = "error_" + name

            # Set the path
            paths[unique_name] = self.error_paths[name]

        # Loop over the mask paths
        for name in self.mask_paths:

            # Generate a unique name
            unique_name = "mask_" + name

            # Set the path
            paths[unique_name] = self.mask_paths[name]

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    def set_depending_path(self, label, path):

        """
        This function ...
        :param label:
        :param path:
        :return:
        """

        # Image path
        if label.startswith("image_"):

            # Get name
            name = strings.split_at_first(label, "image_")[1]

            # Check if defined
            if name not in self.paths: raise ValueError("No such image: '" + name + "'")

            # Replace the iamge path
            self.paths[name] = path

        # Error map path
        elif label.startswith("error_"):

            # Get error map name
            name = strings.split_at_first(label, "error_")[1]

            # Check if defined
            if name not in self.error_paths: raise ValueError("No such error map: '" + name + "'")

            # Replace the error map path
            self.error_paths[name] = path

        # Mask path
        elif label.startswith("mask_"):

            # Get the mask name
            name = strings.split_at_first(label, "mask_")[1]

            # Check if defined
            if name not in self.mask_paths: raise ValueError("No such mask: '" + name + "'")

            # Replace the mask path
            self.mask_paths[name] = path

        # Invalid label
        else: raise ValueError("Invalid label")

    # -----------------------------------------------------------------

    def set_depending_paths(self, paths):

        """
        This function ...
        :param paths:
        :return:
        """

        # Set all depending paths
        for label in paths: self.set_depending_path(label, paths[label])

    # -----------------------------------------------------------------

    def copy(self):
        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    @property
    def names(self):
        return self.paths.keys()

    # -----------------------------------------------------------------

    @property
    def path_list(self):
        return self.paths.values()

    # -----------------------------------------------------------------

    def add_path(self, name, path, check=True):

        """
        This function ...
        :param name:
        :param path:
        :param check:
        :return:
        """

        # Check if the file exists
        if check and not fs.is_file(path): raise IOError("File doesn't exist: '" + path + "'")

        # Check if already such a name
        if name in self.paths: raise ValueError("Already a path in the dataset with the name " + name)

        # Add the path
        self.paths[name] = path

    # -----------------------------------------------------------------

    def add_error_path(self, name, path, check=True):

        """
        This function ...
        :param name:
        :param path:
        :param check:
        :return:
        """

        # Check if the file exists
        if check and not fs.is_file(path): raise IOError("File doesn't exist: '" + path + "'")

        # Check if corresponding frame exists
        if name not in self.paths: raise ValueError("Corresponding image with name " + name + " has not been added")

        # Check if already such a name
        if name in self.error_paths: raise ValueError("Already an error path in the dataset with the name " + name)

        # Add the path
        self.error_paths[name] = path

    # -----------------------------------------------------------------

    def add_mask_path(self, name, path, check=True):

        """
        This function ...
        :param name:
        :param path:
        :param check:
        :return:
        """

        # Check if the file exists
        if check and not fs.is_file(path): raise IOError("File does not exist: '" + path + "'")

        # Check if the corresponding frame exists
        if name not in self.paths: raise ValueError("Corresponding image with name " + name + " has not been added")

        # Check if already such a name
        if name in self.mask_paths: raise ValueError("Already a mask path in the dataset with the name " + name)

        # Add the path
        self.mask_paths[name] = path

    # -----------------------------------------------------------------

    @property
    def filters(self):
        return self.get_filters().values()

    # -----------------------------------------------------------------

    @property
    def filters_from_names(self):
        return self.get_filters_from_names().values()

    # -----------------------------------------------------------------

    def get_filters(self):

        """
        This function ...
        :return:
        """

        # Initialize
        fltrs = OrderedDict()

        # Loop over the images
        #for name in self.paths: fltrs[name] = self.get_frame(name, masked=False).filter
        for name in self.paths: fltrs[name] = self.get_filter(name)

        # Return the dictionary with the filters
        return fltrs

    # -----------------------------------------------------------------

    def get_filters_from_names(self):

        """
        Thins function ...
        :return:
        """

        # Initialize
        fltrs = OrderedDict()

        # Loop over the iamges
        for name in self.paths: fltrs[name] = parse_filter(name)

        # Return
        return fltrs

    # -----------------------------------------------------------------

    def get_filter(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        header = self.get_header(name)
        return headers.get_filter(name, header)
        #return self.get_frame(name).filter

    # -----------------------------------------------------------------

    def get_filter_name(self, name):

        """
        Thisnf ucntion ...
        :param name:
        :return:
        """

        return str(self.get_filter(name))

    # -----------------------------------------------------------------

    def get_pixelscale(self, name, average=False, unit=None):

        """
        This function ...
        :param name:
        :param average:
        :param unit:
        :return:
        """

        header = self.get_header(name)
        pixelscale = headers.get_pixelscale(header)

        # Create coordinate system first
        if pixelscale is None: pixelscale = CoordinateSystem.from_header(header).pixelscale

        if average:
            if unit is not None: return pixelscale.average.to(unit)
            else: return pixelscale.average
        else:
            if unit is not None: raise ValueError("Cannot specify unit when average is not True")
            return pixelscale

    # -----------------------------------------------------------------

    def get_wavelength(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_frame(name).wavelength

    # -----------------------------------------------------------------

    def get_header(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Load the header and return it
        return fits.getheader(self.paths[name])

    # -----------------------------------------------------------------

    def get_frame_path(self, name):

        """
        This function ...
        :param name: 
        :return: 
        """

        return self.paths[name]

    # -----------------------------------------------------------------

    def get_frame(self, name, masked=True, mask_value=0.0):

        """
        This function ...
        :param name:
        :param masked:
        :param mask_value:
        :return:
        """

        # Open the frame and return it
        frame = Frame.from_file(self.paths[name])

        # Check if the frame has to be masked
        if masked:

            # Get the mask and set frame pixels to zero
            if name in self.mask_paths:
                mask = self.get_mask(name)
                frame[mask] = mask_value
            #else: log.warning("No mask available for " + name + " frame")

        # Set the name
        frame.name = name

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def get_directory_paths(self, min_wavelength=None, max_wavelength=None, exclude=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :return:
        """

        # Make sure exclude is a list
        if types.is_string_type(exclude): exclude = [exclude]
        elif exclude is None: exclude = []

        # Initialize a dictionary
        paths = OrderedDict()

        # Loop over the frame paths
        for name in self.paths:

            # Skip if name is in the exclude list
            if exclude is not None and name in exclude: continue

            # Check wavelength of the frame if necessary
            if min_wavelength is not None or max_wavelength is not None:

                # Load the frame
                frame = self.get_frame(name)
                #header = self.get_header(name)

                # Skip images of wavelength smaller than the minimum or greater than the maximum
                if min_wavelength is not None and frame.wavelength < min_wavelength: continue
                if max_wavelength is not None and frame.wavelength > max_wavelength: continue

            # Determine the path of the directory where the frame is in
            directory_path = fs.directory_of(self.paths[name])

            # Set the directory path
            paths[name] = directory_path

        # Return the dictionary of directory paths
        return paths

    # -----------------------------------------------------------------

    def get_frames(self, masked=True, mask_value=0.0, min_wavelength=None, max_wavelength=None, exclude=None, returns="dict", keys="name"):

        """
        This function ...
        :param masked:
        :param mask_value:
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :param returns:
        :param keys:
        :return:
        """

        # Make sure exclude is a list
        if types.is_string_type(exclude): exclude = [exclude]
        elif exclude is None: exclude = []

        # Initialize a dictionary for the frames
        frames = OrderedDict()

        # Loop over the frame paths
        for name in self.paths:

            # Skip if name is in the exclude list
            if exclude is not None and name in exclude: continue

            # Load the frame
            frame = self.get_frame(name, masked, mask_value)

            # Skip images of wavelength smaller than the minimum or greater than the maximum
            if min_wavelength is not None and frame.wavelength < min_wavelength: continue
            if max_wavelength is not None and frame.wavelength > max_wavelength: continue

            # Set the key
            if keys == "name": key = name
            elif keys == "filter": key = frame.filter
            elif keys == "filter_name": key = frame.filter_name
            elif keys == "wavelength": key = frame.wavelength
            else: raise ValueError("Invalid value for 'keys'")

            # Check if the key is not None
            if key is None: raise ValueError("Cannot use '" + keys + "' as key, because it cannot be obtained for frame '" + name + "'")

            # Add the frame
            frames[key] = frame

        # Return the dictionary of frames
        if returns == "dict": return frames
        elif returns == "list": return frames.values()
        else: raise ValueError("Invalid value for 'returns'")

    # -----------------------------------------------------------------

    def get_frame_paths(self, min_wavelength=None, max_wavelength=None, exclude=None, returns="dict", keys="name"):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :param returns:
        :param keys: 'name' or 'filter'
        :return:
        """

        # Make sure exclude is a list
        if types.is_string_type(exclude): exclude = [exclude]
        elif exclude is None: exclude = []

        # Initialize a dictionary for the frames
        paths = OrderedDict()

        # Loop over the frame paths
        for name in self.paths:

            # Skip if name is in the exclude list
            if exclude is not None and name in exclude: continue

            # Skip images of wavelength smaller than the minimum or greater then the maximum
            if min_wavelength is not None or max_wavelength is not None: wavelength = self.get_wavelength(name)
            else: wavelength = None
            if min_wavelength is not None and wavelength < min_wavelength: continue
            if max_wavelength is not None and wavelength > max_wavelength: continue

            # Set the key
            if keys == "name": key = name
            elif keys == "filter": key = self.get_filter(name)
            elif keys == "filter_name": key = self.get_filter_name(name)
            elif keys == "wavelength": key = wavelength if wavelength is not None else self.get_wavelength(name)
            else: raise ValueError("Invalid value for 'keys'")

            # Check if the key is not None
            if key is None: raise ValueError("Cannot use '" + keys + "' as key, because it cannot be obtained for frame '" + name + "'")

            # Add the path
            paths[key] = self.paths[name]

        # Return the dictionary of paths
        if returns == "dict": return paths
        elif returns == "list": return paths.values()
        else: raise ValueError("Invalid value for 'returns'")

    # -----------------------------------------------------------------

    def get_errormaps(self, masked=True, mask_value=0.0, min_wavelength=None, max_wavelength=None, exclude=None, returns="dict", keys="name"):

        """
        This function ...
        :param masked:
        :param mask_value:
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :param returns:
        :param keys:
        :return:
        """

        # Make sure exclude is a list
        if types.is_string_type(exclude): exclude = [exclude]
        elif exclude is None: exclude = []

        # Initialize a dictionary for the error maps
        errormaps = OrderedDict()

        # Loop over the frame paths
        for name in self.paths:

            # Skip if name is in the exclude list
            if exclude is not None and name in exclude: continue

            # Load the error map
            errormap = self.get_errormap(name, masked, mask_value)

            # Skip images of wavelength smaller than the minimum or greater than the maximum
            if min_wavelength is not None and errormap.wavelength < min_wavelength: continue
            if max_wavelength is not None and errormap.wavelength > max_wavelength: continue

            # Set the key
            if keys == "name": key = name
            elif keys == "filter": key = errormap.filter
            elif keys == "wavelength": key = errormap.wavelength
            else: raise ValueError("Invalid value for 'keys'")

            # Check if the key is not None
            if key is None: raise ValueError("Cannot use '" + keys + "' as key, because it cannot be obtained for frame '" + name + "'")

            # Add the error map
            errormaps[key] = errormap

        # Return the dictionary of error maps
        if returns == "dict": return errormaps
        elif returns == "list": return errormaps.values()
        else: raise ValueError("Invalid value for 'returns'")

    # -----------------------------------------------------------------

    def get_images(self, exclude=None, min_wavelength=None, max_wavelength=None, returns="dict", keys="name"):

        """
        This function ...
        :param exclude:
        :param min_wavelength:
        :param max_wavelength:
        :param returns:
        :param keys:
        :return: 
        """

        # Make sure exclude is a list
        if types.is_string_type(exclude): exclude = [exclude]
        elif exclude is None: exclude = []

        # Initialize a dictionary for the images
        images = OrderedDict()

        # Loop over the frame paths
        for name in self.paths:

            # Skip if name is in the exclude list
            if exclude is not None and name in exclude: continue

            # Load the image
            image = self.get_image(name)

            # Skip images of wavelengths smaller than the minimum or greater than the maximum
            if min_wavelength is not None and image.wavelength < min_wavelength: continue
            if max_wavelength is not None and image.wavelength > max_wavelength: continue

            # Set the key
            if keys == "name": key = name
            elif keys == "filter": key = image.filter
            elif keys == "wavelength": key = image.wavelength
            else: raise ValueError("Invalid value for 'keys'")

            # Check if the key is not None
            if key is None: raise ValueError("Cannot use '" + keys + "' as key, because it cannot be obtained for frame '" + name + "'")

            # Add the image
            images[key] = image

        # Return the dictionary of images
        if returns == "dict": return images
        elif returns == "list": return images.values()
        else: raise ValueError("Invalid value for 'returns'")

    # -----------------------------------------------------------------

    def get_framelist(self, masked=True, mask_value=0.0, min_wavelength=None, max_wavelength=None, exclude=None, named=True, lazy=True):

        """
        This function ...
        :param masked:
        :param mask_value:
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :param named:
        :param lazy:
        :return: 
        """

        # Lazy frame list
        if lazy:

            log.warning("Applying the masks is not supported (yet) for lazy frame lists")
            if named:
                arguments = self.get_frame_paths(min_wavelength, max_wavelength, exclude, keys="name", returns="dict")
                arguments["lazy"] = True
                return NamedFrameList.from_paths(**arguments)
            else:
                arguments = self.get_frame_paths(min_wavelength, max_wavelength, exclude, keys="filter_name", returns="dict")
                arguments["lazy"] = True
                return FrameList.from_paths(**arguments)

        # Not lazy
        else:

            if named: return NamedFrameList.from_dictionary(self.get_frames(masked, mask_value, min_wavelength, max_wavelength, exclude))
            else: return FrameList(*self.get_frames(masked, mask_value, min_wavelength, max_wavelength, exclude, returns="list"))

    # -----------------------------------------------------------------

    get_frame_list = get_framelist

    # -----------------------------------------------------------------

    def get_frame_path_list(self, min_wavelength=None, max_wavelength=None, exclude=None, named=True):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :param named:
        :return:
        """

        # Keyed on name or filter
        if named: return NamedFileList(**self.get_frame_paths(min_wavelength, max_wavelength, exclude))
        else: return FileList(**self.get_frame_paths(min_wavelength, max_wavelength, exclude, keys="filter"))

    # -----------------------------------------------------------------

    def get_errormaplist(self, masked=True, mask_value=0.0, min_wavelength=None, max_wavelength=None, exclude=None, named=True):

        """
        This function ...
        :param masked:
        :param mask_value:
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :param named:
        :return: 
        """

        if named: return NamedFrameList.from_dictionary(self.get_errormaps(masked, mask_value, min_wavelength, max_wavelength, exclude))
        else: return FrameList(*self.get_errormaps(masked, mask_value, min_wavelength, max_wavelength, exclude, returns="list"))

    # -----------------------------------------------------------------

    get_errormap_list = get_errormaplist

    # -----------------------------------------------------------------

    def get_imagelist(self, exclude=None, min_wavelength=None, max_wavelength=None, named=True):

        """
        This function ...
        :param exclude: 
        :param min_wavelength: 
        :param max_wavelength:
        :param named:
        :return: 
        """

        if named: return NamedImageList.from_dictionary(self.get_images(exclude, min_wavelength, max_wavelength))
        else: return ImageList(*self.get_images(exclude, min_wavelength, max_wavelength, returns="list"))

    # -----------------------------------------------------------------

    get_image_list = get_imagelist

    # -----------------------------------------------------------------

    def get_name_for_filter(self, fltr, return_none=True):

        """
        This function ...
        :param fltr:
        :param return_none:
        :return:
        """

        # Give warning
        log.warning("The 'get_name_for_filter' function is very slow, so its use is discouraged")

        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        filter_string = str(fltr)

        # Get the dictionary of filters
        fltrs = self.get_filters()

        # Loop over the filters dictionary
        for name in fltrs:

            # Get the filter
            frame_filter = self.get_filter(name)

            # Check the filter
            if str(frame_filter) == filter_string: return name

        # No frame found for this filter
        if return_none: return None
        else: raise ValueError("No image for the '" + filter_string + "' in the dataset")

    # -----------------------------------------------------------------

    def get_frame_path_for_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        name = self.get_name_for_filter(fltr)
        return self.get_frame_path(name)

    # -----------------------------------------------------------------

    def get_errormap_path_for_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        name = self.get_name_for_filter(fltr)
        return self.get_errormap_path(name)

    # -----------------------------------------------------------------

    def get_names_for_filters(self, filters, as_dict=False):

        """
        This function ...
        :param filters:
        :param as_dict:
        :return: 
        """

        # FASTER IMPLEMENTATION (NOT USING THE SLOW GET_NAME_FOR_FILTER)
        # ONLY OPENENING EACH IMAGE HEADER ONCE!

        # Make sure that filters are actual filter objects and not filter name
        filters = [parse_filter(fltr) if types.is_string_type(fltr) else fltr for fltr in filters]

        # Initialize
        if as_dict: names = dict()
        else: names = [None] * len(filters)

        # Loop over the names
        for name in self.paths:

            # Get the filter
            frame_filter = self.get_filter(name)

            if as_dict:

                if frame_filter not in filters: continue
                names[frame_filter] = name

            else:

                # Determine the index for this filter
                index = sequences.find_index(filters, frame_filter)
                if index is None: continue

                # Set the element
                names[index] = name

        # Return the names
        return names

    # -----------------------------------------------------------------

    def get_frames_for_filters(self, filters, as_dict=False):

        """
        Thisf unction ...
        :param filters:
        :param as_dict:
        :return:
        """

        # Initialize
        frames = OrderedDict()

        # Get frame names
        names = self.get_names_for_filters(filters)

        # Add the frames
        for name, fltr in zip(names, filters): frames[fltr] = self.get_frame(name)

        # Return
        if as_dict: return frames
        else: return frames.values()

    # -----------------------------------------------------------------

    def get_frame_paths_for_filters(self, filters):

        """
        This function ...
        :param filters: 
        :return: 
        """

        paths = []
        names = self.get_names_for_filters(filters)
        #for fltr in filters: paths.append(self.get_frame_path_for_filter(fltr))
        for name in names: paths.append(self.get_frame_path(name))
        return paths

    # -----------------------------------------------------------------

    def get_frame_names_and_paths_for_filters(self, filters):

        """
        This function ...
        :param filters: 
        :return: 
        """

        names = self.get_names_for_filters(filters, as_dict=True)

        result = OrderedDict()
        for fltr in filters:

            name = names[fltr]
            result[name] = self.get_frame_path(name) #self.get_frame_path_for_filter(fltr)

        return result

    # -----------------------------------------------------------------

    def get_framelist_for_filters(self, filters, named=True):

        """
        This function ...
        :param filters:
        :param named:
        :return: 
        """

        if named: return NamedFrameList.from_paths(**self.get_frame_names_and_paths_for_filters(filters))
        else: return FrameList.from_paths(*self.get_frame_paths_for_filters(filters))

    # -----------------------------------------------------------------

    def get_errormap_paths_for_filters(self, filters):

        """
        This function ...
        :param filters: 
        :return: 
        """

        paths = []
        for fltr in filters: paths.append(self.get_errormap_path_for_filter(fltr))
        return paths

    # -----------------------------------------------------------------

    def get_errormap_names_and_paths_for_filters(self, filters):

        """
        This funtion ...
        :param filters: 
        :return: 
        """

        result = OrderedDict()
        for fltr in filters:
            result[self.get_name_for_filter(fltr)] = self.get_errormap_path_for_filter(fltr)
        return result

    # -----------------------------------------------------------------

    def get_errormaplist_for_filters(self, filters):

        """
        This function ...
        :param filters: 
        :return: 
        """

        # Doesn't work when error frames are in the same FITS file as the primar frames!!
        #return NamedFrameList.from_paths(**self.get_errormap_names_and_paths_for_filters(filters))

        #frame_paths = self.get_frame_names_and_paths_for_filters(filters)
        #errormaps = {name: Frame.from_file(frame_paths[name], plane="errors") for name in frame_paths}

        # Create dict of error maps
        errormaps = {name: self.get_errormap(name) for name in self.get_names_for_filters(filters)}
        return NamedFrameList(**errormaps)

    # -----------------------------------------------------------------

    def has_frame_for_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        # Get the name
        name = self.get_name_for_filter(fltr)

        if name is None: return False
        else: return True

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, fltr, masked=True, mask_value=0.0):

        """
        This function ...
        :param fltr:
        :param masked:
        :param mask_value:
        :return:
        """

        # Get the name
        name = self.get_name_for_filter(fltr)

        # No frame found for this filter
        if name is None: return None

        # Return the frame
        return self.get_frame(name, masked=masked, mask_value=mask_value)

    # -----------------------------------------------------------------

    def has_errormap_for_filter(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        # Get the name
        name = self.get_name_for_filter(fltr)

        # If name, return False
        if name is None: return False

        # Check whether error map is present
        return name in self.error_paths

    # -----------------------------------------------------------------

    def get_errormap_for_filter(self, fltr, masked=True, mask_value=0.0):

        """
        THis function ...
        :param fltr: 
        :param masked: 
        :param mask_value:
        :return: 
        """

        # Get the name
        name = self.get_name_for_filter(fltr)

        # No frame is found for this filter
        if name is None: return None

        # Return the error map
        return self.get_errormap(name, masked=masked, mask_value=mask_value)

    # -----------------------------------------------------------------

    def get_coordinate_system(self, name):

        """
        This function ...
        :param name: 
        :return: 
        """

        return self.get_wcs(name)

    # -----------------------------------------------------------------

    def get_coordinate_system_for_filter(self, fltr, return_none=False):

        """
        This function ...
        :param filter:
        :param return_none:
        :return: 
        """

        return self.get_wcs_for_filter(fltr, return_none=return_none)

    # -----------------------------------------------------------------

    def get_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return CoordinateSystem.from_file(self.paths[name])

    # -----------------------------------------------------------------

    def get_wcs_for_filter(self, fltr, return_none=False):

        """
        This function ...
        :param fltr:
        :param return_none:
        :return:
        """

        name = self.get_name_for_filter(fltr, return_none=return_none)
        if name is None: return None
        return self.get_wcs(name)

    # -----------------------------------------------------------------

    def get_coordinate_box(self, name):

        """
        This function ...
        :param name: 
        :return: 
        """

        return self.get_wcs(name).bounding_box

    # -----------------------------------------------------------------

    def get_bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegionList()

        # Add the bounding boxes as sky rectangles
        for name in self.paths: boxes_region.append(self.get_coordinate_box(name))

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

    # -----------------------------------------------------------------

    def get_bounding_data(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    def get_overlap_box(self):

        """
        This function ...
        :return: 
        """

        min_ra = None
        max_ra = None
        min_dec = None
        max_dec = None

        # Loop over the frames
        for name in self.paths:

            # Get wcs
            wcs = self.get_coordinate_system(name)

            # Adjust
            if min_ra is None or wcs.min_ra > min_ra: min_ra = wcs.min_ra
            if max_ra is None or wcs.max_ra < max_ra: max_ra = wcs.max_ra
            if min_dec is None or wcs.min_dec > min_dec: min_dec = wcs.min_dec
            if max_dec is None or wcs.max_dec < max_dec: max_dec = wcs.max_dec

        center_ra = 0.5 * (min_ra + max_ra)
        center_dec = 0.5 * (min_dec + max_dec)

        # NO:
        #ra_radius = 0.5 * (max_ra - min_ra)
        dec_radius = 0.5 * (max_dec - min_dec)

        # BUT:
        ra_span = abs(coordinates.ra_distance(center_dec.to("deg").value, min_ra.to("deg").value, max_ra.to("deg").value))
        ra_radius = 0.5 * ra_span * u("deg")

        # Determine center and radius
        # Get center and radius of the new bounding box
        center = SkyCoordinate(center_ra, center_dec)
        radius = SkyStretch(ra_radius, dec_radius)

        # Return the bounding box
        return SkyRectangleRegion(center, radius)

    # -----------------------------------------------------------------

    def get_overlap_data(self):

        """
        This function ...
        :return: 
        """

        min_ra = None
        min_ra_name = None
        max_ra = None
        max_ra_name = None
        min_dec = None
        min_dec_name = None
        max_dec = None
        max_dec_name = None

        # Loop over the images
        for name in self.paths:

            # Get wcs
            wcs = self.get_coordinate_system(name)

            #
            if min_ra is None or wcs.min_ra > min_ra:
                min_ra = wcs.min_ra
                min_ra_name = name

            if max_ra is None or wcs.max_ra < max_ra:
                max_ra = wcs.max_ra
                max_ra_name = name

            if min_dec is None or wcs.min_dec > min_dec:
                min_dec = wcs.min_dec
                min_dec_name = name

            if max_dec is None or wcs.max_dec < max_dec:
                max_dec = wcs.max_dec
                max_dec_name = name

        # Create the data
        data = Map()
        data.ra = Map()
        data.dec = Map()
        data.ra.min = min_ra_name
        data.ra.max = max_ra_name
        data.dec.min = min_dec_name
        data.dec.max = max_dec_name

        center_ra = 0.5 * (min_ra + max_ra)
        center_dec = 0.5 * (min_dec + max_dec)

        # NO:
        #ra_radius = 0.5 * (max_ra - min_ra)
        dec_radius = 0.5 * (max_dec - min_dec)

        # BUT:
        ra_span = abs(coordinates.ra_distance(center_dec.to("deg").value, min_ra.to("deg").value, max_ra.to("deg").value))
        ra_radius = 0.5 * ra_span * u("deg")

        # Determine center and radius
        # Get center and radius of the new bounding box
        center = SkyCoordinate(center_ra, center_dec)
        radius = SkyStretch(ra_radius, dec_radius)

        box = SkyRectangleRegion(center, radius)

        # Return the box, and the data
        return box, data

    # -----------------------------------------------------------------

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.paths:

            wcs = self.get_wcs(name)
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.paths:

            wcs = self.get_wcs(name)
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def median_pixelscale(self):

        """
        This function ...
        :return:
        """

        name = self.median_pixelscale_name
        return self.get_pixelscale(name)

    # -----------------------------------------------------------------

    @property
    def pixelscale_range(self):
        return QuantityRange(self.min_pixelscale, self.max_pixelscale)

    # -----------------------------------------------------------------

    @property
    def min_pixelscale_wcs(self):
        return self.get_wcs(self.min_pixelscale_name)

    # -----------------------------------------------------------------

    @property
    def max_pixelscale_wcs(self):
        return self.get_wcs(self.max_pixelscale_name)

    # -----------------------------------------------------------------

    @property
    def median_pixelscale_wcs(self):
        name = self.median_pixelscale_name
        return self.get_wcs(name)

    # -----------------------------------------------------------------

    @property
    def min_pixelscale_name(self):

        """
        This property ...
        :return:
        """

        pixelscale = None
        frame_name = None

        # Loop over the images
        for name in self.paths:

            wcs = self.get_wcs(name)
            if pixelscale is None or wcs.average_pixelscale < pixelscale:
                pixelscale = wcs.average_pixelscale
                frame_name = name

        # Return the name for the frame with the minimum pixelscale
        return frame_name

    # -----------------------------------------------------------------

    @property
    def max_pixelscale_name(self):

        """
        This function ...
        :return:
        """

        pixelscale = None
        frame_name = None

        # Loop over the images
        for name in self.paths:

            wcs = self.get_wcs(name)
            if pixelscale is None or wcs.average_pixelscale > pixelscale:
                pixelscale = wcs.average_pixelscale
                frame_name = name

        # Return the name for the frame with the maximum pixelscale
        return frame_name

    # -----------------------------------------------------------------

    @property
    def median_pixelscale_name(self):

        """
        This function ...
        :return:
        """

        names = self.names
        pixelscales = [self.get_pixelscale(name) for name in names]

        # Determine which frame contains the median pixelscale
        scalar_pixelscales = [pixelscale.average.to("arcsec").value for pixelscale in pixelscales]
        median_index = np.argsort(scalar_pixelscales)[len(pixelscales) // 2]

        # Return
        return names[median_index]

    # -----------------------------------------------------------------

    def get_fwhm(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        header = self.get_header(name)
        header_fwhm = headers.get_fwhm(header)
        return header_fwhm

    # -----------------------------------------------------------------

    @property
    def min_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        for name in self.paths:

            # Get the FWHM
            header_fwhm = self.get_fwhm(name)

            if fwhm is None or header_fwhm < fwhm: fwhm = header_fwhm

        # Return the minimum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        for name in self.paths:

            # Get the FWHM
            header_fwhm = self.get_fwhm(name)

            if fwhm is None or header_fwhm > fwhm: fwhm = header_fwhm

        # Return the minimum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def fwhm_range(self):

        """
        This function ...
        :return:
        """

        return QuantityRange(self.min_fwhm, self.max_fwhm)

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):

        """
        This function ...
        :return:
        """

        wavelength = None

        for name in self.paths:

            # Get the wavelength
            header = self.get_header(name)
            header_wavelength = headers.get_filter(name, header).pivot

            if wavelength is None or header_wavelength < wavelength: wavelength = header_wavelength

        # Return the minimum wavelength
        return wavelength

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):

        """
        This function ...
        :return:
        """

        wavelength = None

        for name in self.paths:

            # Get the wavelength
            header = self.get_header(name)
            header_wavelength = headers.get_filter(name, header).pivot

            if wavelength is None or header_wavelength > wavelength: wavelength = header_wavelength

        # Return the maximum wavelength
        return wavelength

    # -----------------------------------------------------------------

    @property
    def wavelength_range(self):

        """"
        This function ...
        """

        return QuantityRange(self.min_wavelength, self.max_wavelength)

    # -----------------------------------------------------------------

    @property
    def largest_wcs_name(self):

        """
        Thisf unction ...
        :return:
        """

        names = self.names
        #coordinate_systems = [self.get_wcs(name) for name in names]
        areas = [self.get_wcs(name).area.to("sr").value for name in names]
        index = np.argmax(areas)
        return names[index]

    # -----------------------------------------------------------------

    @property
    def largest_wcs_below_median_pixelscale_name(self):

        """
        This function ...
        :return:
        """

        names = self.names

        pixelscales = [self.get_pixelscale(name) for name in names]

        # Determine which frame contains the median pixelscale
        scalar_pixelscales = [pixelscale.average.to("arcsec").value for pixelscale in pixelscales]
        median_index = np.argsort(scalar_pixelscales)[len(pixelscales) // 2]

        median_pixelscale = pixelscales[median_index]

        # Get pixelscales for determining largest wcs
        names_for_largest = []
        for index in range(len(names)):
            name = names[index]
            pixelscale = pixelscales[index]
            if pixelscale <= median_pixelscale: names_for_largest.append(name)

        # Determine the largest wcs
        areas = [self.get_wcs(name).area.to("sr").value for name in names_for_largest]
        index = np.argmax(areas)
        return names_for_largest[index]

    # -----------------------------------------------------------------

    def get_closest_pixelscale_name(self, pixelscale):

        """
        This function ...
        :param pixelscale:
        :return:
        """

        # Check input
        if isinstance(pixelscale, Pixelscale): pixelscale = pixelscale.average
        elif isinstance(pixelscale, Quantity): pass
        else: raise ValueError("Could not interpret pixelscale")

        # Get names and pixelscales
        names = self.names
        pixelscales = [self.get_pixelscale(name).average for name in names]

        # Find index
        index = sequences.find_closest_index(pixelscales, pixelscale)

        # Return
        return names[index]

    # -----------------------------------------------------------------

    def get_errormap_path(self, name):

        """
        This function ...
        :param name: 
        :return: 
        """

        return self.error_paths[name]

    # -----------------------------------------------------------------

    def get_errormap(self, name, masked=True, mask_value=0.0):

        """
        This function ...
        :param name:
        :param masked:
        :param mask_value:
        :return:
        """

        # Check if the name appears in the frames list
        if name not in self.paths: raise ValueError("Invalid name: " + name)

        # Get the errors frame
        if name in self.error_paths: errors = Frame.from_file(self.error_paths[name])
        elif "errors" in get_frame_names(self.paths[name]): errors = Frame.from_file(self.paths[name], plane="errors")
        else: return None

        # Check if the error frame has to be masked
        if masked:

            # Get the mask and set frame pixels to zero
            if name in self.mask_paths:
                mask = self.get_mask(name)
                errors[mask] = mask_value
            #else: log.warning("No mask available for " + name + " frame")

        # Set the name
        errors.name = name

        # Return the errors frame
        return errors

    # -----------------------------------------------------------------

    def get_relative_errors(self, name, masked=True, mask_value=0.0):

        """
        This function ...
        :param name:
        :param masked:
        :param mask_value:
        :return:
        """

        # Check if the name appears in the frames list
        if name not in self.paths: raise ValueError("Invalid name: " + name)

        # Get the frame
        frame = self.get_frame(name)

        # Get the errors
        errors = self.get_errormap(name)
        if errors is None: return None

        # Calculate the relative error map
        rel_errors = errors / frame

        # Check if the error frame has to be masked
        if masked:

            # Get the mask and set frame pixels to zero
            if name in self.mask_paths:
                mask = self.get_mask(name)
                rel_errors[mask] = mask_value
            #else: log.warning("No mask available for " + name + " frame")

        # Return the relative error map
        return rel_errors

    # -----------------------------------------------------------------

    def get_significance(self, name, levels=None, below_levels_value=float("nan")):

        """
        This function ...
        :param name:
        :param levels:
        :param below_levels_value:
        :return:
        """

        # Sort the level bins from small to large
        if levels is not None: levels = sorted(levels)

        # Get the frame
        frame = self.get_frame(name)

        # Get the errors
        errors = self.get_errormap(name)
        if errors is None: return None

        # If level bins are specified
        if levels is not None:

            # Create a frame full of nans
            significance = Frame.filled_like(frame, below_levels_value)

            # Loop over the levels
            for level in levels:
                significance[frame > level * errors] = level

        # No level bins, just calculate the exact significance level in each pixel
        else: significance = frame / errors

        # Return the significance map
        return significance

    # -----------------------------------------------------------------

    def get_significance_mask(self, name, level):

        """
        This function ...
        :param name:
        :param level:
        :return:
        """

        # Get significance map
        significance = self.get_significance(name)

        # Return the mask for a certain significance level
        return significance > level

    # -----------------------------------------------------------------

    def get_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Open the image and return it
        image = Image.from_file(self.paths[name])

        # Set the name
        image.name = name

        # Return the image
        return image

    # -----------------------------------------------------------------

    def get_image_plane(self, name, plane_name):

        """
        This function ...
        :param name:
        :param plane_name:
        :return:
        """

        # Open the requested frame
        frame = Frame.from_file(self.paths[name], plane=plane_name)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def get_image_mask(self, name, mask_name):

        """
        This function ...
        :param name: 
        :param mask_name
        :return: 
        """

        image = self.get_image(name)
        return image.masks[mask_name]

    # -----------------------------------------------------------------

    def get_image_masks(self, name, mask_names, strict=False):

        """
        This function ...
        :param name: 
        :param mask_names: 
        :param strict: 
        :return: 
        """

        present_mask_names = self.masks_in_image(name)

        masks = []
        for mask_name in mask_names:
            if mask_name not in present_mask_names:
                if strict:
                    raise ValueError("Mask '" + mask_name + "' not present in image '" + name + "'")
                else:
                    continue

            # Get the mask
            mask = self.get_image_mask(name, mask_name)

            # Add the mask
            masks.append(mask)

        # Return the list of masks
        return masks

    # -----------------------------------------------------------------

    def get_image_masks_union(self, name, mask_names, strict=False):

        """
        This function ...
        :param name: 
        :param mask_names: 
        :param strict: 
        :return: 
        """

        # Get the masks
        masks = self.get_image_masks(name, mask_names, strict=strict)

        # No masks
        if len(masks) == 0: return None

        # Create the intersection
        return union(*masks)

    # -----------------------------------------------------------------

    def get_image_masks_intersection(self, name, mask_names, strict=False):

        """
        This function ...
        :param name: 
        :param mask_names: 
        :param strict:
        :return: 
        """

        # Get the masks
        masks = self.get_image_masks(name, mask_names, strict=strict)

        # No masks
        if len(masks) == 0: return None

        # Create the intersection
        return intersection(*masks)

    # -----------------------------------------------------------------

    def planes_in_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return get_plane_names(self.paths[name])

    # -----------------------------------------------------------------

    def frames_in_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return get_frame_names(self.paths[name])

    # -----------------------------------------------------------------

    def masks_in_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return get_mask_names(self.paths[name])

    # -----------------------------------------------------------------

    def get_mask_path(self, name):

        """
        This function ...
        :param name: 
        :return: 
        """

        return self.mask_paths[name]

    # -----------------------------------------------------------------

    def get_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Check if the name appears in the frames list
        if name not in self.paths: raise ValueError("Invalid name: " + name)

        # Check if the name appears in the masks list
        if name not in self.mask_paths: raise ValueError("The " + name + " frame has no mask")

        # Otherwise, return the mask
        mask = Mask.from_file(self.mask_paths[name])
        return mask

    # -----------------------------------------------------------------

    def create_datacube(self, min_wavelength=None, max_wavelength=None, exclude=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :return:
        """

        # Get the frames
        frames = self.get_frames(min_wavelength=min_wavelength, max_wavelength=max_wavelength, exclude=exclude)

        # Debugging
        log.debug("Frames to be used for constructing the datacube:")
        for name in frames: log.debug(" - " + name)

        # Inform the user
        log.info("Determining which image will be used as the reference for rebinning all other images ...")

        # Get the name of the frame with the lowest spatial resolution
        reference_name = get_lowest_resolution_name(frames)

        # Inform the user
        log.info("The reference image used for the rebinning is the " + reference_name + " image")

        # Loop over all images
        for name in frames:

            # Don't rebin the reference image
            if name == reference_name: continue

            # Debugging
            log.debug("Rebinning the " + name + " frame to the pixel grid of the " + reference_name + " image ...")

            # Rebin this frame to the lower resolution pixel grid
            frames[name].rebin(frames[reference_name].wcs)

        # Create the datacube and return it
        # NOTE: the frames are sorted by wavelength by the DataCube constructor
        return DataCube.from_frames(frames.values())

    # -----------------------------------------------------------------

    def create_errorcube(self, min_wavelength=None, max_wavelength=None, exclude=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param exclude:
        :return:
        """

        # Get the error maps
        maps = self.get_errormaps(min_wavelength=min_wavelength, max_wavelength=max_wavelength, exclude=exclude)

        # Debugging
        log.debug("Error maps to be used for constructing the errorcube:")
        for name in maps: log.debug(" - " + name)

        # Inform the user
        log.info("Determining which error map will be used as the reference for rebinning all other error maps ...")

        # Get the name of the error map with the lowest spatial resolution
        reference_name = get_lowest_resolution_name(maps)

        # Inform the user
        log.info("The reference error map used for the rebinning is the " + reference_name + " error map")

        # Loop over all images
        for name in maps:

            # Don't rebin the reference map
            if name == reference_name: continue

            # Debugging
            log.debug("Rebinning the " + name + " error map to the pixel grid of the " + reference_name + " image ...")

            # Rebin this error map to the lower resolution pixel grid
            maps[name].rebin(maps[reference_name].wcs)

        # Create the datacube and return it
        # NOTE: the maps are sorted by wavelength by the DataCube constructor
        return DataCube.from_frames(maps.values())

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check if file exists
        if self.path is None: raise ValueError("The dataset file does not exist yet")

        # Save to the original path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the dataset to " + path + " ...")

        # Create the table
        table = DataSetTable()

        for name in self.paths:
            frame_path = self.paths[name]
            error_path = self.error_paths[name] if name in self.error_paths else None
            mask_path = self.mask_paths[name] if name in self.mask_paths else None
            table.add_entry(name, frame_path, error_path, mask_path)

        # Save the table
        table.saveto(path)

        # Update the path
        self.path = path

    # -----------------------------------------------------------------

    def save_as_fits(self, path):

        """
        This function saves the dataset as a single FITS files with multiple HDUs
        :return:
        """

        # Inform the user
        log.info("Saving the dataset to " + path + " ...")

        # Create the HDUList
        hdulist = fits.HDUList()

        # Add the HDUs
        for name in self.paths:

            # Open the frame and error map
            frame = self.get_frame(name)
            errors = self.get_errormap(name)

            # Create the header
            header = frame.header

            # The HDU data
            if errors is None: data = frame._data
            else:

                # Set plane names and set data
                header["PLANE0"] = "primary [frame]"
                header["PLANE1"] = "errors [frame]"
                header["NAXIS"] = 3
                header["NAXIS3"] = 2
                data = np.array([frame._data, errors._data])

            # Create an imageHDU
            hdu = fits.ImageHDU(data, header=header, name=name)

            # Add the HDU to the HDUList
            hdulist.append(hdu)

        # Write the HDU list
        hdulist.writeto(path)

# -----------------------------------------------------------------

from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter

# -----------------------------------------------------------------

class DataSetCreator(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DataSetCreator, self).__init__(*args, **kwargs)

        # The list of image paths
        self.image_paths = None

        # The list of error map paths
        self.error_paths = None

        # The data set
        self.dataset = DataSet()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 2. Load the image paths
        if self.image_paths is None: self.load_paths()

        # 3. Create the dataset
        self.create()

        # 4. Write the dataset
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataSetCreator, self).setup(**kwargs)

        # Get the image and error map paths
        self.image_paths = kwargs.pop("image_paths", None)
        self.error_paths = kwargs.pop("error_paths", None)

    # -----------------------------------------------------------------

    def load_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the image paths ...")

        # Interactive or from file
        if self.config.interactive: self.load_interactive()
        else: self.load_from_cwd()

    # -----------------------------------------------------------------

    def load_interactive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading interactive prompt for image paths ...")

        # The configuration setter
        setter = InteractiveConfigurationSetter("datasetcreator_images", add_cwd=False, add_logging=False)

        definition = ConfigurationDefinition()
        definition.add_required("image_paths", "filepath_list", "paths to the images")

        # Create the configuration and get the paths
        config = setter.run(definition, prompt_optional=False)
        self.image_paths = config.image_paths

        self.error_paths = dict()

        for path in self.image_paths:

            name = fs.strip_extension(fs.name(path))

            definition = ConfigurationDefinition()
            definition.add_optional("error_path", "file_path", "path to the error map for " + name)

            setter = InteractiveConfigurationSetter("datasetcreator_errors", add_cwd=False, add_logging=False)

            # Get the error path
            config = setter.run(definition, prompt_optional=True)
            error_path = config.error_path

            if error_path is not None: self.error_paths[name] = error_path

    # -----------------------------------------------------------------

    def load_from_cwd(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading image paths from '" + self.config.path + "' ...")

        self.image_paths = []
        self.error_paths = dict()

        # Loop over the FITS files in the current directory
        for image_path, image_name in fs.files_in_path(self.config.path, extension="fits", contains=self.config.contains,
                                                       not_contains=self.config.not_contains, returns=["path", "name"],
                                                       recursive=self.config.recursive):

            # Skip error maps
            if self.config.error_suffix is not None and image_name.endswith(self.config.error_suffix): continue

            # Add the image path
            self.image_paths.append(image_path)

            # Look for the error map, if error suffix is specified
            if self.config.error_suffix is not None:

                error_path = fs.join(fs.directory_of(image_path), image_name + self.config.error_suffix + ".fits")
                if fs.is_file(error_path): self.error_paths[image_name] = error_path

    # -----------------------------------------------------------------

    def create(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dataset ...")

        # Loop over the image paths
        for path in self.image_paths:

            # Determine the filename
            filename = fs.strip_extension(fs.name(path))

            # Open the image frame
            frame = Frame.from_file(path)

            # Determine the preparation name
            if frame.filter is not None: name = str(frame.filter)
            else: name = filename

            # Add the path to the dataset
            self.dataset.add_path(name, path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.config.path, "dataset.dat")
        self.dataset.saveto(path)

# -----------------------------------------------------------------

StaticDataSet = create_lazified_class(DataSet, "StaticDataSet")

# -----------------------------------------------------------------

def get_lowest_resolution_name(frames):

    """
    This function ...
    :param frames:
    :return:
    """

    # Get the name of the image with the lowest resolution
    lowest_resolution = None
    reference_name = None
    for name in frames:
        pixelscale = frames[name].pixelscale.average
        if lowest_resolution is None or pixelscale > lowest_resolution:
            lowest_resolution = pixelscale
            reference_name = name

    return reference_name

# -----------------------------------------------------------------
