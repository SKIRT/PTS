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
import numpy as np
from collections import OrderedDict

# Import astronomical modules
from astropy.io import fits

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .frame import Frame
from .io import get_frame_names, get_mask_names, get_plane_names
from ...core.tools.logging import log
from .datacube import DataCube
from .image import Image
from .mask import Mask
from ..basics.coordinatesystem import CoordinateSystem
from ...core.basics.configurable import Configurable
from ...core.tools import tables
from ..basics.skyregion import SkyRegion
from ..tools import headers

# -----------------------------------------------------------------

class DataSet(object):

    """
    This class...
    """

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
    def from_directory(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Set the paths
        paths, names = fs.files_in_path(path, extension="fits", returns=["path", "name"])

        # Create a new dataset instance
        dataset = cls()

        # Add the paths
        for path, name in zip(paths, names):

            # Add the image path
            dataset.add_path(name, path)

        # Return the dataset instance
        return dataset

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Construct a new dataset
        dataset = cls()

        # Load the table
        table = tables.from_file(path)

        # Loop over the entries in the table
        for i in range(len(table)):

            # Get the paths
            name = table["Name"][i]
            path = table["Path"][i]
            error_path = table["Error path"][i] if not table["Error path"].mask[i] else None
            mask_path = table["Mask path"][i] if not table["Mask path"].mask[i] else None

            # Add the paths to the dataset
            dataset.add_path(name, path)
            if error_path is not None: dataset.add_error_path(name, error_path)
            if mask_path is not None: dataset.add_mask_path(name, mask_path)

        # Return the dataset
        return dataset

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.paths.keys()

    # -----------------------------------------------------------------

    def add_path(self, name, path):

        """
        This function ...
        :param name:
        :param path:
        :return:
        """

        # Check if the file exists
        if not fs.is_file(path): raise IOError("File doesn't exist: '" + path + "'")

        # Check if already such a name
        if name in self.paths: raise ValueError("Already a path in the dataset with the name " + name)

        # Add the path
        self.paths[name] = path

    # -----------------------------------------------------------------

    def add_error_path(self, name, path):

        """
        This function ...
        :param name:
        :param path:
        :return:
        """

        # Check if the file exists
        if not fs.is_file(path): raise IOError("File doesn't exist: '" + path + "'")

        # Check if corresponding frame exists
        if name not in self.paths: raise ValueError("Corresponding image with name " + name + " has not been added")

        # Check if already such a name
        if name in self.error_paths: raise ValueError("Already an error path in the dataset with the name " + name)

        # Add the path
        self.error_paths[name] = path

    # -----------------------------------------------------------------

    def add_mask_path(self, name, path):

        """
        This function ...
        :param name:
        :param path:
        :return:
        """

        # Check if the file exists
        if not fs.is_file(path): raise IOError("File does not exist: '" + path + "'")

        # Check if the corresponding frame exists
        if name not in self.paths: raise ValueError("Corresponding image with name " + name + " has not been added")

        # Check if already such a name
        if name in self.mask_paths: raise ValueError("Already a mask path in the dataset with the name " + name)

        # Add the path
        self.mask_paths[name] = path

    # -----------------------------------------------------------------

    def get_filters(self):

        """
        This function ...
        :return:
        """

        # Initialize
        fltrs = dict()

        # Loop over the images
        for name in self.paths: fltrs[name] = self.get_frame(name, masked=False).filter

        # Return the dictionary with the filters
        return fltrs

    # -----------------------------------------------------------------

    def get_filter(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_frame(name).filter

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
            else: log.warning("No mask available for " + name + " frame")

        # Set the name
        frame.name = name

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, fltr, masked=True):

        """
        This function ...
        :param fltr:
        :param masked:
        :return:
        """

        filter_string = str(fltr)

        # Get the dictionary of filters
        fltrs = self.get_filters()

        # Loop over the filters dictionary
        for name in fltrs:

            # Get the frame
            frame = self.get_frame(name, masked=masked)

            # Check the filter
            if str(frame.filter) == filter_string: return frame

        # No frame found
        return None

    # -----------------------------------------------------------------

    def get_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return CoordinateSystem.from_file(self.paths[name])

    # -----------------------------------------------------------------

    def get_bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegion()

        # Add the bounding boxes as sky rectangles
        for name in self.paths: boxes_region.append(self.get_wcs(name).bounding_box)

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

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
    def min_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        for name in self.paths:

            # Get the FWHM
            header = self.get_header(name)
            header_fwhm = headers.get_fwhm(header)

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
            header = self.get_header(name)
            header_fwhm = headers.get_fwhm(header)

            if fwhm is None or header_fwhm > fwhm: fwhm = header_fwhm

        # Return the minimum FWHM
        return fwhm

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

    def get_errors(self, name, masked=True, mask_value=0.0):

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
            else: log.warning("No mask available for " + name + " frame")

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
        errors = self.get_errors(name)
        if errors is None: return None

        # Calculate the relative error map
        rel_errors = errors / frame

        # Check if the error frame has to be masked
        if masked:

            # Get the mask and set frame pixels to zero
            if name in self.mask_paths:
                mask = self.get_mask(name)
                rel_errors[mask] = mask_value
            else: log.warning("No mask available for " + name + " frame")

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
        errors = self.get_errors(name)
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

        # Inform the user
        log.info("Determining which image will be used as the reference for rebinning all other images ...")

        # Make sure exclude is a list
        if isinstance(exclude, basestring): exclude = [exclude]

        # The image frames
        frames = dict()

        # Open the image frames
        for name in self.paths:

            # Skip if name is in the exclude list
            if exclude is not None and name in exclude: continue

            # Open the frame
            frame = self.get_frame(name)

            # Skip images of wavelength smaller than the minimum or greater than the maximum
            if min_wavelength is not None and frame.wavelength < min_wavelength: continue
            if max_wavelength is not None and frame.wavelength > max_wavelength: continue

            # Add the frame
            frames[name] = frame

        # Get the name of the image with the lowest resolution
        lowest_resolution = None
        reference_name = None
        for name in frames:
            pixelscale = frames[name].pixelscale.average
            if lowest_resolution is None or pixelscale > lowest_resolution:
                lowest_resolution = pixelscale
                reference_name = name

        # Inform the user
        log.info("The reference image used for the rebinning is the " + reference_name + " image")

        # Loop over all images
        for name in frames:

            # Don't rebin the reference image
            if name == reference_name: continue

            # Rebin this frame to the lower resolution pixel grid
            frames[name].rebin(frames[reference_name].wcs)

        # Create the datacube and return it
        return DataCube.from_frames(frames.values())

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

        # Create a table
        column_names = ["Name", "Path", "Error path", "Mask path"]

        # Set the entries
        names = []
        paths = []
        error_paths = []
        mask_paths = []

        for name in self.paths:

            names.append(name)
            paths.append(self.paths[name])
            error_paths.append(self.error_paths[name] if name in self.error_paths else None)
            mask_paths.append(self.mask_paths[name] if name in self.mask_paths else None)

        # Construct the table
        data = [names, paths, error_paths, mask_paths]
        table = tables.new(data, column_names)

        # Save the table
        tables.write(table, path)

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
            errors = self.get_errors(name)

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

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(DataSetCreator, self).__init__(config)

        # The list of image paths
        self.image_paths = None

        # The list of error map paths
        self.error_paths = None

        # The data set
        self.dataset = DataSet()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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
        setter = InteractiveConfigurationSetter("datasetcreator_images")

        definition = ConfigurationDefinition()
        definition.add_required("image_paths", "filepath_list", "paths to the images")

        # Create the configuration and get the paths
        config = setter.run(definition)
        self.image_paths = config.image_paths

        self.error_paths = dict()

        for path in self.image_paths:

            name = fs.strip_extension(fs.name(path))

            definition = ConfigurationDefinition()
            definition.add_optional("error_path", "filepath", "path to the error map for " + name)

            setter = InteractiveConfigurationSetter("datasetcreator_errors")

            # Get the error path
            config = setter.run(definition)
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
