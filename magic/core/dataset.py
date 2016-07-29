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
from .io import get_frame_names
from ...core.tools.logging import log
from .datacube import DataCube
from .image import Image

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

    def get_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Open the frame and return it
        return Frame.from_file(self.paths[name])

    # -----------------------------------------------------------------

    def get_errors(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Check if the name appears in the frames list
        if name not in self.paths: raise ValueError("Invalid name: " + name)

        # Get the errors frame
        if name in self.error_paths: return Frame.from_file(self.error_paths[name])
        elif "errors" in get_frame_names(self.paths[name]): return Frame.from_file(self.paths[name], plane="errors")
        else: return None

    # -----------------------------------------------------------------

    def get_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Open the image and return it
        return Image.from_file(self.paths[name])

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

        # Update the path
        self.path = path

# -----------------------------------------------------------------
