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
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .frame import Frame
from ...core.tools.logging import log
from .datacube import DataCube

# -----------------------------------------------------------------

class DataSet(object):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # The paths to the images
        self.paths = OrderedDict()

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

    def get_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Open the frame and return it
        return Frame.from_file(self.paths[name])

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
