#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.fluxes Contains the ObservedImageMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools.logging import log
from ..tools import filesystem
from ..basics.filter import Filter
from ...magic.core.image import Image
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

class ObservedImageMaker(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(ObservedImageMaker, self).__init__()

        # -- Attributes --

        # The paths to the FITS files produced by SKIRT
        self.fits_paths = None

        # The wavelengths of the simulation
        self.wavelengths = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # The filters for which the images should be created
        self.filters = None

        # The dictionary containing the images
        self.images = dict()

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None, filter_names=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :param filter_names:
        :return:
        """

        # Obtain the paths to the 'total' FITS files created by the simulation
        self.fits_paths = simulation.totalfitspaths()

        # Get the list of wavelengths for the simulation
        self.wavelengths = simulation.wavelengths()

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Create the filters
        self.create_filters()

        # Make the observed images
        self.make_images()

        # Write the results
        if output_path is not None: self.write(output_path)

    # -----------------------------------------------------------------

    def create_filters(self):

        """
        This function ...
        :return:
        """

        # Initialize the list
        self.filters = []

        # Loop over the different filter names
        for filter_name in self.filter_names:

            # Create the filter
            filter = Filter.from_string(filter_name)

            # Add the filter to the list
            self.filters.append(filter)

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the observed images ...")

        # Loop over the different simulated images
        for path in self.fits_paths:

            # Load the simulated image
            datacube = Image.from_file(path)

            # Convert the datacube to a numpy array where wavelength is the third dimension
            fluxdensities = datacube.as_array()

            # Loop over the different filters
            for filter in self.filters:

                # Calculate the observed image
                data = filter.convolve(self.wavelengths, fluxdensities)
                image = Frame(data)

                # Add the observed image to the dictionary
                self.images[filter.name] = image

    # -----------------------------------------------------------------

    def write(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Loop over the different images
        for image_name in self.images:

            # Determine the path to the output image
            path = filesystem.join(output_path, image_name + ".fits")

            # Save the image
            self.images[image_name].save(path)

# -----------------------------------------------------------------
