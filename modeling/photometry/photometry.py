#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.photometry Contains the PhotoMeter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant AstroMagic classes and modules
from ...magic.core import Image
from ...magic.tools import headers

# Import the relevant PTS classes and modules
from ..core import ModelingComponent
from ...core.tools import filesystem, tables
from ...core.tools.logging import log

# -----------------------------------------------------------------

class PhotoMeter(ModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PhotoMeter, self).__init__(config)

        # The list of images
        self.images = []

        # The photometry table
        self.table = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new PhotoMeter instance
        photometer = cls(arguments.config)

        # Set the input and output path
        photometer.config.path = arguments.path
        photometer.config.input_path = os.path.join(arguments.path, "prep")
        photometer.config.output_path = os.path.join(arguments.path, "phot")

        # A single image can be specified so the photometry is only calculated for that image
        photometer.config.single_image = arguments.image

        # Return the new instance
        return photometer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the prepared images
        self.load_images()

        # 3. Do the photometry
        self.do_photometry()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PhotoMeter, self).__init__()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in filesystem.directories_in_path(self.prep_path, names=True):

            # Determine the filter
            filter_name = directory_name
            if filter_name == "Ha": filter_name = "H alpha"
            filter = headers.get_filter(filter_name)

            # Look for a file called 'result.fits'
            image_path = os.path.join(directory_path, "result.fits")
            if not filesystem.is_file(image_path):

                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Open the prepared image
            image = Image.from_file(image_path)

            # Set the image name
            image.name = directory_name

            # Set the filter
            image.filter = filter

            # Add the image to the list
            self.images.append(image)

    # -----------------------------------------------------------------

    def do_photometry(self):

        """
        This function ...
        :return:
        """

        wavelength_column = []
        flux_column = []

        for image in self.images:

            wavelength_column.append(image.wavelength)
            flux_column.append(np.sum(image.frames.primary))

        data = [wavelength_column, flux_column]
        names = ["Wavelength", "Flux"]
        self.table = tables.new(data, names)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        self.write_table()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the output file
        path = self.full_output_path("fluxes.dat")

        # Write the photometry table
        tables.write(self.table, path)

# -----------------------------------------------------------------
