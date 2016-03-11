#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.initialization Contains the DataInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .component import PreparationComponent
from ...magic.sources.finder import SourceFinder
from ...core.tools import filesystem, tables
from ...core.tools.logging import log
from ...magic.misc.imageimporter import ImageImporter
from ...magic.catalog.importer import CatalogImporter

# -----------------------------------------------------------------

class DataInitializer(PreparationComponent):
    
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
        super(DataInitializer, self).__init__(config)

        # -- Attributes --

        # The list of images
        self.images = []

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the images
        self.load_images()

        # 3. Identify the sources in the images
        self.find_sources()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # -- Children --
        self.add_child("importer", ImageImporter, self.config.importation)
        self.add_child("catalog_importer", CatalogImporter, self.config.catalogs)
        self.add_child("source_finder", SourceFinder, self.config.sources)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(DataInitializer, self).setup()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Load image information
        info_path = filesystem.join(self.data_path, "info.dat")
        info_table = tables.from_file(info_path)

        # Loop over all files found in the data directory
        for image_path, image_name in filesystem.files_in_path(self.data_path, extension="fits", not_contains="error", returns="both"):

            # If only a single image must be prepared, check if this image matches the specified image name
            if self.config.single_image is not None and image_name != self.config.single_image: continue

            # Determine the output path for this image
            image_output_path = filesystem.join(self.prep_path, image_name)

            # Check whether this image already has an initialized image
            final_path = filesystem.join(image_output_path, "initialized.fits")
            if filesystem.is_file(final_path): continue

            # Get the corresponding index in the information table
            info_index = tables.find_index(info_table, image_name)
            if info_index is None:
                # TEMP: skip image if not defined in table !!
                log.warning("No information about " + image_name + ": skipping")
                continue

            # Get image properties such as the unit and the FWHM of the PSF
            unit = Unit(info_table["Unit"][info_index]) if not info_table["Unit"].mask[info_index] else None
            fwhm = info_table["FWHM"][info_index] * Unit(info_table["FWHM unit"][info_index]) if not info_table["FWHM"].mask[info_index] else None

            # Set the FWHM of the reference image so that we can set it for the convolution kernel
            #if image_name == self.config.reference_image: self.reference_fwhm = fwhm

            # Set the path to the region of bad pixels
            bad_region_path = filesystem.join(self.data_path, "bad", image_name + ".reg")
            if not filesystem.is_file(bad_region_path): bad_region_path = None

            # Import the image
            importer = ImageImporter()
            importer.run(image_path, bad_region_path, unit=unit, fwhm=fwhm)

            # Add the image that has to be processed to the list
            self.images.append(importer.image)

            # Clear the image importer
            self.importer.clear()

        #if self.reference_fwhm is None:

            #reference_index = tables.find_index(info_table, self.config.reference_image)
            #self.reference_fwhm = info_table["FWHM"][reference_index] * Unit(info_table["FWHM unit"][reference_index])

        #assert self.reference_fwhm is not None

    # -----------------------------------------------------------------

    def get_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Importing the galactic and stellar catalogs ...")

        # Find the coordinate box of the largest image

        # Run the catalog importer
        self.catalog_importer.run(self.image.frames.primary)

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Run the source finder on each image

        # Loop over the different images
        for image in self.images:

            # Run the source finder on this image
            self.source_finder.run(image, self.catalog_importer.galactic_catalog, self.catalog_importer.stellar_catalog)

    # -----------------------------------------------------------------

    def save_images(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
