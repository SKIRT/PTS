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
from ...core.tools import filesystem
from ...core.tools.logging import log
from ...magic.misc.imageimporter import ImageImporter
from ...magic.catalog.importer import CatalogImporter
from ...magic.core.image import Image
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

# Reference: Common-Resolution Convolution Kernels for Space- and Ground-Based Telescopes (G. Aniano et. al)
fwhms = {"GALEX FUV": 4.48 * Unit("arcsec"),
         "GALEX NUV": 5.05 * Unit("arcsec"),
         "Mosaic Halpha": 2.0 * Unit("arcsec"),
         "IRAC I1": 1.90 * Unit("arcsec"),
         "IRAC I2": 1.81 * Unit("arcsec"),
         "IRAC I3": 2.11 * Unit("arcsec"),
         "IRAC I4": 2.82 * Unit("arcsec"),
         "WISE W1": 5.79 * Unit("arcsec"),
         "WISE W2": 6.37 * Unit("arcsec"),
         "WISE W3": 6.60 * Unit("arcsec"),
         "WISE W4": 11.89 * Unit("arcsec"),
         "MIPS 24mu": 6.43 * Unit("arcsec"),
         "MIPS 70mu": 18.74 * Unit("arcsec"),
         "MIPS 160mu": 38.78 * Unit("arcsec"),
         "Pacs blue": 5.67 * Unit("arcsec"),
         "Pacs green": 7.04 * Unit("arcsec"),
         "Pacs red": 11.18 * Unit("arcsec"),
         "SPIRE PSW_ext": 18.15 * Unit("arcsec"),
         "SPIRE PMW_ext": 24.88 * Unit("arcsec"),
         "SPIRE PLW_ext": 36.09 * Unit("arcsec")}

# -----------------------------------------------------------------

# The total H-alpha flux (reference: FAR-ULTRAVIOLET AND Ha IMAGING OF NEARBY SPIRAL GALAXIES: THE OB STELLAR,
# POPULATION IN THE DIFFUSE IONIZED GAS (Hoopes et. al 2001)
halpha_flux = 7.8e40 * Unit("erg/s")

# -----------------------------------------------------------------

reference_image = "Pacs red"

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

        # The paths of the images to be processed
        self.paths = []

        # The reference image frame
        self.reference = None

        # The galactic and stellar catalogs
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new DataInitializer instance
        initializer = cls()

        # Return the data initializer
        return initializer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the images
        self.load_reference_image()

        # 3. Get the galactic and stellar catalogs
        self.get_catalogs()

        # 4. Check which images should be processed
        self.check_images()

        # 4. Process the images (identify sources, create error frames)
        self.process_images()

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

    def load_reference_image(self):

        """
        This function ...
        :return:
        """

        # Get the path to the reference image
        reference_path = self.original_paths[reference_image]
        self.reference = Frame.from_file(reference_path)

    # -----------------------------------------------------------------

    def get_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Importing the galactic and stellar catalogs ...")

        # Run the catalog importer
        self.catalog_importer.run(self.reference)

        # Get the galactic and stellar catalogs
        self.galactic_catalog = self.catalog_importer.galactic_catalog
        self.stellar_catalog = self.catalog_importer.stellar_catalog

    # -----------------------------------------------------------------

    def check_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the input images ...")

        # Loop over all subdirectories of the data directory
        for path in filesystem.directories_in_path(self.data_path, not_contains="bad", returns="path"):

            # Debugging
            log.debug("Opening " + path + " ...")

            # Loop over all FITS files found in the current subdirectory
            for image_path in filesystem.files_in_path(path, extension="fits", not_contains="_Error", returns="path"):

                # Ignore the Planck data (for now)
                if "Planck" in image_path: continue

                # Name
                image_name = filesystem.strip_extension(filesystem.name(image_path))

                # Debugging
                log.debug("Checking " + image_path + " ...")

                # Determine the output path for this image
                output_path = self.prep_paths[image_name]

                # Check whether this image already has an initialized image
                final_path = filesystem.join(output_path, "initialized.fits")
                if filesystem.is_file(final_path): continue

                # Add the image path to the list
                self.paths.append(image_path)

    # -----------------------------------------------------------------

    def process_images(self):

        """
        This function ...
        :return:
        """

        # Loop over all image paths
        for image_path in self.paths:

            # Name
            image_name = filesystem.strip_extension(filesystem.name(image_path))

            # Determine the name used to identify this image for the preparation routines
            prep_name = self.prep_names[image_name]

            # Determine the output path for this image
            output_path = self.prep_paths[image_name]

            # Set the path to the region of bad pixels
            bad_region_path = filesystem.join(self.data_path, "bad", prep_name + ".reg")
            if not filesystem.is_file(bad_region_path): bad_region_path = None

            # Set the FWHM if the instrument has a fixed PSF
            if prep_name in fwhms: fwhm = fwhms[prep_name]
            else: fwhm = None

            # Debugging
            log.debug("Loading image " + image_path + " as " + prep_name + " ...")

            # Import the image
            importer = ImageImporter()
            importer.run(image_path, bad_region_path, fwhm=fwhm)

            # Get the imported image
            image = importer.image

            # Get the mask of bad pixels
            bad_mask = image.masks.bad if "bad" in image.masks else None

            # -----------------------------------------------------------------

            # Don't look for stars in the Halpha image
            if "Halpha" in prep_name: self.source_finder.config.find_stars = False
            else: self.source_finder.config.find_stars = True # still up to the SourceFinder to decide whether stars should be found (based on the filter)

            # Run the source finder on this image
            self.source_finder.run(image.frames.primary, self.galactic_catalog, self.stellar_catalog, bad_mask=bad_mask)

            # Determine the path to the "sources" directory within the output path for this image
            sources_output_path = filesystem.join(output_path, "sources")
            filesystem.create_directory(sources_output_path)

            # Save the galaxy region
            galaxy_region = self.source_finder.galaxy_region
            path = filesystem.join(sources_output_path, "galaxies.reg")
            galaxy_region.save(path)

            # Save the star region
            star_region = self.source_finder.star_region
            path = filesystem.join(sources_output_path, "stars.reg")
            if star_region is not None: star_region.save(path)

            # Save the saturation region
            saturation_region = self.source_finder.saturation_region
            path = filesystem.join(sources_output_path, "saturation.reg")
            if saturation_region is not None: saturation_region.save(path)

            # Save the region of other sources
            other_region = self.source_finder.other_region
            path = filesystem.join(sources_output_path, "other_sources.reg")
            other_region.save(path)

            # -----------------------------------------------------------------

            # Create an image with the segmentation maps
            segments = Image("segments")

            # Add the segmentation map of the galaxies
            segments.add_frame(self.source_finder.galaxy_segments, "galaxies")

            # Add the segmentation map of the saturated stars
            if self.source_finder.star_segments is not None: segments.add_frame(self.source_finder.star_segments, "stars")

            # Add the segmentation map of the other sources
            segments.add_frame(self.source_finder.other_segments, "other_sources")

            # Save the FITS file with the segmentation maps
            path = filesystem.join(sources_output_path, "segments.fits")
            segments.save(path)

            # -----------------------------------------------------------------

            # Set the FWHM of the image
            if image.fwhm is None: image.fwhm = self.source_finder.fwhm

            # -----------------------------------------------------------------

            # Determine the path to the initialized image
            path = filesystem.join(self.prep_paths[image.name], "initialized.fits")

            # Save the image
            image.save(path)

            # -----------------------------------------------------------------

            # Clear the source finder
            self.source_finder.clear()

            # Clear the image importer
            self.importer.clear()

# -----------------------------------------------------------------
