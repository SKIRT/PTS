#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.extraction Contains the Extractor class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant AstroMagic classes and modules
from ..core import Frame
from ..basics import Mask, Region
from .galaxyextraction import GalaxyExtractor
from .starextraction import StarExtractor
from .trainedextractor import TrainedExtractor
from ..catalog import CatalogBuilder, CatalogImporter, CatalogSynchronizer
from ..tools import wavelengths

# Import the relevant PTS classes and modules
from ...core.tools import filesystem, tables
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log

# -----------------------------------------------------------------

class Extractor(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(Extractor, self).__init__(config, "magic")

        # -- Attributes --

        # The image on which to perform the extraction
        self.image = None

        # The mask covering pixels that should be ignored throughout the entire extraction procedure
        self.special_mask = None
        self.ignore_mask = None

        # The output mask
        self.mask = None

        # The name of the principal galaxy
        self.galaxy_name = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        """

        # Create a new Extractor instance
        if arguments.config is not None: extractor = cls(arguments.config)
        elif arguments.settings is not None: extractor = cls(arguments.settings)
        else: extractor = cls()

        # Ignore and special region
        if arguments.ignore is not None: extractor.config.ignore_region = arguments.ignore
        if arguments.special is not None: extractor.config.special_region = arguments.special

        # Set the input and output path
        if arguments.input_path is not None: extractor.config.input_path = arguments.input_path
        if arguments.output_path is not None: extractor.config.output_path = arguments.output_path

        # Set options for writing out regions or masks
        if arguments.regions: extractor.config.write_regions = True
        if arguments.masks: extractor.config.write_masked_frames = True

        # Writing output catalogs
        if arguments.catalogs:

            extractor.config.write_galactic_catalog = True
            extractor.config.writing.galactic_catalog_path = "galaxies.cat"

            extractor.config.write_stellar_catalog = True
            extractor.config.writing.stellar_catalog_path = "stars.cat"

        # Writing segmentation map
        if arguments.segments: extractor.config.write_segments = True

        # Synchronize catalog
        if arguments.synchronize: extractor.config.synchronize_catalogs = True

        # If used from the command line, the result and mask should be written
        extractor.config.write_result = True
        extractor.config.write_mask = True

        # Set the paths to the resulting frame and the total mask
        extractor.config.writing.result_path = "result.fits"
        extractor.config.writing.mask_path = "mask.fits"

        # Manual indices
        if arguments.not_stars is not None: extractor.config.stars.manual_indices.not_stars = arguments.not_stars
        if arguments.remove_stars is not None: extractor.config.stars.manual_indices.remove_stars = arguments.remove_stars
        if arguments.not_saturation is not None: extractor.config.stars.manual_indices.not_saturation = arguments.not_saturation
        if arguments.only_saturation is not None: extractor.config.stars.manual_indices.only_saturation = arguments.only_saturation

        # Options for using a file as input catalog
        if arguments.filecatalog:

            extractor.config.catalogs.galaxies.use_catalog_file = True
            extractor.config.catalogs.galaxies.catalog_path = "galaxies.cat"

            extractor.config.catalogs.stars.use_catalog_file = True
            extractor.config.catalogs.stars.catalog_path = "stars.cat"

        # The interpolation method
        if arguments.interpolation_method is not None: extractor.config.interpolation_method = arguments.interpolation_method

        # Return the new instance
        return extractor

    # -----------------------------------------------------------------

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup(image)

        # 2. Create the directory that will contain the output
        self.create_output_path()

        # 3. Import galactic and stellar catalogs
        self.import_catalogs()

        # 3. Run the galaxy extraction
        self.extract_galaxies()
        
        # 4. Run the star extraction
        self.extract_stars()

        # 5. Look for other sources
        self.find_other_sources()

        # 6. Build and update catalog
        self.build_and_synchronize_catalog()

        # 5. Writing phase
        self.write()

    # -----------------------------------------------------------------

    def setup(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # -- Create children --

        self.add_child("catalog_importer", CatalogImporter, self.config.catalogs)
        self.add_child("galaxy_extractor", GalaxyExtractor, self.config.galaxies)
        self.add_child("star_extractor", StarExtractor, self.config.stars)
        self.add_child("trained_extractor", TrainedExtractor, self.config.other_sources)
        self.add_child("catalog_builder", CatalogBuilder, self.config.building)
        self.add_child("catalog_synchronizer", CatalogSynchronizer, self.config.synchronization)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(Extractor, self).setup()

        # Inform the user
        log.info("Setting up the extractor ...")

        # Make a local reference to the image (mask inside image)
        self.image = image

        # Create a mask with shape equal to the shape of the frame
        self.mask = Mask.empty_like(self.image.frames.primary)
        self.image.add_mask(self.mask, "sources")

        # Set the appropriate configuration settings for writing out the galactic and stellar statistics
        if self.config.write_statistics:

            # Galaxy extractor
            self.galaxy_extractor.config.write_statistics = True
            self.galaxy_extractor.config.writing.statistics_path = "galaxies.stat"

            # Star extractor
            self.star_extractor.config.write_statistics = True
            self.star_extractor.config.writing.statistics_path = "saturation.stat"

        # Set the appropriate configuration settings for writing out the region files
        if self.config.write_regions:

            # Galaxy extractor
            self.galaxy_extractor.config.write_region = True
            self.galaxy_extractor.config.writing.region_path = "galaxies.reg"

            # Star extractor
            self.star_extractor.config.write_star_region = True
            self.star_extractor.config.write_saturation_region = True
            self.star_extractor.config.writing.star_region_path = "stars.reg"
            self.star_extractor.config.writing.saturation_region_path = "saturation.reg"

        # Set the appropriate configuration settings for writing out the masked frames
        if self.config.write_masked_frames:

            # Galaxy extractor
            self.galaxy_extractor.config.write_masked_frame = True
            self.galaxy_extractor.config.writing.masked_frame_path = "masked_galaxies.fits"

            # Star extractor
            self.star_extractor.config.write_masked_frame = True
            self.star_extractor.config.writing.masked_frame_path = "masked_stars.fits"

        # Write segmentation maps
        if self.config.write_segments:

            self.galaxy_extractor.config.write_segments = True
            self.galaxy_extractor.config.writing.segments_path = "galaxy_segments.fits"

            self.star_extractor.config.write_star_segments = True
            self.star_extractor.config.writing.star_segments_path = "star_segments.fits"

            self.star_extractor.config.write_saturation_segments = True
            self.star_extractor.config.writing.saturation_segments_path = "saturation_segments.fits"

            self.trained_extractor.config.write_segments = True
            self.trained_extractor.config.writing.segments_path = "other_segments.fits"

        # Set the interpolation method wherever it is appropriate
        self.galaxy_extractor.config.removal.interpolation_method = self.config.interpolation_method
        self.galaxy_extractor.config.manual.interpolation_method = self.config.interpolation_method
        self.star_extractor.config.removal.interpolation_method = self.config.interpolation_method
        self.star_extractor.config.saturation.interpolation_method = self.config.interpolation_method
        self.star_extractor.config.saturation.aperture_removal.interpolation_method = self.config.interpolation_method
        self.star_extractor.config.manual.interpolation_method = self.config.interpolation_method
        self.trained_extractor.config.removal.interpolation_method = self.config.interpolation_method

        # Create special and ignore mask
        self.set_ignore_mask()
        self.set_special_mask()

    # -----------------------------------------------------------------

    def create_output_path(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the output directory if necessary ...")

        # Create the directory if necessary
        if self.config.output_path is not None: filesystem.create_directory(self.config.output_path)

    # -----------------------------------------------------------------

    def import_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Importing the catalogs ...")

        # Run the catalog importer
        self.catalog_importer.run(self.image.frames.primary)

        # Inform the user
        log.success("Catalog imported")

    # -----------------------------------------------------------------

    def extract_galaxies(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Extracting the galaxies ...")

        # Run the galaxy extractor
        self.galaxy_extractor.run(self.image, self.catalog_importer.galactic_catalog, special=self.special_mask, ignore=self.ignore_mask)

        # Set the name of the principal galaxy
        self.galaxy_name = self.galaxy_extractor.principal.name

        # Inform the user
        log.success("Galaxies extracted")

    # -----------------------------------------------------------------
    
    def extract_stars(self):
        
        """
        This function ...
        """

        # Run the star extraction if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.image.wavelength is None or self.image.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Extracting the stars ...")

            # Run the star extractor
            self.star_extractor.run(self.image, self.galaxy_extractor, self.catalog_importer.stellar_catalog, special=self.special_mask, ignore=self.ignore_mask)

            # Inform the user
            log.success("Stars extracted")

        # No star subtraction for this image
        else: log.info("Star subtraction will not be performed for this image")

    # -----------------------------------------------------------------

    def find_other_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for sources in the frame not in the catalog ...")

        # If the wavelength of this image is greater than 25 micron, don't classify the sources that are found
        if self.image.wavelength is not None and self.image.wavelength > wavelengths.ranges.ir.mir.max: self.trained_extractor.config.classify = False
        else: self.trained_extractor.config.classify = True

        # Run the trained extractor just to find sources
        self.trained_extractor.run(self.image, self.galaxy_extractor, self.star_extractor, special=self.special_mask, ignore=self.ignore_mask)

        # Inform the user
        log.success("Other sources extracted")

    # -----------------------------------------------------------------

    def build_and_synchronize_catalog(self):

        """
        This function ...
        :return:
        """

        # Build the catalog
        self.build_catalog()

        # Synchronize the catalog
        if self.config.synchronize_catalogs: self.synchronize_catalog()

    # -----------------------------------------------------------------

    def build_catalog(self):

        """
        This function ...
        :return:
        """

        # Build the stellar catalog if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.image.wavelength is None or self.image.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Building the stellar catalog ...")

            # Run the catalog builder
            self.catalog_builder.run(self.image.frames.primary, self.galaxy_extractor, self.star_extractor, self.trained_extractor)

            # Inform the user
            log.success("Stellar catalog built")

    # -----------------------------------------------------------------

    def synchronize_catalog(self):

        """
        This function ...
        :return:
        """

        # Synchronize the catalog if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.image.wavelength is None or self.image.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Synchronizing with the DustPedia catalog ...")

            # Run the catalog synchronizer
            self.catalog_synchronizer.run(self.image.frames.primary, self.galaxy_name, self.catalog_builder.galactic_catalog, self.catalog_builder.stellar_catalog)

            # Inform the user
            log.success("Catalog synchronization done")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # If requested, write out the result
        if self.config.write_result: self.write_result()

        # If requested, write out the total mask
        if self.config.write_mask: self.write_mask()

        # If requested, write out the galactic region
        if self.config.write_galactic_region: self.write_galactic_region()

        # If requested, write out the stellar region
        if self.config.write_stellar_region: self.write_stellar_region()

        # If requested, write out the galactic catalog
        if self.config.write_galactic_catalog: self.write_galactic_catalog()

        # If requested, write out the compiled stellar catalog
        if self.config.write_stellar_catalog: self.write_stellar_catalog()

    # -----------------------------------------------------------------

    def write_galactic_region(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_stellar_region(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.galactic_catalog_path)
        if path is None:
            log.error("Galactic catalog path is not defined, skipping writing galactic catalog ...")
            return

        # Inform the user
        log.info("Writing galactic catalog to " + path + " ...")

        # Write the catalog to file
        tables.write(self.catalog_builder.galactic_catalog, path)

    # -----------------------------------------------------------------

    def write_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.stellar_catalog_path)
        if path is None:
            log.error("Stellar catalog path is not defined, skipping writing stellar catalog ...")
            return

        # Inform the user
        log.info("Writing stellar catalog to " + path + " ...")

        # Write the catalog to file
        tables.write(self.catalog_builder.stellar_catalog, path)

    # -----------------------------------------------------------------

    def write_result(self, header=None):

        """
        This function ...
        :param header:
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(self.config.writing.result_path)
        if path is None:
            log.error("Result path is not defined, skipping writing result ...")
            return

        # Inform the user
        log.info("Writing resulting frame to " + path + " ...")

        # Write out the resulting frame
        self.image.frames.primary.save(path, header, origin=self.name)

    # -----------------------------------------------------------------

    def write_mask(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the mask file
        path = self.full_output_path(self.config.writing.mask_path)
        if path is None:
            log.error("Mask path is not defined, skipping writing mask ...")
            return

        # Inform the user
        log.info("Writing the total mask to " + path + " ...")

        # Create a frame for the total mask
        frame = Frame(self.mask.astype(int))

        # Write out the total mask
        frame.save(path, origin=self.name)

    # -----------------------------------------------------------------

    def set_special_mask(self):

        """
        This function ...
        :param path:
        :return:
        """

        # If no special region is defined
        if self.config.special_region is None: return

        # Determine the full path to the special region file
        path = self.full_input_path(self.config.special_region)

        # Inform the user
        log.info("Creating mask covering objects that require special attention from " + path + " ...")

        # Load the region and create a mask from it
        region = Region.from_file(path)
        special_mask = Mask.from_region(region, self.image.xsize, self.image.ysize)

        # Create the mask
        self.special_mask = special_mask

    # -----------------------------------------------------------------

    def set_ignore_mask(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # If no ignore region is defined
        if self.config.ignore_region is None: return

        # Determine the full path to the ignore region file
        path = self.full_input_path(self.config.ignore_region)

        # Inform the user
        log.info("Creating mask covering objects that should be ignored from " + path + " ...")

        # Load the region and create a mask from it
        region = Region.from_file(path)
        ignore_mask = Mask.from_region(region, self.image.xsize, self.image.ysize)

        # Create the mask
        self.ignore_mask = ignore_mask

# -----------------------------------------------------------------
