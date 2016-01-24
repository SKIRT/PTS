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

# Import standard modules
import os

# Import astronomical modules
import astropy.units as u

# Import the relevant AstroMagic classes and modules
from ..core import Frame
from ..basics import Mask, Region
from .galaxyextraction import GalaxyExtractor
from .starextraction import StarExtractor
from .trainedextractor import TrainedExtractor
from ..misc import CatalogBuilder
from ..basics import CatalogCoverage

# Import the relevant PTS classes and modules
from ...core.tools import filesystem, tables, inspection
from ...core.basics.configurable import Configurable

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

        # The image frame on which to perform the extraction
        self.frame = None

        # The mask covering pixels that should be ignored throughout the entire extraction procedure
        self.input_mask = None
        self.special_mask = None
        self.ignore_mask = None

        # The output mask
        self.mask = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        """

        # Create a new Extractor instance
        extractor = cls(arguments.config)

        # Report file
        extractor.config.logging.path = "log.txt"

        # Ignore and special region
        extractor.config.ignore_region = arguments.ignore
        extractor.config.special_region = arguments.special

        # Set the input and output path
        extractor.config.input_path = arguments.input_path
        extractor.config.output_path = arguments.output_path

        # Set options for writing out regions or masks
        extractor.config.write_regions = arguments.regions
        extractor.config.write_masked_frames = arguments.masks
        extractor.config.write_catalogs = arguments.catalogs

        # Build catalog
        extractor.config.build_catalog = arguments.build

        # If used from the command line, the result and mask should be written
        extractor.config.write_result = True
        extractor.config.write_mask = True

        # Manual indices
        if arguments.not_stars is not None: extractor.config.stars.manual_indices.not_stars = arguments.not_stars
        if arguments.remove_stars is not None: extractor.config.stars.manual_indices.remove_stars = arguments.remove_stars
        if arguments.not_saturation is not None: extractor.config.stars.manual_indices.not_saturation = arguments.not_saturation

        if arguments.filecatalog:

            extractor.config.galaxies.fetching.use_catalog_file = True
            extractor.config.galaxies.fetching.catalog_path = "galaxies.cat"

            extractor.config.stars.fetching.use_catalog_file = True
            extractor.config.stars.fetching.catalog_path = "stars.cat"

        # Return the new instance
        return extractor

    # -----------------------------------------------------------------

    def run(self, frame, input_mask):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, input_mask)

        # 2. Create the directory that will contain the output
        self.create_output_path()

        # 3. Import cached catalogs
        self.import_catalogs()

        # 3. Run the galaxy extraction
        self.extract_galaxies()
        
        # 4. Run the star extraction
        self.extract_stars()

        # 5. Look for other sources
        self.find_other_sources()

        # 6. Build catalogs
        if self.config.build_catalog: self.build_catalogs()

        # 5. Writing phase
        self.write()

    # -----------------------------------------------------------------

    def setup(self, frame, input_mask):

        """
        This function ...
        :return:
        """

        # -- Create children --

        self.add_child("galaxy_extractor", GalaxyExtractor, self.config.galaxies)
        self.add_child("star_extractor", StarExtractor, self.config.stars)
        self.add_child("trained_extractor", TrainedExtractor, self.config.other_sources)
        self.add_child("catalog_builder", CatalogBuilder)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(Extractor, self).setup()

        # Inform the user
        self.log.info("Setting up the extractor ...")

        # Make a local reference to the frame and mask
        self.frame = frame
        self.input_mask = input_mask

        # Set the paths to the resulting frame and the total mask
        self.config.writing.result_path = "result.fits"
        self.config.writing.mask_path = "mask.fits"

        # Create a mask with shape equal to the shape of the frame
        self.mask = Mask.from_shape(self.frame.shape)

        # Set the appropriate configuration settings for writing out the galactic and stellar catalogs
        if self.config.write_catalogs:

            # Galaxy extractor
            self.galaxy_extractor.config.write_catalog = True
            self.galaxy_extractor.config.writing.catalog_path = "galaxies.cat"

            # Star extractor
            self.star_extractor.config.write_catalog = True
            self.star_extractor.config.writing.catalog_path = "stars.cat"

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
        self.log.info("Creating the output directory if necessary ...")

        # Create the directory if necessary
        if self.config.output_path is not None: filesystem.create_directory(self.config.output_path)

    # -----------------------------------------------------------------

    def import_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Checking whether cached catalogs can be imported ...")

        # Determine the path to the user catalogs directory
        catalogs_user_path = os.path.join(inspection.pts_user_dir, "magic", "catalogs")

        # Get bounding box of the frame
        bounding_box = self.frame.bounding_box()

        # Loop over all directories within the catalogs directory (different galaxies)
        for galaxy_path in filesystem.directories_in_path(catalogs_user_path):

            # Get the galaxy name
            galaxy_name = os.path.basename(galaxy_path)

            # Get the catalog coverage for this galaxy
            coverage = CatalogCoverage(galaxy_name)

            # If this is the galaxy matching the frame, check if the current catalog covers the entire frame
            # Second part is probably time-consuming .. But if first is False, won't be executed
            if coverage.matches(bounding_box) and coverage.covers(bounding_box):

                # Debug info
                self.log.debug("DustPedia catalog will be imported for " + galaxy_name + " galaxy")

                galactic_catalog_path = os.path.join(galaxy_path, "galaxies.cat")
                stellar_catalog_path = os.path.join(galaxy_path, "stars.cat")

                self.galaxy_catalog = tables.from_file(galactic_catalog_path)
                self.star_catalog = tables.from_file(stellar_catalog_path)

                break

        # If no break is encountered
        else:

            self.galaxy_catalog = None
            self.star_catalog = None

    # -----------------------------------------------------------------

    def extract_galaxies(self):
        
        """
        This function ...
        """

        # Inform the user
        self.log.info("Extracting the galaxies ...")

        # Run the galaxy extractor
        self.galaxy_extractor.run(self.frame, self.input_mask, self.galaxy_catalog, special=self.special_mask, ignore=self.ignore_mask)

        # Add to the total mask
        self.mask += self.galaxy_extractor.mask

    # -----------------------------------------------------------------
    
    def extract_stars(self):
        
        """
        This function ...
        """

        # Inform the user
        self.log.info("Extracting the stars ...")

        # ...
        if self.frame.wavelength is None or self.frame.wavelength < 10.0 * u.Unit("micron"):

            # Run the star extractor
            self.star_extractor.run(self.frame, self.input_mask, self.galaxy_extractor, self.star_catalog, special=self.special_mask, ignore=self.ignore_mask)

            # Add star mask to the total mask
            self.mask += self.star_extractor.mask

        else: self.log.info("No star extraction for this image")

    # -----------------------------------------------------------------

    def find_other_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Looking for sources in the frame not in the catalog ...")

        # Run the trained extractor just to find sources
        # As the mask, pass the input mask (bad pixels) + self.mask (at this point, combined galaxy and star extraction masks)
        self.trained_extractor.run(self.frame, self.input_mask + self.mask, self.galaxy_extractor, self.star_extractor, special=self.special_mask, ignore=self.ignore_mask)

        # Add sources to the total mask
        self.mask += self.trained_extractor.mask

    # -----------------------------------------------------------------

    def build_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Building the catalog for later ...")

        # Run the catalog builder
        self.catalog_builder.run(self.frame, self.galaxy_extractor, self.star_extractor, self.trained_extractor)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing ...")

        # If requested, write out the result
        if self.config.write_result: self.write_result()

        # If requested, write out the total mask
        if self.config.write_mask: self.write_mask()

    # -----------------------------------------------------------------

    def write_result(self, header=None):

        """
        This function ...
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(self.config.writing.result_path)

        # Inform the user
        self.log.info("Writing resulting frame to " + path + " ...")

        # Write out the resulting frame
        self.frame.save(path, header)

    # -----------------------------------------------------------------

    def write_mask(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the mask file
        path = self.full_output_path(self.config.writing.mask_path)

        # Inform the user
        self.log.info("Writing the total mask to " + path + " ...")

        # Create a frame for the total mask
        frame = Frame(self.mask.astype(int))

        # Write out the total mask
        frame.save(path)

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
        self.log.info("Setting special region from " + path)

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)
        special_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Return the mask
        #return special_mask
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
        self.log.info("Setting region to ignore for subtraction from " + path)

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)
        ignore_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Return the mask
        #return ignore_mask
        self.ignore_mask = ignore_mask

# -----------------------------------------------------------------
