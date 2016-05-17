#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.finder Contains the SourceFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .galaxyfinder import GalaxyFinder
from .starfinder import StarFinder
from .trainedfinder import TrainedFinder
from ..basics.mask import Mask
from ..catalog.builder import CatalogBuilder
from ..catalog.synchronizer import CatalogSynchronizer
from ..tools import wavelengths
from ...core.tools import tables
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log

# -----------------------------------------------------------------

class SourceFinder(Configurable):

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
        super(SourceFinder, self).__init__(config, "magic")

        # -- Attributes --

        # The image frame
        self.frame = None

        # The original WCS
        self.original_wcs = None

        # The galactic and stellar catalog
        self.galactic_catalog = None
        self.stellar_catalog = None

        # The mask covering pixels that should be ignored throughout the entire extraction procedure
        self.special_mask = None
        self.ignore_mask = None
        self.bad_mask = None

        # The animation
        self.animation = None

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

        # Create a new SourceFinder instance
        if arguments.config is not None: finder = cls(arguments.config)
        elif arguments.settings is not None: finder = cls(arguments.settings)
        else: finder = cls()

        # Options for using a file as input catalog
        if arguments.filecatalog:

            finder.config.catalogs.galaxies.use_catalog_file = True
            finder.config.catalogs.galaxies.catalog_path = "galaxies.cat"

            finder.config.catalogs.stars.use_catalog_file = True
            finder.config.catalogs.stars.catalog_path = "stars.cat"

        # Set the downsample factor
        if arguments.downsample is not None: finder.config.downsample_factor = arguments.downsample

        # Return the new instance
        return finder

    # -----------------------------------------------------------------

    @property
    def galaxy_region(self):

        """
        This function ...
        :return:
        """

        if self.downsampled:
            sky_region = self.galaxy_sky_region
            return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        return self.galaxy_finder.region

    # -----------------------------------------------------------------

    @property
    def galaxy_sky_region(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_finder.region.to_sky(self.frame.wcs)

    # -----------------------------------------------------------------

    @property
    def star_region(self):

        """
        This function ...
        :return:
        """

        if self.downsampled:
            sky_region = self.star_sky_region
            return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        else: return self.star_finder.star_region

    # -----------------------------------------------------------------

    @property
    def star_sky_region(self):

        """
        This function ...
        :return:
        """

        return self.star_finder.star_region.to_sky(self.frame.wcs)

    # -----------------------------------------------------------------

    @property
    def saturation_region(self):

        """
        This function ...
        :return:
        """

        if self.downsampled:
            sky_region = self.saturation_sky_region
            return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        else: return self.star_finder.saturation_region

    # -----------------------------------------------------------------

    @property
    def saturation_sky_region(self):

        """
        This function ...
        :return:
        """

        return self.star_finder.saturation_region.to_sky(self.frame.wcs)

    # -----------------------------------------------------------------

    @property
    def other_region(self):

        """
        This function ...
        :return:
        """

        if self.downsampled:
            sky_region = self.other_sky_region
            return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        else: return self.trained_finder.region

    # -----------------------------------------------------------------

    @property
    def other_sky_region(self):

        """
        This function ...
        :return:
        """

        return self.trained_finder.region.to_sky(self.frame.wcs) if self.trained_finder.region is not None else None

    # -----------------------------------------------------------------

    @property
    def galaxy_segments(self):

        """
        This property ...
        :return:
        """

        if self.downsampled: return self.galaxy_finder.segments.upsampled(self.config.downsample_factor, integers=True) if self.galaxy_finder.segments is not None else None
        else: return self.galaxy_finder.segments

    # -----------------------------------------------------------------

    @property
    def star_segments(self):

        """
        This property ...
        :return:
        """

        if self.downsampled: return self.star_finder.segments.upsampled(self.config.downsample_factor, integers=True) if self.star_finder.segments is not None else None
        else: return self.star_finder.segments

    # -----------------------------------------------------------------

    @property
    def other_segments(self):

        """
        This property ...
        :return:
        """

        if self.downsampled: return self.trained_finder.segments.upsampled(self.config.downsample_factor, integers=True) if self.trained_finder.segments is not None else None
        return self.trained_finder.segments

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        return self.star_finder.fwhm

    # -----------------------------------------------------------------

    def run(self, frame, galactic_catalog, stellar_catalog, special_region=None, ignore_region=None, bad_mask=None, animation=None):

        """
        This function ...
        :param frame:
        :param galactic_catalog:
        :param stellar_catalog:
        :param special_region:
        :param ignore_region:
        :param bad_mask:
        :param animation:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galactic_catalog, stellar_catalog, special_region, ignore_region, bad_mask, animation)

        # 2. Find the galaxies
        self.find_galaxies()
        
        # 3. Find the stars
        if self.config.find_stars: self.find_stars()

        # 4. Look for other sources
        if self.config.find_other_sources: self.find_other_sources()

        # 5. Build and update catalog
        self.build_and_synchronize_catalog()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Base class implementation removes the children
        super(SourceFinder, self).clear()

        # Set default values for all attributes
        self.frame = None
        self.original_wcs = None
        self.galactic_catalog = None
        self.stellar_catalog = None
        self.special_mask = None
        self.ignore_mask = None
        self.bad_mask = None
        self.animation = None
        self.mask = None
        self.galaxy_name = None

    # -----------------------------------------------------------------

    def setup(self, frame, galactic_catalog, stellar_catalog, special_region, ignore_region, bad_mask=None, animation=None):

        """
        This function ...
        :param frame:
        :param galactic_catalog:
        :param stellar_catalog:
        :param special_region:
        :param ignore_region:
        :param bad_mask:
        :param animation:
        :return:
        """

        # -- Create children --

        self.add_child("galaxy_finder", GalaxyFinder, self.config.galaxies)
        self.add_child("star_finder", StarFinder, self.config.stars)
        self.add_child("trained_finder", TrainedFinder, self.config.other_sources)
        self.add_child("catalog_builder", CatalogBuilder, self.config.building)
        self.add_child("catalog_synchronizer", CatalogSynchronizer, self.config.synchronization)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(SourceFinder, self).setup()

        # Inform the user
        log.info("Setting up the source finder ...")

        # Downsample or just make a local reference to the image frame
        if self.downsampled:
            self.frame = frame.downsampled(self.config.downsample_factor)
            self.original_wcs = frame.wcs
        else: self.frame = frame

        # Set the galactic and stellar catalog
        self.galactic_catalog = galactic_catalog
        self.stellar_catalog = stellar_catalog

        # Set the special and ignore mask
        self.special_mask = Mask.from_region(special_region, self.frame.xsize, self.frame.ysize) if special_region is not None else None
        self.ignore_mask = Mask.from_region(ignore_region, self.frame.xsize, self.frame.ysize) if ignore_region is not None else None

        # Set a reference to the mask of bad pixels
        self.bad_mask = bad_mask

        # Make a reference to the animation
        self.animation = animation

    # -----------------------------------------------------------------

    @property
    def downsampled(self):

        """
        This function ...
        :return:
        """

        return self.config.downsample_factor is not None and self.config.downsample != 1

    # -----------------------------------------------------------------

    def find_galaxies(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the galaxies ...")

        # Run the galaxy finder
        self.galaxy_finder.run(self.frame, self.galactic_catalog, special=self.special_mask, ignore=self.ignore_mask, bad=self.bad_mask)

        # Set the name of the principal galaxy
        self.galaxy_name = self.galaxy_finder.principal.name

        # Inform the user
        log.success("Finished finding the galaxies")

    # -----------------------------------------------------------------
    
    def find_stars(self):
        
        """
        This function ...
        """

        # Run the star finder if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.frame.wavelength is None or self.frame.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Finding the stars ...")

            # Run the star finder
            self.star_finder.run(self.frame, self.galaxy_finder, self.stellar_catalog, special=self.special_mask, ignore=self.ignore_mask, bad=self.bad_mask)

            # Inform the user
            log.success("Finished finding the stars")

        # No star subtraction for this image
        else: log.info("Finding stars will not be performed on this frame")

    # -----------------------------------------------------------------

    def find_other_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources in the frame not in the catalog ...")

        # If the wavelength of this image is greater than 25 micron, don't classify the sources that are found
        if self.frame.wavelength is not None and self.frame.wavelength > wavelengths.ranges.ir.mir.max: self.trained_finder.config.classify = False
        else: self.trained_finder.config.classify = True

        # Run the trained finder just to find sources
        self.trained_finder.run(self.frame, self.galaxy_finder, self.star_finder, special=self.special_mask, ignore=self.ignore_mask, bad=self.bad_mask)

        # Inform the user
        log.success("Finished finding other sources")

    # -----------------------------------------------------------------

    def build_and_synchronize_catalog(self):

        """
        This function ...
        :return:
        """

        # Build the catalog
        if self.config.build_catalogs: self.build_catalog()

        # Synchronize the catalog
        if self.config.build_catalogs and self.config.synchronize_catalogs: self.synchronize_catalog()

    # -----------------------------------------------------------------

    def build_catalog(self):

        """
        This function ...
        :return:
        """

        # Build the stellar catalog if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.frame.wavelength is None or self.frame.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Building the stellar catalog ...")

            # Run the catalog builder
            self.catalog_builder.run(self.frame, self.galaxy_finder, self.star_finder, self.trained_finder)

            # Inform the user
            log.success("Stellar catalog built")

    # -----------------------------------------------------------------

    def synchronize_catalog(self):

        """
        This function ...
        :return:
        """

        # Synchronize the catalog if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.frame.wavelength is None or self.frame.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Synchronizing with the DustPedia catalog ...")

            # Run the catalog synchronizer
            self.catalog_synchronizer.run(self.frame.frames.primary, self.galaxy_name, self.catalog_builder.galactic_catalog, self.catalog_builder.stellar_catalog)

            # Inform the user
            log.success("Catalog synchronization done")

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

    def write_statistics(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Writing statistics to '" + path + "' ...")

        # Open the file, write the info
        with open(path, 'w') as statistics_file:
            statistics_file.write("FWHM: " + str(self.fwhm) + "\n")

# -----------------------------------------------------------------
