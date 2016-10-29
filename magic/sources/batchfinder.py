#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.batchfinder Contains the BatchSourceFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules

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
from ..core.dataset import DataSet
from ..catalog.importer import CatalogImporter
from ...core.tools import filesystem as fs
from ..basics.skyregion import SkyRegion

# -----------------------------------------------------------------

class BatchSourceFinder(Configurable):

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
        super(BatchSourceFinder, self).__init__(config)

        # -- Attributes --

        # The data set
        self.dataset = None

        # The galactic and stellar catalog
        self.galactic_catalog = None
        self.stellar_catalog = None

        # The mask covering pixels that should be ignored throughout the entire extraction procedure
        self.special_mask = None
        self.ignore_mask = None
        self.bad_mask = None

        # The name of the principal galaxy
        self.galaxy_name = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get catalogs
        self.get_catalogs()

        # 3. Find the galaxies
        self.find_galaxies()
        
        # 3. Find the stars
        if self.config.find_stars: self.find_stars()

        # 4. Look for other sources
        if self.config.find_other_sources: self.find_other_sources()

        # 5. Build and update catalog
        #self.build_and_synchronize_catalog()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BatchSourceFinder, self).setup(**kwargs)

        # Load the dataset
        self.dataset = DataSet.from_file(self.config.dataset)

        # Load special region
        self.special_region = SkyRegion.from_file(self.config.special_region) if self.config.special_region is not None else None

        # Load ignore region
        self.ignore_region = SkyRegion.from_file(self.config.ignore_region) if self.config.ignore_region is not None else None

        # Create the finders
        self.galaxy_finder = GalaxyFinder(self.config.galaxies)
        self.star_finder = StarFinder(self.config.stars)
        self.trained_finder = TrainedFinder(self.config.other_sources)

        #self.galactic_catalog = kwargs.pop("galactic_catalog")
        #self.stellar_catalog = kwargs.pop("stellar_catalog")




        #self.add_child("galaxy_finder", GalaxyFinder, self.config.galaxies)
        #self.add_child("star_finder", StarFinder, self.config.stars)
        #self.add_child("trained_finder", TrainedFinder, self.config.other_sources)
        #self.add_child("catalog_builder", CatalogBuilder, self.config.building)
        #self.add_child("catalog_synchronizer", CatalogSynchronizer, self.config.synchronization)

        # -- Setup of the base class --

        # Call the setup function of the base class
        #super(BatchSourceFinder, self).setup()

        # Inform the user
        #log.info("Setting up the batch source finder ...")

        # Set the galactic and stellar catalog
        #self.galactic_catalog = galactic_catalog
        #self.stellar_catalog = stellar_catalog

        # Set the special and ignore mask
        #self.special_mask = Mask.from_region(special_region, self.frame.xsize, self.frame.ysize) if special_region is not None else None
        #self.ignore_mask = Mask.from_region(ignore_region, self.frame.xsize, self.frame.ysize) if ignore_region is not None else None


        # Set a reference to the mask of bad pixels
        #self.bad_mask = bad_mask

    # -----------------------------------------------------------------

    def get_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting galactic and stellar catalogs ...")

        # Create a CatalogImporter instance
        catalog_importer = CatalogImporter()

        # Set configuration
        if self.config.galactic_catalog_file is not None:
            catalog_importer.config.galaxies.use_catalog_file = True
            catalog_importer.config.galaxies.catalog_path = self.config.galactic_catalog_file

        if self.config.stellar_catalog_file is not None:
            catalog_importer.config.stars.use_catalog_file = True
            catalog_importer.config.stars.catalog_path = self.config.stellar_catalog_file

        # Get the coordinate box and minimum pixelscale
        coordinate_box = self.dataset.get_bounding_box()
        min_pixelscale = self.dataset.min_pixelscale

        # Run the catalog importer
        catalog_importer.run(coordinate_box=coordinate_box, pixelscale=min_pixelscale)  # work with coordinate box instead ? image.coordinate_box ?

        # Set the catalogs
        self.galactic_catalog = catalog_importer.galactic_catalog
        self.stellar_catalog = catalog_importer.stellar_catalog

    # -----------------------------------------------------------------

    def find_galaxies(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the galaxies ...")

        # Loop over the images
        for name in self.dataset.names:

            # Get the frame
            frame = self.dataset.get_frame(name)

            # Run the galaxy finder
            self.galaxy_finder.run(frame, self.galactic_catalog, special=self.special_mask, ignore=self.ignore_mask, bad=self.bad_mask)

            # Set the name of the principal galaxy
            #self.galaxy_name = self.galaxy_finder.principal.name

            # Inform the user
            log.success("Finished finding the galaxies for '" + name + "' ...")

    # -----------------------------------------------------------------
    
    def find_stars(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the stars ...")

        # Loop over the images
        for name in self.dataset.names:

            # Get the frame
            frame = self.dataset.get_frame(name)

            # Run the star finder if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
            if frame.wavelength is None or frame.wavelength < wavelengths.ranges.ir.mir.max:

                # Inform the user
                log.info("Finding the stars ...")

                # Get the frame
                frame = self.dataset.get_frame(name)

                # Run the star finder
                self.star_finder.run(frame, self.galaxy_finder, self.stellar_catalog, special=self.special_mask, ignore=self.ignore_mask, bad=self.bad_mask)

                # Inform the user
                log.success("Finished finding the stars for '" + name + "' ...")

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

        # Loop over the frames
        for name in self.dataset.names:

            # Get the frame
            frame = self.dataset.get_frame(name)

            # If the wavelength of this image is greater than 25 micron, don't classify the sources that are found
            if frame.wavelength is not None and frame.wavelength > wavelengths.ranges.ir.mir.max: self.trained_finder.config.classify = False
            else: self.trained_finder.config.classify = True

            # Run the trained finder just to find sources
            self.trained_finder.run(frame, self.galaxy_finder, self.star_finder, special=self.special_mask, ignore=self.ignore_mask, bad=self.bad_mask)

            # Inform the user
            log.success("Finished finding other sources for '" + name + "' ...")

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

    def write(self):

        """
        This function ...
        :return:
        """

        # 1. Write
        self.write_galactic_catalogs()

        # 2. Write
        self.write_stellar_catalogs()

        # 3. Write ...
        self.write_statistics()

    # -----------------------------------------------------------------

    def write_galactic_catalogs(self):

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

    def write_stellar_catalogs(self):

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
