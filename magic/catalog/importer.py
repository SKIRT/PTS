#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalogimporter Contains the CatalogImporter class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
from config import Sequence

# Import the relevant AstroMagic classes and modules
from ..basics.catalogcoverage import CatalogCoverage
from ..tools import catalogs

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import inspection, tables, filesystem
from ...core.tools.logging import log

# -----------------------------------------------------------------

class CatalogImporter(Configurable):

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
        super(CatalogImporter, self).__init__(config, "magic")

        # The image frame
        self.frame = None

        # The name of the galaxy
        self._has_dustpedia_catalog = None
        self.galaxy_name = None

        # The galactic and stellar catalogs
        self.galactic_catalog = None
        self.stellar_catalog = None

        # Determine the path to the catalogs user directory (where the DustPedia catalogs are now stored)
        self.catalogs_user_path = os.path.join(inspection.pts_user_dir, "magic", "catalogs")

    # -----------------------------------------------------------------

    def run(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame)
        
        # 2. Import the galactic catalog
        self.import_galactic_catalog()

        # 3. Import the stellar catalog
        self.import_stellar_catalog()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set attributes to None
        self.frame = None
        self.galaxy_name = None
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def setup(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Call the setup function of the base class
        super(CatalogImporter, self).setup()

        # Set the frame
        self.frame = frame

    # -----------------------------------------------------------------

    def import_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Importing galactic catalog ...")

        # Use catalog file if specified
        if self.config.galaxies.use_catalog_file: self.import_galactic_catalog_from_file()
        elif self.has_dustpedia_catalog: self.import_galactic_dustpedia_catalog()
        else: self.fetch_galactic_catalog()

    # -----------------------------------------------------------------

    def import_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Importing stellar catalog ...")

        # Use catalog file if specified
        if self.config.stars.use_catalog_file: self.import_stellar_catalog_from_file()
        elif self.has_dustpedia_catalog: self.import_stellar_dustpedia_catalog()
        else: self.fetch_stellar_catalog()

    # -----------------------------------------------------------------

    @property
    def has_dustpedia_catalog(self):

        """
        This function ...
        :return:
        """

        # Avoid overhead when calling this function twice
        if self._has_dustpedia_catalog is not None: return self._has_dustpedia_catalog

        # Get bounding box of the frame
        bounding_box = self.frame.bounding_box()

        # Loop over all directories within the catalogs directory (different galaxies)
        for galaxy_path in filesystem.directories_in_path(self.catalogs_user_path):

            # Get the galaxy name
            galaxy_name = os.path.basename(galaxy_path)

            # Get the catalog coverage for this galaxy
            coverage = CatalogCoverage(galaxy_name)

            # Check whether this galaxy matches the current frame
            if coverage.matches(bounding_box):

                # Debug info
                log.debug("Found a match with DustPedia catalog for galaxy " + galaxy_name)

                self._has_dustpedia_catalog = True
                self.galaxy_name = galaxy_name

                break

        # If a break is not encountered
        else: self._has_dustpedia_catalog = False

        # Return the answer
        return self._has_dustpedia_catalog

    # -----------------------------------------------------------------

    def import_galactic_catalog_from_file(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_input_path(self.config.galaxies.catalog_path)

        # Inform the user
        log.info("Importing galactic catalog from file " + path + " ...")

        # Load the catalog
        self.galactic_catalog = tables.from_file(path)

    # -----------------------------------------------------------------

    def import_stellar_catalog_from_file(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_input_path(self.config.stars.catalog_path)

        # Inform the user
        log.info("Importing stellar catalog from file " + path + " ...")

        # Load the catalog
        self.stellar_catalog = tables.from_file(path)

    # -----------------------------------------------------------------

    def import_galactic_dustpedia_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Importing galactic DustPedia catalog for " + self.galaxy_name)

        # Determine the path to the DustPedia galactic catalog for this galaxy
        galaxy_path = os.path.join(self.catalogs_user_path, self.galaxy_name)
        galactic_catalog_path = os.path.join(galaxy_path, "galaxies.cat")

        # Load the galactic catalog
        self.galactic_catalog = tables.from_file(galactic_catalog_path)

    # -----------------------------------------------------------------

    def import_stellar_dustpedia_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Importing stellar DustPedia catalog for " + self.galaxy_name)

        # Determine the path to the DustPedia stellar catalog for this galaxy
        galaxy_path = os.path.join(self.catalogs_user_path, self.galaxy_name)
        stellar_catalog_path = os.path.join(galaxy_path, "stars.cat")

        # Load the stellar catalog
        self.stellar_catalog = tables.from_file(stellar_catalog_path)

    # -----------------------------------------------------------------

    def fetch_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy positions from an online catalog ...")

        # Create the galaxy catalog
        self.galactic_catalog = catalogs.create_galaxy_catalog(self.frame)

        # Inform the user
        log.debug("Number of galaxies: " + str(len(self.galactic_catalog)))

    # -----------------------------------------------------------------

    def fetch_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching star positions from online catalogs ...")

        # Check whether the 'catalogs' setting defines a single catalog name or a list of such names
        if isinstance(self.config.stars.fetching.catalogs, basestring): catalog_list = [self.config.stars.fetching.catalogs]
        elif isinstance(self.config.stars.fetching.catalogs, Sequence): catalog_list = self.config.stars.fetching.catalogs
        else: raise ValueError("Invalid option for 'catalogs', should be a string or a list of strings")

        # Create the star catalog
        self.stellar_catalog = catalogs.create_star_catalog(self.frame, catalog_list)

        # Inform the user
        log.debug("Number of stars: " + str(len(self.stellar_catalog)))

# -----------------------------------------------------------------
