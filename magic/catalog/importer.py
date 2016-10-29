#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalog.importer Contains the CatalogImporter class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from config import Sequence

# Import the relevant PTS classes and modules
from ..basics.catalogcoverage import CatalogCoverage
from ..tools import catalogs
from ...core.basics.configurable import Configurable
from ...core.tools import introspection, tables
from ...core.tools import filesystem as fs
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
        super(CatalogImporter, self).__init__(config)

        # The catalog box
        self.coordinate_box = None

        # The name of the galaxy
        self._has_dustpedia_catalog = None
        self.galaxy_name = None

        # The galactic and stellar catalogs
        self.galactic_catalog = None
        self.stellar_catalog = None

        # Determine the path to the catalogs user directory (where the DustPedia catalogs are now stored)
        self.catalogs_user_path = fs.join(introspection.pts_user_dir, "magic", "catalogs")

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)
        
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
        self.coordinate_box = None
        self.galaxy_name = None
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(CatalogImporter, self).setup(**kwargs)

        # Set the coordinate box
        self.coordinate_box = kwargs.pop("coordinate_box")

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
        bounding_box = self.coordinate_box

        # Loop over all directories within the catalogs directory (different galaxies)
        for galaxy_path in fs.directories_in_path(self.catalogs_user_path):

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
        self.galactic_catalog = tables.from_file(path, format="ascii.commented_header")

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
        self.stellar_catalog = tables.from_file(path, format="ascii.commented_header")

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
        self.galactic_catalog = catalogs.create_galaxy_catalog(self.coordinate_box)

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
        self.stellar_catalog = catalogs.create_star_catalog(self.coordinate_box, catalog_list)

        # Inform the user
        log.debug("Number of stars: " + str(len(self.stellar_catalog)))

    # -----------------------------------------------------------------

    def write_galactic_catalog_to(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Writing galactic catalog to " + path + " ...")

        # Write the catalog to file
        tables.write(self.galactic_catalog, path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def write_stellar_catalog_to(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Writing stellar catalog to " + path + " ...")

        # Write the catalog to file
        tables.write(self.stellar_catalog, path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def write_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.galactic_catalog_path)

        # Write
        self.write_galactic_catalog_to(path)

    # -----------------------------------------------------------------

    def write_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.stellar_catalog_path)

        # Write
        self.write_stellar_catalog_to(path)

# -----------------------------------------------------------------
