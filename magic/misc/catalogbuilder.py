#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalogbuilder Contains the CatalogBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant AstroMagic classes and modules
from ..basics import Position, Extent, Rectangle

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import inspection

# -----------------------------------------------------------------

class CatalogBuilder(Configurable):

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
        super(CatalogBuilder, self).__init__(config, "magic")

        self.galaxy_extractor = None
        self.star_extractor = None

    # -----------------------------------------------------------------

    def run(self, frame, galaxy_extractor, star_extractor):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galaxy_extractor, star_extractor)
        
        # 2. 
        self.build()

    # -----------------------------------------------------------------

    def setup(self, frame, galaxy_extractor, star_extractor):

        """
        This function ...
        :param frame:
        :param galaxy_extractor:
        :param star_extractor:
        :return:
        """

        # Call the setup function of the base class
        super(CatalogBuilder, self).setup()

        self.frame = frame
        self.galaxy_extractor = galaxy_extractor
        self.star_extractor = star_extractor

    # -----------------------------------------------------------------

    def build(self):

        """
        This function ...
        :return:
        """

        # Get coordinate range
        center, ra_span, dec_span = self.frame.coordinate_range(silent=True)

        # Right ascension and declination in degrees
        ra = center.ra.degree
        dec = center.dec.degree

        # Right ascension and declination span in degrees
        ra_span = ra_span.value
        dec_span = dec_span.value

        # Create rectangle
        center = Position(ra, dec)
        radius = Extent(0.5 * ra_span, 0.5 * dec_span)
        box = Rectangle(center, radius)

        ## First check for identical entry



        # Name of the principal galaxy
        galaxy_name = self.galaxy_extractor.principal.name

        # Determine the path to the user catalogs directory
        catalogs_user_path = os.path.join(inspection.pts_user_dir, "magic", "catalogs")

        # Determine the path to the directory to contain the catalogs for this galaxy
        galaxy_user_path = os.path.join(catalogs_user_path, galaxy_name)

        # Cache the galaxy and stellar catalog
        if os.path.isdir(galaxy_user_path):

            old_galactic_catalog_path = os.path.join(galaxy_user_path, "galaxies.dat")
            old_stellar_catalog_path = os.path.join(galaxy_user_path, "stars.dat")

            if os.path.isfile(old_galactic_catalog_path):

                # Open the 'old' galaxy catalog
                old_galaxy_catalog = tables.from_file(old_galactic_catalog_path)

                # Create merged galaxy catalog
                galaxy_catalog = catalogs.merge_galactic_catalogs(galactic_catalog, old_galaxy_catalog)

                # Save the merged catalog
                path = os.path.join(galaxy_user_path, "galaxies.dat")
                tables.write(galaxy_catalog, path)

            # If a galactic catalog file does not exist yet
            else: filesystem.copy_file(galactic_catalog_path, galaxy_user_path, "galaxies.dat")

            # Check whether there is a stellar catalog file in the galaxy's directory
            if os.path.isfile(old_stellar_catalog_path):

                # Open the new stellar catalog
                stellar_catalog = tables.from_file(stellar_catalog_path)

                # Open the 'old' stellar catalog
                old_stellar_catalog = tables.from_file(old_stellar_catalog_path)

                # Create merged stellar catalog
                stellar_catalog = catalogs.merge_stellar_catalogs(stellar_catalog, old_stellar_catalog)

                # Save the merged catalog
                path = os.path.join(galaxy_user_path, "stars.dat")
                tables.write(stellar_catalog, path)

            # If a stellar catalog file does not exist yet
            else: filesystem.copy_file(stellar_catalog_path, galaxy_user_path, "stars.dat")

        else:

            # Create the directory to contain the catalogs for this galaxy
            filesystem.create_directory(galaxy_user_path)

            # Copy the galaxy and stellar catalog files into the new directory
            #filesystem.copy_file(galactic_catalog_path, galaxy_user_path, "galaxies.dat")
            #filesystem.copy_file(stellar_catalog_path, galaxy_user_path, "stars.dat")

            # Save galactic catalog
            galactic_catalog_path = os.path.join(galaxy_user_path, "galaxies.cat")
            tables.write(self.galaxy_extractor.catalog, galactic_catalog_path)

            # Save stellar catalog
            stellar_catalog_path = os.path.join(galaxy_user_path, "stars.cat")
            tables.write(self.star_extractor.catalog, stellar_catalog_path)

# -----------------------------------------------------------------
