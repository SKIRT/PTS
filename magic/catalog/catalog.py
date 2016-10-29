#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalog.catalog Contains the GalacticCatalog and StellarCatalog classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import catalogs
from ...core.tools import introspection, tables
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

catalogs_user_path = fs.join(introspection.pts_user_dir, "catalogs")

# -----------------------------------------------------------------

class GalacticCatalog(object):

    """
    This class ...
    """

    def __init__(self, frame_or_wcs):

        """
        The constructor ...
        :param frame_or_wcs:
        :return:
        """

        # Create the catalogs user directory if necessary
        if not fs.is_directory(catalogs_user_path): fs.create_directory(catalogs_user_path)

        # Determine the path to the 'galaxies' catalog path
        galaxies_catalog_path = fs.join(catalogs_user_path, "galaxies")

        # Create the catalogs/galaxies directory is necessary
        if not fs.is_directory(galaxies_catalog_path): fs.create_directory(galaxies_catalog_path)

        # Get the center coordinate and the range of RA and DEC
        center, ra_span, dec_span = frame_or_wcs.coordinate_range

        # Generate a unique string for the coordinate range
        name = str(center) + "_" + str(ra_span) + "_" + str(dec_span)

        # Determine the path to the catalog file
        self.path = fs.join(galaxies_catalog_path, name + ".cat")

        # Check whether the local file exists
        if not fs.is_file(self.path):

            # Get the table
            self.table = catalogs.create_galaxy_catalog(frame_or_wcs.bounding_box)

            # Save the table
            tables.write(self.table, self.path, format="ascii.ecsv")

        # Load the table
        else: self.table = tables.from_file(self.path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        tables.write(self.table, path, format="ascii.ecsv")

# -----------------------------------------------------------------

class StellarCatalog(object):

    """
    This class ...
    """

    def __init__(self, frame_or_wcs, catalog_names="II/246"):

        """
        This function ...
        :param frame_or_wcs:
        :param catalog_names:
        :return:
        """

        # Create the catalogs user directory if necessary
        if not fs.is_directory(catalogs_user_path): fs.create_directory(catalogs_user_path)

        # Determine the path to the 'galaxies' catalog path
        stars_catalog_path = fs.join(catalogs_user_path, "stars")

        # Create the catalogs/stars directory is necessary
        if not fs.is_directory(stars_catalog_path): fs.create_directory(stars_catalog_path)

        # Get the center coordinate and the range of RA and DEC
        center, ra_span, dec_span = frame_or_wcs.coordinate_range

        # Generate a unique string for the coordinate range
        name = str(center) + "_" + str(ra_span) + "_" + str(dec_span)

        # Determine the path to the catalog file
        self.path = fs.join(stars_catalog_path, name + ".cat")

        # Check whether the local file exists
        if not fs.is_file(self.path):

            # Get the table
            self.table = catalogs.create_star_catalog(frame_or_wcs.bounding_box, frame_or_wcs.pixelscale, catalog_names)

            # Save the table
            tables.write(self.table, self.path, format="ascii.ecsv")

        # Load the table
        else: self.table = tables.from_file(self.path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        tables.write(self.table, path, format="ascii.ecsv")

# -----------------------------------------------------------------
