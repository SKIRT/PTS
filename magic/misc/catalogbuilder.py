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
import numpy as np

# Import the relevant AstroMagic classes and modules
from ..basics import Position, Extent, Rectangle, CatalogCoverage
from ..tools import catalogs

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import inspection, tables, filesystem

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

        self.frame = None
        self.galaxy_extractor = None
        self.star_extractor = None
        self.trained_extractor = None

        self.galaxy_user_path = None

        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def run(self, frame, galaxy_extractor, star_extractor, trained_extractor):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galaxy_extractor, star_extractor, trained_extractor)
        
        # 2. Build the catalog
        self.build()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        self.frame = None
        self.galaxy_extractor = None
        self.star_extractor = None
        self.trained_extractor = None

        self.galaxy_name = None
        self.galaxy_user_path = None

    # -----------------------------------------------------------------

    def setup(self, frame, galaxy_extractor, star_extractor, trained_extractor):

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
        self.trained_extractor = trained_extractor

        # Name of the principal galaxy
        self.galaxy_name = self.galaxy_extractor.principal.name

        # Determine the path to the user catalogs directory
        catalogs_user_path = os.path.join(inspection.pts_user_dir, "magic", "catalogs")

        # Determine the path to the directory to contain the catalogs for this galaxy
        self.galaxy_user_path = os.path.join(catalogs_user_path, self.galaxy_name)

    # -----------------------------------------------------------------

    def build(self):

        """
        This function ...
        :return:
        """

        # Get bounding box
        coordinate_box = self.frame.bounding_box()

        # Get current catalog coverage
        coverage = CatalogCoverage(self.galaxy_name)

        # Cache the galaxy and stellar catalog
        if os.path.isdir(self.galaxy_user_path):

            # If the coordinate range of the current frame extents that of the previous frames
            if not coverage.covers(coordinate_box):

                # -- MERGE GALAXY CATALOG --

                # Check whether there is a galactic catalog file in the galaxy's directory
                old_galactic_catalog_path = os.path.join(self.galaxy_user_path, "galaxies.cat")
                assert os.path.isfile(old_galactic_catalog_path)

                # Open the 'old' galactic catalog
                old_galaxy_catalog = tables.from_file(old_galactic_catalog_path)

                # Create merged galactic catalog
                galaxy_catalog = catalogs.merge_galactic_catalogs(self.galaxy_extractor.catalog, old_galaxy_catalog)

                # Save the merged catalog
                path = os.path.join(self.galaxy_user_path, "galaxies.cat")
                tables.write(galaxy_catalog, path)


                # -- MERGE STELLAR CATALOG --

                # Check whether there is a stellar catalog file in the galaxy's directory
                old_stellar_catalog_path = os.path.join(self.galaxy_user_path, "stars.cat")
                assert os.path.isfile(old_stellar_catalog_path)

                # Open the 'old' stellar catalog
                old_stellar_catalog = tables.from_file(old_stellar_catalog_path)

                # Find the maximal DustPedia star ID
                maximal_id = -1
                for j in range(len(old_stellar_catalog)):
                    if old_stellar_catalog["Id"][j] > maximal_id:
                        maximal_id = int(old_stellar_catalog["Id"][j].split("/")[1])

                # Create merged stellar catalog
                self.stellar_catalog = self.create_merged_stellar_catalog(maximal_id, old_stellar_catalog, coverage)

                # Not good: wasted many DustPedia star id's for stars that then are not added after all because
                # they are matched to a star already in the old catalog
                #stellar_catalog = self.create_stellar_catalog(maximal_id)
                #stellar_catalog = catalogs.merge_stellar_catalogs(stellar_catalog, old_stellar_catalog)

                if self.config.write:

                    # Save the merged catalog
                    path = os.path.join(self.galaxy_user_path, "stars.cat")
                    tables.write(self.stellar_catalog, path)


                    # -- UPDATE RANGES TABLE --

                    coverage.add_box(coordinate_box)
                    coverage.save()

        else:

            # Create the directory to contain the catalogs for this galaxy
            filesystem.create_directory(self.galaxy_user_path)

            coverage.add_box(coordinate_box)

            # Save galactic catalog
            galactic_catalog_path = os.path.join(self.galaxy_user_path, "galaxies.cat")
            tables.write(self.galaxy_extractor.catalog, galactic_catalog_path)

            # Append stars from the trained extractor
            # Loop over the stars found by the trained extractor
            maximal_id = -1
            self.stellar_catalog = self.create_stellar_catalog(maximal_id)

            if self.config.write:

                # Save stellar catalog
                stellar_catalog_path = os.path.join(self.galaxy_user_path, "stars.cat")
                tables.write(self.stellar_catalog, stellar_catalog_path)

                # Save the coverage or range table
                coverage.save()

    # -----------------------------------------------------------------

    def create_merged_stellar_catalog(self, maximal_id, old_stellar_catalog, old_coverage):

        """
        This function ...
        :param maximal_id:
        :return:
        """

        catalog_column = []
        id_column = []
        ra_column = []
        dec_column = []
        ra_error_column = []
        dec_error_column = []
        confidence_level_column = []
        on_galaxy_column = []
        original_id_column = []

        # Keep track of which stars in catalog a have been encountered in the 'old' catalog (to avoid unnecessary checking)
        encountered = [False] * len(old_stellar_catalog)

        # Append stars from the star extractor
        # Loop over the stellar statistics
        for i in range(len(self.star_extractor.statistics)):

            id = self.galaxy_name.replace(' ', '') + "/" + str(maximal_id + 1)

            index = self.star_extractor.statistics["Star index"][i]

            # Skip undetected stars
            if not self.star_extractor.statistics["Detected"][i]: continue

            # IMPORTANT: Skip stars that are already in the cached catalog
            for j in range(len(old_stellar_catalog)):

                # Skip this entry (star) if it has already been matched with a star from the new catalog
                if encountered[j]: continue

                old_catalog_original_catalog_star = old_stellar_catalog["Original catalog and id"][j].split("**")[0]
                old_catalog_original_id_star = old_stellar_catalog["Original catalog and id"][j].split("**")[1]

                # Check whether it's a match
                if self.star_extractor.catalog["Catalog"][index] == old_catalog_original_catalog_star and self.star_extractor.catalog["Id"][index] == old_catalog_original_id_star:

                    # This star from the old catalog has been matched
                    encountered[j] = True
                    break

            # If a break is not encountered, no match is found -> add the new star to the catalog
            else:

                catalog_column.append("DustPedia")
                id_column.append(id)
                ra_column.append(self.star_extractor.catalog["Right ascension"][index])
                dec_column.append(self.star_extractor.catalog["Declination"][index])
                ra_error_column.append(self.star_extractor.catalog["Right ascension error"][index])
                dec_error_column.append(self.star_extractor.catalog["Declination error"][index])
                confidence_level_column.append(self.star_extractor.catalog["Confidence level"][index])
                on_galaxy_column.append(self.star_extractor.catalog["On galaxy"][index])
                original_id_column.append(self.star_extractor.catalog["Catalog"][index] + "**" + self.star_extractor.catalog["Id"][index])

                maximal_id += 1

        position_error = 0.5 * self.frame.pixelscale * 1000  # in mas !!

        # Append stars from the trained extractor
        # Loop over the stars found by the trained extractor
        for star in self.trained_extractor.stars:

            ra_deg = star.position.ra.value
            dec_deg = star.position.dec.value

            # Check if the 'other source' lies outside of what the old catalog currently covers
            position = Position(ra_deg, dec_deg)
            if not old_coverage.covers(position):

                id = self.galaxy_name.replace(' ', '') + "/" + str(maximal_id + 1)

                catalog_column.append("DustPedia")
                id_column.append(id)
                ra_column.append(ra_deg)
                dec_column.append(dec_deg)
                ra_error_column.append(position_error) # in mas
                dec_error_column.append(position_error)  # in mas
                confidence_level_column.append(star.confidence_level)
                on_galaxy_column.append(False)
                original_id_column.append(None)

                maximal_id += 1

        data = [catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, confidence_level_column,
                on_galaxy_column, original_id_column]
        names = ['Catalog', 'Id', 'Right ascension', 'Declination', 'Right ascension error', 'Declination error', 'Confidence level',
                 'On galaxy', 'Original catalog and id']
        stellar_catalog = tables.new(data, names)

        return stellar_catalog

    # -----------------------------------------------------------------

    def create_stellar_catalog(self, maximal_id):

        """
        This function ...
        :param maximal_id:
        :return:
        """

        catalog_column = []
        id_column = []
        ra_column = []
        dec_column = []
        ra_error_column = []
        dec_error_column = []
        confidence_level_column = []
        on_galaxy_column = []
        original_id_column = []

        # Append stars from the star extractor
        # Loop over the stellar statistics
        for i in range(len(self.star_extractor.statistics)):

            id = self.galaxy_name.replace(' ', '') + "/" + str(maximal_id + 1)

            index = self.star_extractor.statistics["Star index"][i]

            # Skip undetected stars
            if not self.star_extractor.statistics["Detected"][i]: continue

            catalog_column.append("DustPedia")
            id_column.append(id)
            ra_column.append(self.star_extractor.catalog["Right ascension"][index])
            dec_column.append(self.star_extractor.catalog["Declination"][index])
            ra_error_column.append(self.star_extractor.catalog["Right ascension error"][index])
            dec_error_column.append(self.star_extractor.catalog["Declination error"][index])
            confidence_level_column.append(self.star_extractor.catalog["Confidence level"][index])
            on_galaxy_column.append(self.star_extractor.catalog["On galaxy"][index])
            original_id_column.append(self.star_extractor.catalog["Catalog"][index] + "**" + self.star_extractor.catalog["Id"][index])

            maximal_id += 1

        position_error = 0.5 * self.frame.pixelscale * 1000  # in mas !!

        # Append stars from the trained extractor
        # Loop over the stars found by the trained extractor
        for star in self.trained_extractor.stars:

            id = self.galaxy_name.replace(' ', '') + "/" + str(maximal_id + 1)

            catalog_column.append("DustPedia")
            id_column.append(id)
            ra_column.append(star.position.ra.value)
            dec_column.append(star.position.dec.value)
            ra_error_column.append(position_error)
            dec_error_column.append(position_error)
            confidence_level_column.append(star.confidence_level)
            on_galaxy_column.append(False)
            original_id_column.append(None)

            maximal_id += 1

        data = [catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, confidence_level_column,
                on_galaxy_column, original_id_column]
        names = ['Catalog', 'Id', 'Right ascension', 'Declination', 'Right ascension error', 'Declination error', 'Confidence level',
                 'On galaxy', 'Original catalog and id']
        stellar_catalog = tables.new(data, names)

        return stellar_catalog

    # -----------------------------------------------------------------

    def extents_current_ranges(self, table, ra, dec, ra_span, dec_span):

        """
        This function ...
        :param table:
        :param box:
        :return:
        """

        # Create rectangle
        center = Position(ra, dec)
        radius = Extent(0.5 * ra_span, 0.5 * dec_span)
        box = Rectangle(center, radius)

        # List of boxes from ranges in the table
        boxes = []

        # Loop over all entries in the table
        for i in range(len(table)):

            entry_ra = table["Central right ascension"][i]
            entry_dec = table["Central declination"][i]
            entry_ra_span = table["Right ascension span"][i]
            entry_dec_span = table["Declination span"][i]

            # Check whether the box properties are the same (neglecting small roundoff errors)
            same_ra = np.isclose(entry_ra, ra)
            same_dec = np.isclose(entry_dec, dec)
            same_ra_span = np.isclose(entry_ra_span, ra_span)
            same_dec_span = entry_dec_span == dec_span

            ## First check for identical entry
            identical_box_already_processed = same_ra and same_dec and same_ra_span and same_dec_span

            # If a frame with identical ranges has already been processed before
            if identical_box_already_processed: return False

            entry_center = Position(entry_ra, entry_dec)
            entry_radius = Extent(0.5 * entry_ra_span, 0.5 * entry_dec_span)

            entry_box = Rectangle(entry_center, entry_radius)

            boxes.append(entry_box)

        # Loop over corners of the box
        for corner in box.corners:

            # Loop over all boxes
            for box in boxes:

                # If at least one of the boxes contains this corner, this corner is OK
                if box.contains(corner): break

            # If a break is not encountered: this corner is not covered by any of the boxes
            else: return True

        # If a break IS encountered for each corner, each corner is covered by one or more boxes, so this box
        # does not extent the range of these other boxes
        return False

# -----------------------------------------------------------------
