#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalogsynchronizer Contains the CatalogSynchronizer class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np

# Import the relevant PTS classes and modules
from ..basics.vector import Position, Extent
from ..region.rectangle import PixelRectangleRegion
from ..basics.catalogcoverage import CatalogCoverage
from ..tools import catalogs
from ...core.basics.configurable import Configurable
from ...core.tools import introspection, tables
from ...core.tools import filesystem as fs
from ...core.units.parsing import parse_unit

# -----------------------------------------------------------------

class CatalogSynchronizer(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(CatalogSynchronizer, self).__init__(*args, **kwargs)

        # The image frame
        self.frame = None

        # The name of the principal galaxy
        self.galaxy_name = None
        self.galaxy_user_path = None

        # The galactic and stellar catalog
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def run(self, frame, galaxy_name, galactic_catalog, stellar_catalog):

        """
        This function ...
        :param frame:
        :param galaxy_name:
        :param galactic_catalog:
        :param stellar_catalog:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galaxy_name, galactic_catalog, stellar_catalog)
        
        # 2. Synchronize the catalog
        self.synchronize()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set attributes to None
        self.frame = None
        self.galaxy_name = None
        self.galaxy_user_path = None
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def setup(self, frame, galaxy_name, galactic_catalog, stellar_catalog):

        """
        This function ...
        :param frame:
        :param galaxy_name:
        :param galactic_catalog:
        :param stellar_catalog:
        :return:
        """

        # Call the setup function of the base class
        super(CatalogSynchronizer, self).setup()

        # Set the image frame
        self.frame = frame

        # Name of the principal galaxy
        self.galaxy_name = galaxy_name

        # Determine the path to the user catalogs directory
        catalogs_user_path = os.path.join(introspection.pts_user_dir, "magic", "catalogs")

        # Determine the path to the directory to contain the catalogs for this galaxy
        self.galaxy_user_path = os.path.join(catalogs_user_path, self.galaxy_name)

        # Set the catalogs
        self.galactic_catalog = galactic_catalog
        self.stellar_catalog = stellar_catalog

    # -----------------------------------------------------------------

    def synchronize(self):

        """
        This function ...
        :return:
        """

        # Get bounding box
        coordinate_box = self.frame.bounding_box

        # Get current catalog coverage
        coverage = CatalogCoverage(self.galaxy_name)

        # Check if a directory exists for this galaxy
        if fs.is_directory(self.galaxy_user_path):

            # Add the coordinate box to the catalog coverage object
            if not coverage.covers(coordinate_box): coverage.add_box(coordinate_box)

            # Update the catalogs
            self.update_galactic_catalog()
            self.update_stellar_catalog()

        # If a directory does not exist yet for this galaxy
        else:

            # Create a directory for this galaxy
            fs.create_directory(self.galaxy_user_path)

            # Add the coordinate box to the catalog coverage object
            coverage.add_box(coordinate_box)

            # Create the catalogs
            self.create_new_galactic_catalog()
            self.create_new_stellar_catalog()

        # Update the catalog coverage
        coverage.save()

    # -----------------------------------------------------------------

    def create_new_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        # Write the galactic catalog to the appropriate location
        galactic_catalog_path = os.path.join(self.galaxy_user_path, "galaxies.cat")
        tables.write(self.galactic_catalog, galactic_catalog_path)

    # -----------------------------------------------------------------

    def create_new_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        last_id = -1

        # Loop over all entries in the stellar catalog
        for i in range(len(self.stellar_catalog)):

            # Skip stars already in the DustPedia catalog ## NOT POSSIBLE HERE BECAUSE DIRECTORY DID NOT EXIST WHEN WE ARE HERE
            #if self.stellar_catalog["Catalog"][i] == "DustPedia": continue

            # Set the original ID, if the 'Catalog' entry is not masked
            if not self.stellar_catalog["Catalog"].mask[i]:

                original_id = self.stellar_catalog["Catalog"][i] + "//" + self.stellar_catalog["Id"][i]
                self.stellar_catalog["Original catalog and id"][i] = original_id

            # Create a new Id
            self.stellar_catalog["Catalog"][i] = "DustPedia"
            self.stellar_catalog["Id"][i] = self.galaxy_name.replace(' ', '') + "/" + str(last_id + 1)

            # Increment the last Id
            last_id += 1

        # Write the stellar catalog to the appropriate location
        stellar_catalog_path = os.path.join(self.galaxy_user_path, "stars.cat")
        tables.write(self.stellar_catalog, stellar_catalog_path)

    # -----------------------------------------------------------------

    def update_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        # Check whether there is a galactic catalog file in the galaxy's directory
        old_galactic_catalog_path = os.path.join(self.galaxy_user_path, "galaxies.cat")
        assert os.path.isfile(old_galactic_catalog_path)

        # Open the 'old' galactic catalog
        old_galaxy_catalog = tables.from_file(old_galactic_catalog_path)

        # Create merged galactic catalog
        self.galactic_catalog = catalogs.merge_galactic_catalogs(self.galactic_catalog, old_galaxy_catalog)

        # Write the galactic catalog to the appropriate location
        galactic_catalog_path = os.path.join(self.galaxy_user_path, "galaxies.cat")
        tables.write(self.galactic_catalog, galactic_catalog_path)

    # -----------------------------------------------------------------

    def update_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        # Check whether there is a stellar catalog file in the galaxy's directory
        old_stellar_catalog_path = os.path.join(self.galaxy_user_path, "stars.cat")
        assert os.path.isfile(old_stellar_catalog_path)

        # Open the 'old' stellar catalog
        old_stellar_catalog = tables.from_file(old_stellar_catalog_path)

        # Find the maximal DustPedia star ID
        last_id = -1
        for j in range(len(old_stellar_catalog)):
            entry_id = int(old_stellar_catalog["Id"][j].split("/")[1])
            if entry_id > last_id: last_id = entry_id

        encountered = [False] * len(old_stellar_catalog)

        import copy
        new_stellar_catalog = copy.deepcopy(old_stellar_catalog)

        # Loop over all entries in the stellar catalog
        for i in range(len(self.stellar_catalog)):

            # Skip stars already in the DustPedia catalog
            if self.stellar_catalog["Catalog"][i] == "DustPedia": continue
            if self.matches_dustpedia_star(i, old_stellar_catalog, encountered): continue

            # Set the original ID, if the 'Catalog' entry is not masked
            if not self.stellar_catalog["Catalog"].mask[i]:

                original_id = self.stellar_catalog["Catalog"][i] + "//" + self.stellar_catalog["Id"][i]
                self.stellar_catalog["Original catalog and id"] = original_id

            # Create a new Id
            self.stellar_catalog["Catalog"] = "DustPedia"
            self.stellar_catalog["Id"] = self.galaxy_name.replace(' ', '') + "/" + str(last_id + 1)

            # Add the star to the old catalog
            new_stellar_catalog.add_row([self.stellar_catalog["Catalog"], self.stellar_catalog["Id"], self.stellar_catalog["Right ascension"], self.stellar_catalog["Declination"],
                                         self.stellar_catalog["Right ascension error"], self.stellar_catalog["Declination error"], self.stellar_catalog["Confidence level"],
                                         self.stellar_catalog["On galaxy"], self.stellar_catalog["Original catalog and id"]])

            # Increment the last Id
            last_id += 1

        # Write the stellar catalog to the appropriate location
        stellar_catalog_path = os.path.join(self.galaxy_user_path, "stars.cat")
        tables.write(new_stellar_catalog, stellar_catalog_path)

    # -----------------------------------------------------------------

    def matches_dustpedia_star(self, i, old_stellar_catalog, encountered):

        """
        This function ...
        :param i
        :param old_stellar_catalog:
        :param encountered:
        :return:
        """

        # If an original catalog exists for this star
        if not self.stellar_catalog["Catalog"].mask[i]:

            catalog_name = self.stellar_catalog["Catalog"][i]
            star_id = self.stellar_catalog["Id"][i]

            # Skip stars already in the DustPedia catalog
            for j in range(len(old_stellar_catalog)):

                if encountered[j]: continue
                if old_stellar_catalog["Original catalog and id"].mask[j]: continue
                if old_stellar_catalog["Original catalog and id"][j] == catalog_name + "//" + star_id:
                    encountered[j] = True
                    return True

        # Else (stars without original catalog), cross-reference with old catalog based on position
        else:

            # Loop over the stars in the DustPedia catalog who have no original catalog either
            for j in range(len(old_stellar_catalog)):

                if encountered[j]: continue
                if not old_stellar_catalog["Original catalog and id"].mask[j]: continue

                # Calculate the distance between the star in the old catalog and the new star
                old_position = Position(old_stellar_catalog["Right ascension"][j], old_stellar_catalog["Declination"][j])
                new_position = Position(self.stellar_catalog["Right ascension"][i], self.stellar_catalog["Declination"][i])
                distance = (old_position - new_position).norm

                # Calculate the error on the star's position in the old catalog
                old_error_ra_mas = old_stellar_catalog["Right ascension error"][j] * parse_unit("mas")
                old_error_ra_deg = old_error_ra_mas.to("deg").value
                old_error_dec_mas = old_stellar_catalog["Declination error"][j] * parse_unit("mas")
                old_error_dec_deg = old_error_dec_mas.to("deg").value
                old_error = Extent(old_error_ra_deg, old_error_dec_deg).norm

                # Calculate the error on the star's position in the new catalog
                new_error_ra_mas = self.stellar_catalog["Right ascension error"][i] * parse_unit("mas")
                new_error_ra_deg = new_error_ra_mas.to("deg").value
                new_error_dec_mas = self.stellar_catalog["Declination error"][i] * parse_unit("mas")
                new_error_dec_deg = new_error_dec_mas.to("deg").value
                new_error = Extent(new_error_ra_deg, new_error_dec_deg).norm

                # If the distance is smaller than the error on the old star's position and the error on the
                # new star's position, we have a match
                if distance < old_error and distance < new_error:
                    encountered[j] = True
                    return True

        # If no match is found, return False
        return False

    # -----------------------------------------------------------------

    def old_synchronize(self):

        """
        This function ...
        :return:
        """

        # Get bounding box
        coordinate_box = self.frame.bounding_box

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

                if self.config.write:

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

            # The coordinate range of the current frame does not extent that of previous frames
            #else:

                #maximal_id = -1
                #self.create_stellar_catalog(maximal_id)

        else:

            # Create the directory to contain the catalogs for this galaxy
            if self.config.write:

                fs.create_directory(self.galaxy_user_path)

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
        box = PixelRectangleRegion(center, radius)

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

            entry_box = PixelRectangleRegion(entry_center, entry_radius)

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
