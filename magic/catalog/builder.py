#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalogbuilder Contains the CatalogBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import tables

# -----------------------------------------------------------------

class CatalogBuilder(Configurable):

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
        super(CatalogBuilder, self).__init__(*args, **kwargs)

        # The image frame
        self.frame = None

        # References to the extractors
        self.galaxy_extractor = None
        self.star_extractor = None
        self.trained_extractor = None

        # The output catalogs
        self.galactic_catalog = None
        self.stellar_catalog = None

    # -----------------------------------------------------------------

    def run(self, frame, galaxy_extractor, star_extractor, trained_extractor):

        """
        This function ...
        :param frame:
        :param galaxy_extractor:
        :param star_extractor:
        :param trained_extractor:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galaxy_extractor, star_extractor, trained_extractor)
        
        # 2. Build the catalog
        self.build()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set attributes to None
        self.frame = None
        self.galaxy_extractor = None
        self.star_extractor = None
        self.trained_extractor = None

    # -----------------------------------------------------------------

    def setup(self, frame, galaxy_extractor, star_extractor, trained_extractor):

        """
        This function ...
        :param frame:
        :param galaxy_extractor:
        :param star_extractor:
        :param trained_extractor:
        :return:
        """

        # Call the setup function of the base class
        super(CatalogBuilder, self).setup()

        # The frame
        self.frame = frame

        # References to the extractors
        self.galaxy_extractor = galaxy_extractor
        self.star_extractor = star_extractor
        self.trained_extractor = trained_extractor

    # -----------------------------------------------------------------

    def build(self):

        """
        This function ...
        :return:
        """

        # Build the galactic catalog
        self.build_galactic_catalog()

        # Build the stellar catalog
        self.build_stellar_catalog()

    # -----------------------------------------------------------------

    def build_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        # Set galactic catalog (no merging with trained extractor (yet) and undetected galaxies are included anyway)
        self.galactic_catalog = self.galaxy_extractor.catalog

    # -----------------------------------------------------------------

    def build_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        # Initialize columns
        catalog_column = []
        id_column = []
        ra_column = []
        dec_column = []
        ra_error_column = []
        dec_error_column = []
        confidence_level_column = []
        on_galaxy_column = []
        original_id_column = []

        # Append stars from the star extractor; loop over the stellar statistics
        for i in range(len(self.star_extractor.statistics)):

            # Get the index of this star in the input catalog used by the star extractor
            index = self.star_extractor.statistics["Star index"][i]

            # Skip undetected stars
            if not self.star_extractor.statistics["Detected"][i]: continue

            # Add the appropriate values in the columns
            catalog_column.append(self.star_extractor.catalog["Catalog"][index] if not (hasattr(self.star_extractor.catalog["Catalog"], "mask") and self.star_extractor.catalog["Catalog"].mask[index]) else None)
            id_column.append(self.star_extractor.catalog["Id"][index] if not (hasattr(self.star_extractor.catalog["Id"], "mask") and self.star_extractor.catalog["Id"].mask[index]) else None)
            ra_column.append(self.star_extractor.catalog["Right ascension"][index])
            dec_column.append(self.star_extractor.catalog["Declination"][index])
            ra_error_column.append(self.star_extractor.catalog["Right ascension error"][index])
            dec_error_column.append(self.star_extractor.catalog["Declination error"][index])
            confidence_level_column.append(self.star_extractor.catalog["Confidence level"][index])
            on_galaxy_column.append(self.star_extractor.catalog["On galaxy"][index])
            original_id_column.append(None)

        #position_error = 0.5 * self.frame.average_pixelscale.to("mas").value  # in mas !!
        x_position_error = 0.5 * self.frame.pixelscale.x.to("mas").value
        y_position_error = 0.5 * self.frame.pixelscale.y.to("mas").value

        # Append stars from the trained extractor; loop over the stars found by the trained extractor
        for star in self.trained_extractor.stars:

            # Add the appropriate values in the columns
            catalog_column.append(None)
            id_column.append(None)
            ra_column.append(star.position.ra.value)
            dec_column.append(star.position.dec.value)
            ra_error_column.append(x_position_error)
            dec_error_column.append(y_position_error)
            confidence_level_column.append(star.confidence_level)
            on_galaxy_column.append(False)
            original_id_column.append(None)

        data = [catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, confidence_level_column,
                on_galaxy_column, original_id_column]
        names = ['Catalog', 'Id', 'Right ascension', 'Declination', 'Right ascension error', 'Declination error', 'Confidence level',
                 'On galaxy', 'Original catalog and id']

        # Create the merged stellar catalog
        self.stellar_catalog = tables.new(data, names)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write the galactic catalog
        self.write_galactic_catalog()

        # Write the stellar catalog
        self.write_stellar_catalog()

    # -----------------------------------------------------------------

    def write_galactic_catalog(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_stellar_catalog(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
