#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalog.point Contains the PointSourceCatalog class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.tools.logging import log
from ..core.pointsource import PointSource
from ..basics.coordinate import SkyCoordinate

# -----------------------------------------------------------------

class PointSourceCatalog(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PointSourceCatalog, self).__init__(*args, **kwargs)

        # Add columns
        self.add_column_info("Catalog", str, None, "Original catalog")
        self.add_column_info("ID", str, None, "ID in the catalog")
        self.add_column_info("RA", float, "deg", "Right ascension")
        self.add_column_info("DEC", float, "deg", "Declination")
        self.add_column_info("RA error", float, "mas", "Error on right ascension")
        self.add_column_info("DEC error", float, "mas", "Error on declination")
        self.add_column_info("Confidence", int, None, "Confidence level")

    # -----------------------------------------------------------------

    def add_coordinate(self, coordinate, catalog=None, id=None, ra_error=None, dec_error=None):

        """
        This function ...
        :param coordinate:
        :param catalog:
        :param id:
        :param ra_error:
        :param dec_error:
        :return:
        """

        self.add_entry(catalog, id, coordinate.ra, coordinate.dec, ra_error, dec_error, confidence_level=None)

    # -----------------------------------------------------------------

    def add_entry(self, catalog, id, ra, dec, ra_error, dec_error, confidence_level):

        """
        This function ...
        :param catalog:
        :param id:
        :param ra:
        :param dec:
        :param ra_error:
        :param dec_error:
        :param confidence_level
        :return:
        """

        # Add
        values = [catalog, id, ra, dec, ra_error, dec_error, confidence_level]
        self.add_row(values)

    # -----------------------------------------------------------------

    def coordinates(self):

        """
        This function ...
        :return:
        """

        coordinates = []

        unit = self.column_unit("RA")
        assert unit == self.column_unit("DEC")

        for index in range(len(self)):

            ra = self["RA"][index]
            dec = self["DEC"][index]

            coordinate = SkyCoordinate(ra=ra, dec=dec, unit=unit)
            coordinates.append(coordinate)

        # Return the coordinates
        return coordinates

    # -----------------------------------------------------------------

    def get_position(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return SkyCoordinate(ra=self["RA"][index], dec=self["DEC"][index], unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def get_catalog(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Catalog", index)

    # -----------------------------------------------------------------

    def get_id(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("ID", index)

    # -----------------------------------------------------------------

    def get_ra_error(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("RA error", index)

    # -----------------------------------------------------------------

    def get_dec_error(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("DEC error", index)

    # -----------------------------------------------------------------

    def get_confidence(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Confidence", index)

    # -----------------------------------------------------------------

    def create_source(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Debugging
        log.debug("Creating a point source for entry " + str(index) + " in the catalog ...")

        #print(index, len(self))

        # Get the row
        row = self.get_row(index)

        #print(row)

        # Get properties
        catalog = row["Catalog"]
        id = row["ID"]
        position = self.get_position(index)
        ra_error = row["RA error"]
        dec_error = row["DEC error"]
        confidence = row["Confidence"]

        # Check for which bands magnitudes are defined
        magnitudes = dict()
        magnitude_errors = {}
        for name in row:
            if "magnitude" in name:
                band = name.split(" magnitude")[0]
                magnitudes[band] = row[name]
                magnitude_errors[band] = row[name + " error"]

        # Create a new PointSource instance
        source = PointSource(index=index, position=position, catalog=catalog, id=id, ra_error=ra_error,
                             dec_error=dec_error, confidence=confidence, magnitudes=magnitudes,
                             magnitude_errors=magnitude_errors)

        # Return the source
        return source

# -----------------------------------------------------------------
