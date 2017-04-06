#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalog.extended Contains the ExtendedSourceCatalog class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ..basics.coordinate import SkyCoordinate
from ...core.tools.logging import log
from ..core.extendedsource import ExtendedSource

# -----------------------------------------------------------------

class ExtendedSourceCatalog(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ExtendedSourceCatalog, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Name", str, None, "Name of the galaxy")
        self.add_column_info("RA", float, "deg", "Right ascension")
        self.add_column_info("DEC", float, "deg", "Declination")
        self.add_column_info("Redshift", float, None, "Redshift")
        self.add_column_info("Type", str, None, "Galaxy type")
        self.add_column_info("Names", str, None, "Alternative names")
        self.add_column_info("Distance", float, "Mpc", "distance")
        self.add_column_info("Incl", float, "deg", "inclination")
        self.add_column_info("D25", float, "arcmin", "D25")
        self.add_column_info("Major", float, "arcmin", "Major axis length")
        self.add_column_info("Minor", float, "arcmin", "Minor axis length")
        self.add_column_info("Posangle", float, "deg", "Position angle")
        self.add_column_info("Principal", bool, None, "Is principal galaxy")
        self.add_column_info("Companions", str, None, "Companion galaxies")
        self.add_column_info("Parent", str, None, "Parent galaxy (is companion galaxy)")

    # -----------------------------------------------------------------

    def add_entry(self, name, ra, dec, z, galtype, alternative_names, distance, inclination, d25, major, minor,
                  posangle, principal, companions, parent):

        """
        This function ...
        :param name:
        :param ra:
        :param dec:
        :param z:
        :param galtype:
        :param alternative_names:
        :param distance:
        :param inclination:
        :param d25:
        :param major:
        :param minor:
        :param posangle:
        :param principal:
        :param companions:
        :param parent:
        :return:
        """

        values = [name, ra, dec, z, galtype, alternative_names, distance, inclination, d25, major, minor, posangle,
                  principal, companions, parent]
        self.add_row(values)

    # -----------------------------------------------------------------

    def get_position(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return SkyCoordinate(ra=self["RA"][index], dec=self["DEC"][index], unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def get_name(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Name", index)

    # -----------------------------------------------------------------

    def get_redshift(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Redshift", index)

    # -----------------------------------------------------------------

    def get_type(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Type", index)

    # -----------------------------------------------------------------

    def get_names(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        names = self.get_quantity("Names", index)
        if names is not None: names = names.split(",")
        return names

    # -----------------------------------------------------------------

    def get_distance(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Distance", index)

    # -----------------------------------------------------------------

    def get_inclination(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Incl", index)

    # -----------------------------------------------------------------

    def get_d25(self, index):

        """
        This function ...
        :return:
        """

        return self.get_quantity("D25", index)

    # -----------------------------------------------------------------

    def get_major(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Major", index)

    # -----------------------------------------------------------------

    def get_minor(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Minor", index)

    # -----------------------------------------------------------------

    def get_position_angle(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Posangle", index)

    # -----------------------------------------------------------------

    def is_principal(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Principal", index)

    # -----------------------------------------------------------------

    def get_companions(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        companions = self.get_quantity("Companions", index)
        if companions is not None: companions = companions.split(",")
        return companions

    # -----------------------------------------------------------------

    def get_parent(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Parent", index)

    # -----------------------------------------------------------------

    def create_source(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Debugging
        log.debug("Creating an extended source for entry " + str(index) + " in the catalog ...")

        # Get the row
        row = self.get_row(index)

        # Get properties
        name = row["Name"]
        redshift = row["Redshift"]
        galaxy_type = row["Type"]
        position = self.get_position(index)
        names = row["Names"].split(",") if row["Names"] is not None else None
        distance = row["Distance"]
        inclination = row["Incl"]
        d25 = row["D25"]
        major = row["Major"]
        minor = row["Minor"]
        position_angle = row["Posangle"]
        principal = row["Principal"]
        companions = row["Companions"].split(",") if row["Companions"] is not None else None
        parent = row["Parent"]

        # Create a new ExtendedSource instance
        source = ExtendedSource(index=index, name=name, position=position, redshift=redshift, galaxy_type=galaxy_type,
                                names=names, distance=distance, inclination=inclination, d25=d25, major=major,
                                minor=minor, position_angle=position_angle, principal=principal, companions=companions,
                                parent=parent)

        # Return the source
        return source

# -----------------------------------------------------------------
