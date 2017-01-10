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

    column_info = [("Name", str, None, "Name of the galaxy"),
                    ("RA", float, "deg", "Right ascension"),
                    ("DEC", float, "deg", "Declination"),
                    ("Redshift", float, None, "Redshift"),
                    ("Type", str, None, "Galaxy type"),
                    ("Names", str, None, "Alternative names"),
                    ("Distance", float, "Mpc", "distance"),
                    ("Incl", float, "deg", "inclination"),
                    ("D25", float, "arcmin", "D25"),
                    ("Major", float, "arcmin", "Major axis length"),
                    ("Minor", float, "arcmin", "Minor axis length"),
                    ("Posangle", float, "deg", "Position angle"),
                    ("Principal", bool, None, "Is principal galaxy"),
                    ("Companions", str, None, "Companion galaxies"),
                    ("Parent", str, None, "Parent galaxy (is companion galaxy)")]

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
        names = row["Nmes"].split(",")
        distance = row["Distance"]
        inclination = row["Incl"]
        d25 = row["D25"]
        major = row["Major"]
        minor = row["Minor"]
        position_angle = row["Posangle"]
        principal = row["Principal"]
        companions = row["Companions"].split(",")
        parent = row["Parent"]

        # Create a new ExtendedSource instance
        source = ExtendedSource(index=index, name=name, position=position, redshift=redshift, galaxy_type=galaxy_type,
                                names=names, distance=distance, inclination=inclination, d25=d25, major=major,
                                minor=minor, position_angle=position_angle, principal=principal, companions=companions,
                                parent=parent)

        # Return the source
        return source

# -----------------------------------------------------------------
