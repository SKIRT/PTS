#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.sample Contains the DustPediaSample class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...magic.basics.coordinate import SkyCoordinate
from ...magic.tools import catalogs
from ...core.units.parsing import parse_unit as u
from ...core.tools import tables

# -----------------------------------------------------------------

dustpedia_dat_path = fs.join(introspection.pts_dat_dir("dustpedia"))

# -----------------------------------------------------------------

dustpedia_data_path = fs.join(dustpedia_dat_path, "data")

# -----------------------------------------------------------------

ledawise_table_path = fs.join(dustpedia_data_path, "DustPedia_LEDAWISE_Herschel.csv")

# -----------------------------------------------------------------

filter_names = {"GALEX FUV": "GALEX_FUV",
                "GALEX NUV": "GALEX_NUV",
                "SDSS u": "SDSS_u",
                "SDSS g": "SDSS_g",
                "SDSS r": "SDSS_r",
                "SDSS i": "SDSS_i",
                "SDSS z": "SDSS_z",
                "2MASS J": "2MASS_J",
                "2MASS H": "2MASS_H",
                "2MASS Ks": "2MASS_Ks",
                "IRAC I1": "Spitzer_3.6",
                "IRAC I2": "Spitzer_4.5",
                "IRAC I3": "Spitzer_5.8",
                "IRAC I4": "Spitzer_8.0",
                "MIPS 24mu": "Spitzer_24",
                "MIPS 70mu": "Spitzer_70",
                "MIPS 160mu": "Spitzer_160",
                "WISE W1": "WISE_3.4",
                "WISE W2": "WISE_4.6",
                "WISE W3": "WISE_12",
                "WISE W4": "WISE_22",
                "Pacs blue": "PACS_70",
                "Pacs green": "PACS_100",
                "Pacs red": "PACS_160",
                "SPIRE PSW": "SPIRE_250",
                "SPIRE PMW": "SPIRE_350",
                "SPIRE PLW": "SPIRE_500",
                "HFI 857": "Planck_350",
                "HFI 545": "Planck_550",
                "HFI 353": "Planck_850",
                "HFI 217": "Planck_1380",
                "HFI 143": "Planck_2100",
                "HFI_100": "Planck_3000",
                "LFI 070": "Planck_4260",
                "LFI 044": "Planck_6810",
                "LFI 030": "Planck_10600"}

# -----------------------------------------------------------------

# naming convention: [galaxy] [telescope] [band].fits.

# -----------------------------------------------------------------

def resolve_name(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    sample = DustPediaSample()
    name = sample.get_name(galaxy_name)
    return name

# -----------------------------------------------------------------

properties_table_path = fs.join(dustpedia_data_path, "DustPedia_HyperLEDA_Herschel.csv")

# -----------------------------------------------------------------

def get_center(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    sample = DustPediaSample()
    return sample.get_position(galaxy_name)

# -----------------------------------------------------------------

def get_distance(galaxy_name):

    """
    This fucntion ...
    :param galaxy_name:
    :return:
    """

    name = resolve_name(galaxy_name)

    # Load the table
    properties = Table.read(properties_table_path)

    # Find galaxy
    index = tables.find_index(properties, name)

    # Return the distance
    return properties["dist_best"][index] * u("Mpc")

# -----------------------------------------------------------------

def get_inclination(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    name = resolve_name(galaxy_name)

    # Load the table
    properties = Table.read(properties_table_path)

    # Find galaxy
    index = tables.find_index(properties, name)

    # Return the inclination
    return Angle(properties["incl"][index], "deg")

# -----------------------------------------------------------------

class DustPediaSample(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Load the table
        leda_wise_herschel = Table.read(ledawise_table_path) # 877 galaxies

        # Get the object names
        names = list(leda_wise_herschel["objname"])
        self.primary_sample = names

        # The Vizier querying object
        self.vizier = Vizier()
        self.vizier.ROW_LIMIT = -1

    # -----------------------------------------------------------------

    def get_names(self):

        """
        This function ...
        :return:
        """

        return self.primary_sample

    # -----------------------------------------------------------------

    def get_position(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        result = self.vizier.query_object(galaxy_name, catalog=["VII/237"])
        table = result[0]

        if len(table) > 1: raise ValueError("Ambiguous result")
        if len(table) == 0: raise ValueError("No result")

        # Get coordinate
        ra, dec = catalogs.get_ra_dec_degrees_from_table(table, 0)

        # Create sky coordinate and return it
        return SkyCoordinate(ra=ra, dec=dec, unit="deg")

    # -----------------------------------------------------------------

    def get_d25(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        result = self.vizier.query_object(galaxy_name, catalog=["VII/237"])
        table = result[0]

        if len(table) > 1: raise ValueError("Ambiguous result")

        # Calculate the diameter
        diameter = np.power(10.0, table["logD25"][0]) * 0.1 * u("arcmin") if table["logD25"][0] else None

        # Return
        return diameter

    # -----------------------------------------------------------------

    def get_r25(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        result = self.vizier.query_object(galaxy_name, catalog=["VII/237"])
        table = result[0]

        if len(table) > 1: raise ValueError("Ambiguous result")

        # Calculate the ratio
        ratio = np.power(10.0, table["logR25"][0]) if table["logR25"][0] else None

        # Return the ratio
        return ratio

    # -----------------------------------------------------------------

    def get_name(self, galaxy_name): # gets name = HYPERLEDA name, and checks whether in DustPedia sample

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the HYPERLEDA name
        objname = catalogs.get_hyperleda_name(galaxy_name)

        # CHECK WHETHER IN DUSTPEDIA SAMPLE
        if objname not in self.primary_sample: raise ValueError("This galaxy is not in the DustPedia primary sample")

        # Return the name of the galaxy
        return objname

    # -----------------------------------------------------------------

    def get_filter_name(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        fltr_string = str(fltr)

        if fltr_string not in filter_names: raise ValueError("Filter not recognized")

        return filter_names[fltr_string]

# -----------------------------------------------------------------
