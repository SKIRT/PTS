#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.dataprocessing Contains the DustPediaDataProcessing class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import tempfile

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.tools import tables
from ...core.tools import time
from ...magic.tools import mosaicing
from ...core.units.parsing import parse_unit as u
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

dustpedia_dat_path = fs.join(introspection.pts_dat_dir("dustpedia"))

# -----------------------------------------------------------------

dustpedia_data_path = fs.join(dustpedia_dat_path, "data")

# -----------------------------------------------------------------

# Paths to Chris' tables
galex_url_table_path = fs.join(dustpedia_data_path, "GALEX_DustPedia_Herschel_Tile_URLs.dat")
galex_observations_table_path = fs.join(dustpedia_data_path, "DustPedia_Herschel_GALEX_Results.csv")
sdss_fields_table_path = fs.join(dustpedia_data_path, "SDSS_DR12_Primary_Fields.dat")
ledawise_table_path = fs.join(dustpedia_data_path, "DustPedia_LEDAWISE_Herschel.csv")

# -----------------------------------------------------------------

dustpedia_final_pixelsizes = {"GALEX": 3.2 * u("arcsec"), "SDSS": 0.45 * u("arcsec"),
                              "2MASS": 1. * u("arcsec"), "WISE": 1.375 * u("arcsec")} # in arcsec, 2MASS and WISE pixel sizes are original pixelsizes

# -----------------------------------------------------------------

class DustPediaDataProcessing(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Determine the path to a temporary directory
        self.temp_path = fs.join(tempfile.gettempdir(), time.unique_name("DustPedia"))
        fs.create_directory(self.temp_path)

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        if fs.is_directory(self.temp_path): fs.remove_directory(self.temp_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def galex_observations_table(self):

        """
        This property ...
        :return:
        """

        return get_galex_observations_table()

    # -----------------------------------------------------------------

    @lazyproperty
    def galex_url_table(self):

        """
        This property ...
        :return:
        """

        return get_galex_url_table()

    # -----------------------------------------------------------------

    @lazyproperty
    def leda_wise_table(self):

        """
        This function ...
        :return:
        """

        return get_leda_wise_table()

    # -----------------------------------------------------------------

    @lazyproperty
    def sdss_field_table(self):

        """
        This function ...
        :return:
        """

        return get_sdss_field_table()

    # -----------------------------------------------------------------

    def get_pixelscale_for_instrument(self, instrument):

        """
        This function ...
        :param instrument:
        :return:
        """

        return dustpedia_final_pixelsizes[instrument]

    # -----------------------------------------------------------------

    def get_header_for_galaxy(self, galaxy_name, instrument, returns="header"):

        """
        This function ...
        :param galaxy_name:
        :param instrument
        :param returns:
        :return:
        """

        # Get cutout range
        ra, dec, width = self.get_cutout_range_for_galaxy(galaxy_name)

        # Get the pixelscale
        pixelscale = dustpedia_final_pixelsizes[instrument]

        # Get the header
        return mosaicing.make_header(ra, dec, width, pixelscale, returns=returns)

    # -----------------------------------------------------------------

    def get_cutout_range_for_galaxy(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Find the index of the galaxy in the LEDA - WISE table
        index = tables.find_index(self.leda_wise_table, galaxy_name)

        # Get the RA and DEC
        ra = self.leda_wise_table["ra2000"][index] * u("deg")
        dec = self.leda_wise_table["de2000"][index] * u("deg")

        # Get the D25
        d25_arcmin = self.leda_wise_table["d25"][index]

        if d25_arcmin < 6: width = .5 * u("deg")
        else: width = 1. * u("deg")

        # Return the RA, DEC and width (all in degrees)
        return ra, dec, width

# -----------------------------------------------------------------

def get_leda_wise_table():

    """
    This function ...
    :return:
    """

    # Load the table
    table = Table.read(ledawise_table_path)

    # Return the table
    return table

# -----------------------------------------------------------------

def get_sdss_field_table():

    """
    This function ...
    :return:
    """

    # Load the table
    table = Table.read(sdss_fields_table_path, format="ascii.commented_header")

    # Return the table
    return table

# -----------------------------------------------------------------

def get_galex_url_table():

    """
    This function ...
    :return:
    """

    # Load the table
    table = Table.read(galex_url_table_path, format="ascii.no_header")

    # Set column name
    table.rename_column("col1", "URL")

    # Return the table
    return table

# -----------------------------------------------------------------

def get_galex_observations_table():

    """
    This function ...
    :return:
    """

    # Load the table
    table = Table.read(galex_observations_table_path)

    # Return the table
    return table

# -----------------------------------------------------------------
