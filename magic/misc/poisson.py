#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.poisson Contains the PoissonErrorCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import urllib
import tempfile

# Import astronomical modules
import montage_wrapper as montage
from astropy.table import Table
from astropy.utils import lazyproperty
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import inspection
from ...core.tools import filesystem as fs
from ...core.tools import tables
from ...magic.tools import plotting

# -----------------------------------------------------------------

dustpedia_dat_path = fs.join(inspection.pts_dat_dir("magic"), "dustpedia")

# -----------------------------------------------------------------

# Paths to Chris' tables
galex_url_table_path = fs.join(dustpedia_dat_path, "GALEX_DustPedia_Herschel_Tile_URLs.dat")
galex_observations_table_path = fs.join(dustpedia_dat_path, "DustPedia_Herschel_GALEX_Results.csv")
sdss_fields_table_path = fs.join(dustpedia_dat_path, "SDSS_DR12_Primary_Fields.dat")
ledawise_table_path = fs.join(dustpedia_dat_path, "DustPedia_LEDAWISE_Herschel.csv")

# -----------------------------------------------------------------

dustpedia_final_pixelsizes = {"GALEX": 3.2 * Unit("arcsec"), "SDSS": 0.45 * Unit("arcsec"),
                              "2MASS": 1. * Unit("arcsec"), "WISE": 1.375 * Unit("arcsec")} # in arcsec, 2MASS and WISE pixel sizes are original pixelsizes

# -----------------------------------------------------------------

# GALEX:

# For GALEX, I should be able to give you exactly what you're after immediately.
# I have attached a file called DustPedia_Herschel_GALEX_Results.csv containing the details of every GALEX
# observation that I used in the cutout-making process (such as exposure time), and which target galaxy each
# observation was used for (ie, a given observation can appear more than once in that file, as it could be associated
# with more than one target galaxy - and vice-a-versa).

# However, note that some of these observations are not incorporated into the final cutouts, for various reasons.
# For example:
# - Any observations that were not contiguously "attached" to the location of the target galaxy (ie, if there were
#   "gaps" in coverage between observations) were rejected.
# - The outer region of each observation was masked, to remove low-quality data. In some cases, this resulted in the
#   remainder of the observation not providing coverage of the region of interest.
# - The final cutouts had diameters of either 1 degree of 0.5 degrees. However, for all sources I queried a region of
#   sky large enough to produce a 1 degree cutout - I just discarded the observation that were not needed for the smaller
#   cutouts.

# Also, note that the DustPedia GALEX cutouts have been re-gridded to larger pixel sizes than the standard GALEX
# archive data (the ancillary data report provides details).

# Also attached is a file that lists the URLs each GALEX observation can be downloaded from. Annoyingly, the URLs are
# *not* properly standardised; simply knowing the details of am observation is *not* enough to automatically construct
# the correct URL! I had to write a script that crawled through the GALEX GR6/7 archive website page-by-page,
# doing string-matching to find the appropriate URLs... I hate that archive.

# It's easy to work out which URLs correspond to which observations, as the URL for a given observation always
# contains the tilename (note that there will often be multiple observations, and hence URLS, associated with a given tilename).


#SDSS:

# For SDSS, you can use the Montage function mArchiveGet to get a listing of the SDSS DR9 fields that cover a given
# area of sky (there is a nice Python wrapper available for Montage, which can make it easier to interact with).
# This function outputs a table giving the details of the relevant fields, including URLs where they can be downloaded.

#Note that I only used the SDSS primary fields to produce the DustPedia cutouts (the Ancillary Data Report I attached in the previous email explains this in more detail). I have attached a file that lists all the SDSS primary fields. After you have the results table produced by mArchiveGet, you can match the results with the primary fields list to find out the fields I used to make the final cutouts.

#Also, bear in mind that the DustPedia SDSS cutouts are re-gridded to North-East orientation, and 0.45" pixel sizes (again, see the ancillary data report for details).

#The SDSS database website provides extensive information about every SDSS field. It may be possible to use this database to find out the information you want about each field without having to download the actual FITS files, once you have the fields' ID information.


## EXTRA INFO:


# It should be straightforward to produce things to the same pixel grid as I used. I used the Montage command mHdr to
# construct the basic header, including WCS, for my cutouts. And in my experience, giving mHdr a particular set of
# inputs always results in the same output. Specifically, I ran mHdr through the Montage wrapper for Python, as follows:

# montage_wrapper.commands.mHdr( str(ra)+' '+str(dec), width, '/some/output/path/header.hdr', pix_size=pix_size )

# Where:
# - ra and dec are taken from the ra2000 and de2000 of the attached DustPedia_LEDAWISE_Herschel.csv table
# - width is 0.5 degrees for galaxies with D25<6 arcmin, and 1 degree for galaxies with D25>=6 arcmin (as listed in DustPedia_LEDAWISE_Herschel.csv)
# - pix_size is 3.2 for GALEX, and 0.45 for SDSS.

# This should allow you to work on the exact same pixel grid as I did for any given target.


# -----------------------------------------------------------------

class PoissonErrorCalculator(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Determine the path to a temporary directory
        self.temp_path = tempfile.gettempdir()

    # -----------------------------------------------------------------

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        pass

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

        return  get_galex_url_table()

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

    def get_header_for_galaxy(self, galaxy_name, instrument):

        """
        This function ...
        :param galaxy_name:
        :param instrument
        :return:
        """

        ra, dec, width = self.get_cutout_range_for_galaxy(galaxy_name)

        # Get the pixelscale
        pixelscale = dustpedia_final_pixelsizes[instrument]

        # Get the header
        return self.make_header(ra, dec, width, pixelscale)

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
        ra = self.leda_wise_table["ra2000"][index] * Unit("deg")
        dec = self.leda_wise_table["de2000"][index] * Unit("deg")

        # Get the D25
        d25_arcmin = self.leda_wise_table["d25"][index]

        if d25_arcmin < 6: width = .5 * Unit("deg")
        else: width = 1. * Unit("deg")

        # Return the RA, DEC and width (all in degrees)
        return ra, dec, width

    # -----------------------------------------------------------------

    def make_header(self, ra, dec, width, pix_size):

        """
        This function ...
        :return:
        """

        # ra and dec are taken from the ra2000 and de2000 of the attached DustPedia_LEDAWISE_Herschel.csv table
        # width is 0.5 degrees for galaxies with D25<6 arcmin, and 1 degree for galaxies with D25>=6 arcmin (as listed in DustPedia_LEDAWISE_Herschel.csv)
        # pix_size is 3.2 for GALEX, and 0.45 for SDSS.

        #

        ra = ra.to("deg").value
        dec = dec.to("deg").value
        width = width.to("deg").value
        pix_size = pix_size.to("arcsec").value

        #header_path = fs.join(self.temp_path, "header.hdr")

        header_path = fs.join(fs.home(), "header.hdr")

        # Create the header
        montage.commands.mHdr(str(ra) + ' ' + str(dec), width, header_path, pix_size=pix_size)

    # -----------------------------------------------------------------

    def get_sdss_dr9_field_for_sky_area(self, sky_area):

        """
        This function ...
        :param sky_area:
        :return:
        """

        montage.mArchiveList('2MASS', 'K', 'm31', 0.5, 0.5, 'm31.tbl')

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
