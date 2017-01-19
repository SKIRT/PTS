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
import shutil
import tempfile
import numpy as np
from scipy.interpolate import interpn, interp2d

# Import astronomical modules
import montage_wrapper as montage
from astropy.table import Table
from astropy.utils import lazyproperty
from astropy.units import Unit
from astropy.io.fits import Header, getheader
#from astroquery.sdss import SDSS
from astropy.io.fits import open as open_fits

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.tools import tables
from ...core.tools import network
from ...core.tools import time
from ...magic.core.image import Image
from .galex_montage_functions import mosaic_galex
from ...core.tools import archive
from ...magic.core.frame import Frame, sum_frames, sum_frames_quadratically
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.tools.general import split_xyz

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

dustpedia_final_pixelsizes = {"GALEX": 3.2 * Unit("arcsec"), "SDSS": 0.45 * Unit("arcsec"),
                              "2MASS": 1. * Unit("arcsec"), "WISE": 1.375 * Unit("arcsec")} # in arcsec, 2MASS and WISE pixel sizes are original pixelsizes

# -----------------------------------------------------------------

boss_photoobj_url =	"https://data.sdss.org/sas/dr13/eboss/photoObj"

# photoField	$BOSS_PHOTOOBJ/RERUN/RUN
# frame	        $BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL


# frame naming convention:

# frame-[ugriz]-[0-9]{6}-[1-6]-[0-9]{4}\.fits\.bz2, where [ugriz] is the filter, [0-9]{6} is a zero-padded six-digit number containing the run number, [1-6] is the camera column ('camcol') and [0-9]{4} is the zero-padded, four-digit frame sequence number.
# OR THUS:
# frame-[ugriz]-RUN-CAMCOL-FRAMESEQ.fits.bz2

# photoField naming convention:

# photoField-[0-9]{6}-[1-6]\.fits, where [0-9]{6} is the zero-padded, six digit run number and [1-6] is the camcol number.
# OR THUS:
# photoField-RUN-CAMCOL.fits  ## NO .BZ2 !!

# -----------------------------------------------------------------

sdss_filter_number = {"u": 0, "g": 1, "r": 2, "i": 3, "z": 4}

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

# Also, note that the DustPedia GALtEX cutouts have been re-gridded to larger pixel sizes than the standard GALEX
# archive data (the ancillary data report provides details).

# Also attached is a file that lists the URLs each GALEX observation can be downloaded from. Annoyingly, the URLs are
# *not* properly standardised; simply knowing the details of am observation is *not* enough to automatically construct
# the correct URL! I had to write a script that crawled through the GALEX GR6/7 archive website page-by-page,
# doing string-matching to find the appropriate URLs... I hate that archive.

# It's easy to work out which URLs correspond to which observations, as the URL for a given observation always
# contains the tilename (note that there will often be multiple observations, and hence URLS, associated with a given tilename).


#SDSS:

# For SDSS, you can use the Montage function mArchiveList to get a listing of the SDSS DR9 fields that cover a given
# area of sky (there is a nice Python wrapper available for Montage, which can make it easier to interact with).
# This function outputs a table giving the details of the relevant fields, including URLs where they can be downloaded.

# Note that I only used the SDSS primary fields to produce the DustPedia cutouts (the Ancillary Data Report I attached
# in the previous email explains this in more detail). I have attached a file that lists all the SDSS primary fields.
# After you have the results table produced by mArchiveGet, you can match the results with the primary fields list
# to find out the fields I used to make the final cutouts.

# Also, bear in mind that the DustPedia SDSS cutouts are re-gridded to North-East orientation, and 0.45" pixel sizes (again, see the ancillary data report for details).

# The SDSS database website provides extensive information about every SDSS field. It may be possible to use this
# database to find out the information you want about each field without having to download the actual FITS files,
# once you have the fields' ID information.


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


# The reason mArchiveList queries DR9 is that the SDSS imaging products haven't changed since DR9, so Montage haven't
# felt the need to explicitly update which data release they query, as they'd get the same results.
# Hence the lack of specifying the data release.



## NEW INFO:

# Firstly, checking my code, I realised that I didn't implement the "attachment to target location" criterion in the
# GALEX maps; it appears it basically never happens in GALEX. (I mainly brought in that criterion to handle Spitzer
# maps, where it is really common to have multiple unattached observations within a given 0.5x0.5 degree patch of sky.)

# To use Montage to find out what images cover a given part of the sky, you first want to run mImgTable on the folder
# containing your image files, like so:
# montage_wrapper.commands.mImgtbl('/folder/with/maps/', 'Image_Metadata_Table.dat', corners=True)

# The 'Image_Metadata_Table.dat' output file is a table recording the geometry of all the files in the folder in
# question. You then use this file as an input to the function mCoverageCheck, which is used like this:

# montage_wrapper.commands_extra.mCoverageCheck('Image_Metadata_Table.dat', 'Overlap_Table.dat', mode='box', ra=ra,
# dec=dec, width=width)

# The 'Overlap_Table.dat' contains the subset of the rows from 'Image_Metadata_Table.dat' that intersect the defined
# region (not just those that cover the centre coordinate). I then read in this in with:

# overlap_files = np.genfromtxt('Overlap_Table.dat', skip_header=3, usecols=[31], dtype=('S500'))

# Regarding doing the actual mosaicing, it was nice and straightforward with the SDSS data, where I used the
# high-level Montage command mExec:

# montage_wrapper.commands.mExec('SDSS', band, raw_dir='/path/to/files/', level_only=False, corners=False,
# debug_level=0, output_image='Mosaic.fits', region_header='Header.hdr', workspace_dir='/some/temp/dir')

# The nice thing about this command is that if you leave out the set raw_dir=None, then it will automatically
# retrieve the necessary SDSS primary fields from the SDSS server. However, when dealing with large numbers of
# cutouts in succession, I found it somewhat quicker to wget the files myself. The region_header is typically a
# file created by Montage's mHdr function.

# For GALEX, the mosaicing was a bit more involved:
#
# - I had to use the relative-response and sky-background maps to work out which part of each fits file was usable,
#   becuase raw GALEX tiles use 0 to indicate that a pixel recieved no photons... *or* that the pixel is outside
#   the detector footprint... MADNESS!
# - The co-addition had to be weighted by the exposure time. To do this, I created relative weight maps where each
#   pixel was just the square root of the 'exptime' header keyword value (where the relative-response and
#   sky-background maps were used to remove irrelevant pixels).
# - I had to adjust the background levels myself, as Montage hates all the zeros in GALEX maps.

# And whilst I did the re-projection with Montage, I did the actual coaddition using SWarp (as it has a far quicker
# runtime than doing it myself in Python). I have attached the very (very, very, very, very) ugly script I used to do
# all this. This should make it possible for you to work out the process I followed.

# Regarding not needing to make the full maps to work out the poisson noise, you're probably right.

# For the GALEX maps, you probably just need the input maps and the weight maps, re-projected to the output grid.

# For SDSS, all input maps are given equal weighting. By using mExec I skipped all of the 'intermediate' steps of
# the mosaicing. However you could just use mProjExec to regrid to all of the input fields for a given galaxy to the
# output projection mExec will have used. Jjust give the ra, dec, and width to mHdr to produce a header, and then
# use this header as input to mProjExec to regrid all the individual image to the final projection.


## GALEX exposure, background , etc.

# That's nice and straightforwad! You can use the same URLs as for the standard maps, but just replace
# '-int.fits' with '-rr.fits' for the relative response maps, or '-exp.fits' for exposure maps, etc.

# If you can look at the various files available for each tilenum by looking at this page for each tile:

# 'http://galex.stsci.edu/GR6/?page=downloadlist&tilenum='+str(tilenum)+'&type=coaddI&product=Imaging%20Only'

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

    def download_galex_observations_for_galaxy(self, galaxy_name, images_path, response_path, background_path, counts_path):

        """
        This function ...
        :param galaxy_name:
        :param images_path:
        :param response_path:
        :param background_path:
        :param counts_path:
        :return:
        """

        # Inform the user
        log.info("Downloading the GALEX observations for " + galaxy_name + " to '" + images_path + "' ...")

        # Get the urls
        urls = self.get_galex_observation_urls_for_galaxy(galaxy_name)

        # Debugging
        log.debug("Number of observations that will be downloaded: " + str(len(urls)))

        # Download the files
        paths = network.download_files(urls, images_path)

        # Debugging
        log.debug("Decompressing the files ...")

        # Decompress the files and remove the originals
        archive.decompress_files(paths, remove=True)

        ## RESPONSE

        # Inform the user
        log.info("Downloading the GALEX response maps for " + galaxy_name + " to '" + response_path + "' ...")

        # Determine the URLS for the response maps
        response_urls = [url.replace("-int.fits.gz", "-rr.fits.gz") for url in urls]

        # Download the response maps
        response_paths = network.download_files(response_urls, response_path)

        # Debugging
        log.debug("Decompressing the response files ...")

        # Decompress
        archive.decompress_files(response_paths, remove=True)

        ## BACKGROUND

        # Inform the user
        log.info("Downloading the GALEX background maps for " + galaxy_name + " to '" + background_path + "' ...")

        # Determine the URLs for the background maps
        background_urls = [url.replace("-int.fits.gz", "-skybg.fits.gz") for url in urls]

        # Download the background maps
        background_paths = network.download_files(background_urls, background_path)

        # Debugging
        log.debug("Decompressing the background files ...")

        # Decompress
        archive.decompress_files(background_paths, remove=True)

        ## COUNTS

        # Inform the user
        log.info("Downloading the GALEX maps in counts for " + galaxy_name + " to '" + counts_path + "' ...")

        # Determine the URLS for the maps in original detector counts
        counts_urls = [url.replace("-int.fits.gz", "-cnt.fits.gz") for url in urls]

        # Download the count maps
        counts_paths = network.download_files(counts_urls, counts_path)

        # Debugging
        log.debug("Decompressing the count maps ...")

        # Decompress
        archive.decompress_files(counts_paths, remove=True)

    # -----------------------------------------------------------------

    def get_galex_observation_urls_for_galaxy(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Get the tilenames
        tilenames = self.get_galex_tilenames_for_galaxy(galaxy_name)

        # Search through the URL table to get all the URLS that contain one of the tilenames
        urls = []
        for i in range(len(self.galex_url_table)):

            url = self.galex_url_table["URL"][i]

            for tilename in tilenames:

                if tilename in url:
                    urls.append(url)
                    break # URL added, so no need to look at the other tilenames for this URL

        # Return the list of URLS
        return urls

    # -----------------------------------------------------------------

    def get_galex_tilenames_for_galaxy(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Find the indices in the table
        indices = tables.find_indices(self.galex_observations_table, galaxy_name, column_name="uploadID")

        tilenames = set()

        # Loop over the indices corresponding to the specified galaxy
        for index in indices:

            tilename = self.galex_observations_table["tilename"][index]
            tilenames.add(tilename)

        # Return the tilenames as a list
        return list(tilenames)

    # -----------------------------------------------------------------

    def make_galex_mosaic_and_poisson_frame(self, galaxy_name, output_path):

        """
        This function ...
        :param galaxy_name:
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Making GALEX mosaic for " + galaxy_name + " and map of relative poisson errors ...")

        # Create directories
        working_path, download_path, response_path, background_path, counts_path, raw_path, temp_path = self.create_base_directories(galaxy_name)

        # Create temp subdirectories
        temp_fuv_path = fs.join(temp_path, "FUV")
        temp_nuv_path = fs.join(temp_path, "NUV")
        fs.create_directories(temp_fuv_path, temp_nuv_path)

        # Create download subdirectories
        download_images_path, download_response_path, download_background_path, download_counts_path = self.create_download_directories(download_path)

        # Download data and split into FUV and NUV
        self.download_and_split(galaxy_name, download_path, download_images_path, download_response_path, download_background_path, download_counts_path, raw_path, response_path, background_path, counts_path)

        raw_fuv_path = fs.join(raw_path, "FUV")
        raw_nuv_path = fs.join(raw_path, "NUV")

        # Inform the user
        log.info("Getting cutout range ...")

        # Get coordinate range for target image
        ra, dec, width = self.get_cutout_range_for_galaxy(galaxy_name)

        # Check coverage
        bands_dict, bands_in_dict, metadata_paths = self.check_galex_coverage(galaxy_name, ra, dec, working_path, temp_path, raw_fuv_path, raw_nuv_path)

        # Loop over the two bands, perform the mosaicing
        for band in bands_in_dict.keys():

            # Make the mosaic for this and
            mosaic_galex(galaxy_name, ra, dec, width, bands_dict[band], working_path, temp_path, metadata_paths[band], output_path)

    # -----------------------------------------------------------------

    def create_base_directories(self, galaxy_name, create=True):

        """
        This function ...
        :param galaxy_name:
        :param create:
        :return:
        """

        # Inform the user
        log.info("Creating directories ...")

        # Determine the path to the temporary directory for downloading the images
        working_path = fs.join(fs.home(), time.unique_name("GALEX_" + galaxy_name))

        #working_path = fs.join(fs.home(), "GALEX_NGC3031_2016-08-04--10-59-27-047")

        # Create the working directory
        if create: fs.create_directory(working_path)

        # DOWNLOAD PATH
        download_path = fs.join(working_path, "download")

        # Create download directory
        if create: fs.create_directory(download_path)

        # RESPONSE AND BACKGROUND PATH
        response_path = fs.join(working_path, "response")
        background_path = fs.join(working_path, "background")
        if create: fs.create_directories(response_path, background_path)

        # COUNT PATH
        counts_path = fs.join(working_path, "counts")
        if create: fs.create_directory(counts_path)

        # RAW PATH
        raw_path = fs.join(working_path, "raw")
        # Create raw directory
        if create: fs.create_directory(raw_path)

        # TEMP PATH
        temp_path = fs.join(working_path, "temp")

        # Create temp directory
        if create: fs.create_directory(temp_path)

        # Return ...
        return working_path, download_path, response_path, background_path, counts_path, raw_path, temp_path

    # -----------------------------------------------------------------

    def create_download_directories(self, download_path, create=True):

        """
        This function ...
        :param download_path:
        :param create:
        :return:
        """

        # download/images, download/response and download/background
        download_images_path = fs.join(download_path, "images")
        download_response_path = fs.join(download_path, "reponse")
        download_background_path = fs.join(download_path, "background")
        download_counts_path = fs.join(download_path, "counts")

        # Create
        if create: fs.create_directories(download_images_path, download_response_path, download_background_path, download_counts_path)

        # Return
        return download_images_path, download_response_path, download_background_path, download_counts_path

    # -----------------------------------------------------------------

    def download_and_split(self, galaxy_name, download_path, download_images_path, download_response_path, download_background_path,
                           download_counts_path, raw_path, response_path, background_path, counts_path, create=True):

        """
        This function ...
        :return:
        """

        ###

        # 1 and 2 RAW directories
        raw_fuv_path = fs.join(raw_path, "FUV")
        raw_nuv_path = fs.join(raw_path, "NUV")
        if create: fs.create_directories(raw_fuv_path, raw_nuv_path)

        # FUV and NUV response directories
        response_fuv_path = fs.join(response_path, "FUV")
        response_nuv_path = fs.join(response_path, "NUV")
        if create: fs.create_directories(response_fuv_path, response_nuv_path)

        # FUV and NUV background directories
        background_fuv_path = fs.join(background_path, "FUV")
        background_nuv_path = fs.join(background_path, "NUV")
        if create: fs.create_directories(background_fuv_path, background_nuv_path)

        # FUV AND NUV counts directories
        counts_fuv_path = fs.join(counts_path, "FUV")
        counts_nuv_path = fs.join(counts_path, "NUV")
        if create: fs.create_directories(counts_fuv_path, counts_nuv_path)

        ###

        # Download the GALEX observations to the temporary directory  # they are decompressed here also
        if create: self.download_galex_observations_for_galaxy(galaxy_name, download_images_path, download_response_path,
                                                    download_background_path, download_counts_path)

        # Inform the user
        log.info("Splitting observations into FUV and NUV ...")

        # Split downloaded images into FUV and NUV
        if create: self.split_galex_observations(download_images_path, raw_fuv_path, raw_nuv_path)

        # Split response maps into FUV and NUV
        if create: self.split_galex_observations(download_response_path, response_fuv_path, response_nuv_path)

        # Split background maps into FUV and NUV
        if create: self.split_galex_observations(download_background_path, background_fuv_path, background_nuv_path)

        # Split count maps into FUV and NUV
        if create: self.split_galex_observations(download_counts_path, counts_fuv_path, counts_nuv_path)

        ###

        # AFTER SPLIT, THIS CAN BE DONE:

        fs.remove_directory(download_images_path)
        fs.remove_directory(download_response_path)
        fs.remove_directory(download_background_path)
        fs.remove_directory(download_counts_path)
        fs.remove_directory(download_path)

    # -----------------------------------------------------------------

    def check_galex_coverage(self, galaxy_name, ra, dec, working_path, temp_path, raw_fuv_path, raw_nuv_path):

        """
        This function ...
        :param galaxy_name:
        :param ra:
        :param dec:
        :param working_path:
        :param temp_path:
        :param raw_fuv_path:
        :param raw_nuv_path:
        :return:
        """

        ## Generate the meta and then overlap file
        ##meta_path, overlap_path = self.generate_meta_and_overlap_file(temp_path, ra, dec, width)

        ##meta_path = fs.join(temp_path, "meta.dat")
        ## Get the image table of which images cover a given part of the sky
        ##montage.commands.mImgtbl(temp_path, meta_path, corners=True)

        raw_band_paths = {"FUV": raw_fuv_path, "NUV": raw_nuv_path}

        # State band information
        bands_dict = {'FUV': {'band_short': 'fd', 'band_long': 'FUV'},
                      'NUV': {'band_short': 'nd', 'band_long': 'NUV'}}

        metadata_paths = {"FUV": fs.join(working_path, "FUV_Image_metadata_table.dat"),
                          "NUV": fs.join(working_path, "NUV_Image_metadata_table.dat")}

        # If not yet done, produce Montage image table of raw tiles
        for band in bands_dict.keys():

            # Inform the user
            log.info("Creating " + band + " image metadata table ...")

            metadata_path = metadata_paths[band]
            montage.commands.mImgtbl(raw_band_paths[band], metadata_path, corners=True)

        # Check if any GALEX tiles have coverage of source in question; if not, continue
        bands_in_dict = {}
        for band in bands_dict.keys():

            # Inform the user
            log.info("Creating " + band + " overlap file ...")

            # OVERLAP TABLE IS HERE TEMPORARY
            overlap_path = fs.join(temp_path, "overlap_" + band + ".dat")

            # Create overlap file
            montage.commands_extra.mCoverageCheck(metadata_paths[band], overlap_path, mode='point', ra=ra.to("deg").value, dec=dec.to("deg").value)

            # Check if there is any coverage for this galaxy and band
            if sum(1 for line in open(overlap_path)) <= 3:
                log.warning("No GALEX " + band + " coverage for " + galaxy_name)
            else: bands_in_dict[band] = bands_dict[band]

            # Remove overlap file
            fs.remove_file(overlap_path)

        # Check if coverage in any band
        if len(bands_in_dict) == 0: raise RuntimeError("No coverage in any GALEX band")

        # Return ...
        return bands_dict, bands_in_dict, metadata_paths

    # -----------------------------------------------------------------

    def split_galex_observations(self, original_path, path_fuv, path_nuv):

        """
        This function ...
        :param original_path:
        :param path_fuv:
        :param path_nuv:
        :return:
        """

        # Inform the user
        log.info("Splitting GALEX observations from '" + original_path + "' into '" + path_fuv + "' and '" + path_nuv + "' ...")

        # Loop over the files in the path
        for path, name in fs.files_in_path(original_path, extension="fits", returns=["path", "name"]):

            # Get header
            header = getheader(path)

            # Check band
            band = header["BAND"]

            # PREVIOUSLY
            #if band == 1: shutil.copy(path, path_fuv)
            #elif band == 2: shutil.copy(path, path_nuv)
            #else: raise RuntimeError("Invalid band: " + str(band))

            # NEW
            if band == 1: shutil.copy(path, path_nuv)
            elif band == 2: shutil.copy(path, path_fuv)
            else: raise RuntimeError("Invalid band: " + str(band))

    # -----------------------------------------------------------------

    def make_sdss_mosaic_and_poisson_frame(self, galaxy_name, band, output_path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param output_path:
        :return:
        """

        # FROM DR7 page: (http://classic.sdss.org/dr7/algorithms/fluxcal.html#counts2mag)

        # Computing errors on counts (converting counts to photo-electrons)
        # The fpC (corrected frames) and fpObjc (object tables with counts for each object instead of magnitudes)
        # files report counts (or "data numbers", DN). However, it is the number of photo-electrons which is really
        # counted by the CCD detectors and which therefore obeys Poisson statistics. The number of photo-electrons
        # is related to the number of counts through the gain (which is really an inverse gain):
        #
        #      photo-electrons = counts * gain
        #
        # The gain is reported in the headers of the tsField and fpAtlas files (and hence also in the field table in
        # the CAS). The total noise contributed by dark current and read noise (in units of DN2) is also reported in
        # the tsField files in header keyword dark_variance (and correspondingly as darkVariance in the field table in
        # the CAS), and also as dark_var in the fpAtlas header.
        #
        # Thus, the error in DN is given by the following expression:
        #
        #      error(counts) = sqrt([counts+sky]/gain + Npix*(dark_variance+skyErr)),
        #
        # where counts is the number of object counts, sky is the number of sky counts summed over the same area as
        # the object counts, Npix is the area covered by the object in pixels, gain and dark_variance and skyErr are
        # the gain, dark variance, and the error on the estimate of the average sky level in the frame, respectively,
        # from the corresponding tsField file.

        ####
        # BUT: Note: fpC files are no longer produced as of DR8. Instead, equivalent data is stored in frame files.
        # (https://data.sdss.org/datamodel/files/PHOTO_REDUX/RERUN/RUN/objcs/CAMCOL/fpC.html)
        ###

        # EXPLANATION OF (POISSON) NOISE CALCULATION FROM DR8 ONWARDS !!!:

        # on DATA MODEL page for 'frame': https://data.sdss.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
        #
        # Frame = The calibrated, sky-subtracted corrected frame plus associated calibration meta-data.
        # The units of the images are in nanomaggies. It is kept compressed under "bzip2", which we have found is
        # the most efficient compressor of this format. In addition, there is a lossy compression applied to the
        # floating point values (which retains accuracy at the 0.1 percent level). The IDL routine "read_frame.pro"
        # in photoop will back out the calibration and sky-subtraction from this file if necessary, in steps explained
        #  below. Also explained below is how to calculate the noise in the image.
        #

        # HDU0: the corrected frame, what is normally in the "fpC" files
        #
        #  The "image", a 2048x1489 array of floating point values, the calibrated and sky-subtracted version
        #  of the fpC "corrected frame" files produced by photo. Units are in nanomaggies.
        #
        # HDU1: the flat-field and calibration vector
        #
        #  The "calibvec", a 2048-element array of "float32" values, encompassing the flat-field correction to
        #  apply, multiplied by the calibration. Translates the counts in the original image into nanomaggies.
        #
        # HDU2: the sky image
        #
        #

        # ;; 0. find filename of the frame file
        # framename = (sdss_name('frame', run, camcol, field, $
        #               filter=filternum(filter), rerun=rerun))[0]+'.bz2'
        #
        # ;; 1. read in the FITS image from HDU0; the resulting image will be
        # ;;    sky-subtracted as well as calibrated in nanomaggies/pixel
        # img = mrdfits(framename,0,hdr)
        # nrowc= (size(img,/dim))[1]
        #
        # ;; 2. read in sky, and interpolate to full image size; this returns a
        #;;    sky image the same size as the frame image, in units of counts
        # sky = mrdfits(framename,2)
        # simg = interpolate(sky.allsky, sky.xinterp, sky.yinterp, /grid)
        #
        # ;; 3. read in calibration, and expand to full image size; this returns
        # ;;    a calibration image the same size as the frame image, in units of
        # ;;    nanomaggies per count
        # calib= mrdfits(framename,1)
        # cimg= calib#replicate(1.,nrowc)
        #

        # DATA MODEL DIRECTORIES: https://data.sdss.org/datamodel/files/

        #   BOSS_PHOTOOBJ  https://data.sdss.org/sas/dr13/eboss/photoObj
        #   PHOTO_CALIB    https://data.sdss.org/sas/dr13/eboss/calib/dr13_final
        #   PHOTO_SKY      https://data.sdss.org/sas/dr13/eboss/photo/sky

        # DATA MODEL: frame : https://data.sdss.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html  (DR13?)

        # DATA MODEL: photoField table (GAIN AND DARKVARIANCE) : https://data.sdss.org/datamodel/files/BOSS_PHOTOOBJ/RERUN/RUN/photoField.html  (DR13?)

        # SKY FRAMES: https://data.sdss.org/datamodel/files/PHOTO_SKY/RERUN/RUN/sky/skyframes.html

        #

        # MAKING SQL QUERY: http://skyserver.sdss.org/dr12/en/help/docs/realquery.aspx

        # SDSS SOFTWARE UTILS: http://www.sdss.org/dr13/software/products/

        # PHOTOOP IDL PACKAGE LATEST VERSION: https://svn.sdss.org/public/repo/sdss/photoop/tags/v1_9_9/
        # svn export https://svn.sdss.org/public/repo/sdss/photoop/tags/v1_9_9/ photoop

        # IDLUTILS IDL PACKAGE LATEST VERSION: https://svn.sdss.org/public/repo/sdss/idlutils/tags/v5_5_9/
        # svn export https://svn.sdss.org/public/repo/sdss/idlutils/tags/v5_5_9/ idlutils

        ##
        # PyIDL installation: https://pypi.python.org/pypi/pyIDL/

        ###

        # MAKE PATHS

        # Determine the path to the temporary directory
        temp_path = fs.join(fs.home(), time.unique_name("SDSS_" + galaxy_name + "_" + band))

        # Create the temporary directory
        fs.create_directory(temp_path)

        #temp_path = fs.join(fs.home(), "SDSS_NGC3031_u_2016-08-10--16-48-33-775")

        ###

        # RAW PATH
        raw_path = fs.join(temp_path, "raw") # frames with HDU0, HDU1, HDU2 .. (primary, calib, sky, ...) IN NANOMAGGIES PER PIXEL
        fs.create_directory(raw_path)

        # FIELDS PATH
        fields_path = fs.join(temp_path, "fields") # photoField files
        fs.create_directory(fields_path)

        # COUNTS PATH
        counts_path = fs.join(temp_path, "counts") # frames IN COUNTS (DN)
        fs.create_directory(counts_path)

        # POISSON PATH
        poisson_path = fs.join(temp_path, "poisson_nmgy")  # error maps in NANOMAGGIES PER PIXEL
        fs.create_directory(poisson_path)

        # REBINNED PATH
        rebinned_path = fs.join(temp_path, "rebinned") # images with frame and corresponding error map in NANOMAGGIES, REBINNED
        fs.create_directory(rebinned_path)

        # FOOTPRINTS PATH
        footprints_path = fs.join(temp_path, "footprints")
        fs.create_directory(footprints_path)

        # MOSAIC PATH
        mosaics_path = fs.join(temp_path, "mosaics")
        fs.create_directory(mosaics_path)

        # RESULT PATH
        results_path = fs.join(temp_path, "results")
        fs.create_directory(results_path)

        ###

        # GET AND DECALIBRATE FRAMES

        # Inform the user
        log.info("Getting SDSS frames in nanomaggies (per pixel) and converting to counts (DN) ...")

        # Make SDSS frames in counts
        self.make_sdss_frames_in_counts(galaxy_name, band, raw_path, fields_path, counts_path)

        # Make (Poisson) noise maps
        self.make_sdss_noise_maps_in_nanomaggies(raw_path, counts_path, fields_path, poisson_path)

        ## GET THE TARGET HEADER

        # Get the target header
        header = self.get_header_for_galaxy(galaxy_name, "SDSS")

        # To coordinate system
        rebin_wcs = CoordinateSystem(header)

        ##

        # DO REBINNING, CREATE IMAGES WITH REBINNED PRIMARY AND ERROR FRAME IN NANOMAGGIES PER PIXEL, AND FOOTPRINT FILES
        self.rebin_sdss_frames_and_error_maps(rebin_wcs, raw_path, poisson_path, rebinned_path, footprints_path)

        ## Make the footprints: WAS ONLY NECESSARY FOR WHEN FOOTPRINTS WERE NOT CREATED DURING FUNCTION ABOVE
        ##self.make_sdss_footprints(rebinned_path, footprints_path)

        # DO THE COMBINING
        # rebinned_path, footprints_path, mosaics_path, wcs
        self.combine_sdss_frames_and_error_maps(rebinned_path, footprints_path, mosaics_path, rebin_wcs)

        # CONVERT TO JY/PIX
        self.convert_sdss_mosaic_and_error_map_to_jansky(galaxy_name, band, mosaics_path, results_path)

        ## WRITE RESULT TO OUTPUT DIRECTORY

        # Load the image
        id_string = galaxy_name + '_SDSS_' + band
        result_path = fs.join(results_path, id_string + ".fits")
        image = Image.from_file(result_path)

        # Get mosaic and error map
        mosaic = image.frames.primary
        mosaic_errors = image.frames.errors

        # Calculate the relative error map
        relerrors = mosaic_errors / mosaic
        relerrors[relerrors < 0.] = 0.0 # set negative values for relative error map to zero
        relerrors.replace_nans(0.0)     # set NaN values (because mosaic was zero) to zero

        # Save mosaic as FITS file
        mosaic_output_path = fs.join(output_path, id_string + ".fits")
        mosaic.saveto(mosaic_output_path)

        # Save error map as FITS file
        errors_output_path = fs.join(output_path, id_string + "_errors.fits")
        mosaic_errors.saveto(errors_output_path)

        # Save relative error map as FITS file
        relerrors_output_path = fs.join(output_path, id_string + "_relerrors.fits")
        relerrors.saveto(relerrors_output_path)

        ## END

    # -----------------------------------------------------------------

    def make_sdss_footprints(self, rebinned_path, footprints_path):
        
        """
        This function was used for when the rebinning was already performed but we forgot to include saving the footprints 
        back then. Luckily, the footprints are equivalent to the masks of the pixel that are not set to NaN in the rebinned frames
        """
        
        # Loop over the files in the rebinned path
        for path, name in fs.files_in_path(rebinned_path, extension="fits", returns=["path", "name"]):
            
            # Open the rebinned frame
            frame = Frame.from_file(path)
            
            # Create the footprint
            footprint = Frame(frame.nans().inverse().astype(float))
            footprint.wcs = frame.wcs
        
            # Save the footprint
            footprint_path = fs.join(footprints_path, name + ".fits")
            footprint.saveto(footprint_path)
        
    # -----------------------------------------------------------------

    def make_sdss_noise_maps_in_nanomaggies(self, raw_path, counts_path, fields_path, poisson_path):

        """
        This function ...
        :param raw_path:
        :param counts_path:
        :param fields_path:
        :param poisson_path:
        :return:
        """

        # Inform the user
        log.info("Computing error maps for each SDSS frame in nanomaggies per pixel ...")

        # Loop over the files in the counts path
        for path, name in fs.files_in_path(counts_path, extension="fits", returns=["path", "name"]):

            # EXAMPLE FILE NAME: frame-u-004294-4-0231.fits

            # FIELD URL = field_url_start / RERUN / RUN    + / photoField-6digits-CAMCOL.fits
            # example: http://data.sdss3.org/sas/dr12/env/BOSS_PHOTOOBJ/301/4294/photoField-004294-5.fits

            # FRAME URL = $BOSS_PHOTOOBJ / frames / RERUN / RUN / CAMCOL    +   /frame-[ugriz]-6digits-CAMCOL-FRAMESEQ.fits.bz2
            # example: http://data.sdss3.org/sas/dr10/boss/photoObj/frames/301/4294/5/frame-i-004294-5-0229.fits.bz2

            splitted = name.split("-")

            band = splitted[1]
            digits = splitted[2]
            camcol = splitted[3]
            frameseq = splitted[4]

            # Determine the path to the corresponding field file
            field_file_path = fs.join(fields_path, "photoField-" + digits + "-" + camcol + ".fits")

            # Get the gain and dark variance
            gain, dark_variance = self.get_gain_and_dark_variance_from_photofield(field_file_path, band)

            # Calculate the error

            # dn_err = sqrt(dn/gain+darkVariance)   # NOISE IN DN

            # img_err = dn_err*cimg                 # NOISE IN NANOMAGGIES / PIXEL


            # Load the image in DN
            dn_frame = Frame.from_file(path)

            # COMPUTE THE DN ERROR
            dn_error = np.sqrt(dn_frame.data / gain + dark_variance)

            # CONVERT ERROR MAP IN DN TO ERROR MAP IN NANOMAGGIES PER PIXEL

            # COMPUTE THE CALIBRATION FRAME AGAIN
            raw_frame_path = fs.join(raw_path, name + ".fits")
            hdulist = open_fits(raw_frame_path)
            img = hdulist[0].data
            nrowc = img.shape[0]  # ysize    (x = columns, y = rows)
            # Get the calibration HDU (header)
            calib = hdulist[1].data
            replicate = np.ones(nrowc)
            cimg = np.outer(replicate, calib)

            # IMAGE ERROR MAP IN NANOMAGGIES PER PIXEL
            img_error = dn_error * cimg
            error_map = Frame(img_error)

            # SET THE WCS OF THE IMAGE ERROR MAP
            error_map.wcs = dn_frame.wcs

            # SAVE the error map
            poisson_error_map_path = fs.join(poisson_path, name + ".fits")
            error_map.saveto(poisson_error_map_path)

    # -----------------------------------------------------------------

    def rebin_sdss_frames_and_error_maps(self, rebin_wcs, raw_path, poisson_path, rebinned_path, footprints_path):

        """
        This function ...
        :param rebin_wcs:
        :param raw_path:
        :param poisson_path:
        :param rebinned_path:
        :param footprints_path
        :return:
        """

        # Inform the user
        log.info("Rebinning primary SDSS frames and corresponding error maps to the target coordinate system ...")

        # Loop over the files in the raw directory
        for path, name in fs.files_in_path(raw_path, extension="fits", returns=["path", "name"]):

            # Open the frame IN NANOMAGGIES PER PIXEL
            frame = Frame.from_file(path)  # HDU 0 is used

            # Open the corresponding error map IN NANOMAGGIES PER PIXEL
            error_path = fs.join(poisson_path, name + ".fits")
            errormap = Frame.from_file(error_path)

            # Debugging
            log.debug("Converting " + name + " frame and error map to nanomaggy / sr...")

            # NUMBER OF SR PER PIXEL
            pixelsr = frame.pixelarea.to("sr").value

            # CONVERT FRAME
            frame /= pixelsr  # IN NANOMAGGIES PER SR NOW

            # CONVERT ERROR MAP
            errormap /= pixelsr  # IN NANOMAGGIES PER SR NOW

            # Debugging
            log.debug("Rebinning " + name + " frame and error map to the target coordinate system ...")

            # Rebin the frame
            footprint = frame.rebin(rebin_wcs, exact=False)

            # REBIN THE ERROR MAP
            errormap.rebin(rebin_wcs, exact=False)

            # Debugging
            log.debug("Converting " + name + " frame and error map back to nanomaggy ...")

            # NUMBER OF SR PER PIXEL
            new_pixelsr = frame.pixelarea.to("sr").value

            # CONVERT THE REBINNED FRAME BACK TO NANOMAGGY (PER PIXEL)
            frame *= new_pixelsr

            # CONVERT THE REBINNED ERROR MAP BACK TO NANOMAGGY
            errormap *= new_pixelsr

            # CREATE IMAGE
            image = Image()

            image.add_frame(frame, "primary")
            image.add_frame(errormap, "noise")

            # Determine path for the image
            image_path = fs.join(rebinned_path, name + ".fits")
            image.saveto(image_path)

            # SAVE THE FOOTPRINT
            footprint_path = fs.join(footprints_path, name + ".fits")
            footprint.saveto(footprint_path)

    # -----------------------------------------------------------------

    def combine_sdss_frames_and_error_maps(self, rebinned_path, footprints_path, mosaics_path, wcs):

        """
        This function ...
        :param rebinned_path:
        :param footprints_path:
        :param mosaics_path:
        :param wcs:
        :return:
        """

        # Inform the user
        log.info("Combining the frames and error maps of the different SDSS fields ...")

        # MOSAIC[x,y](nanoMaggy) = SUM_i^n_overlapping_frames[x,y] ( d[x,y]_i ) / n_overlapping_frames[x,y]

        # where n_overlapping_frames[x,y] = SUM_i^m ( footprint[x,y]_i )

        # where m = the total number of fields

        # BUT BECAUSE d[x,y]_i = 0 for field that doesn't overlap in pixel [x,y], we can replace the sum from i=0 to n_overlapping_frames[x,y]
        # to a sum over ALL FIELDS; the contribution of the other fields will just be zero

        primary_frames = []
        error_frames = []
        n_overlapping_frames = []

        # Loop over the images in the 'rebinned' directory
        for path, name in fs.files_in_path(rebinned_path, extension="fits", returns=["path", "name"]):

            # Determine the path to the footprint
            footprint_path = fs.join(footprints_path, name + ".fits")

            # Open the footprint
            footprint = Frame.from_file(footprint_path)

            # Open the image
            image = Image.from_file(path)

            # Get the a and b frame
            a = image.frames.primary    # IN NANOMAGGIES PER PIXEL
            b = image.frames.noise      # IN NANOMAGGIES PER PIXEL

            # SET NANS TO ZERO
            a.replace_nans(0.0)
            b.replace_nans(0.0)

            # Add to the list
            primary_frames.append(a)
            error_frames.append(b)
            n_overlapping_frames.append(footprint)

        # Calculate n_overlapping_frames array
        n_overlapping = sum_frames(*n_overlapping_frames) # = n_overlapping_frames[x,y] = 2D ARRAY !

        # CALCULATE THE MOSAIC FRAME IN NANOMAGGIES
        mosaic_frame = sum_frames(*primary_frames) / n_overlapping
        mosaic_frame.wcs = wcs

        # CALCULATE THE MOSAIC ERROR MAP IN NANOMAGGIES
        mosaic_errormap = sum_frames_quadratically(*error_frames) / n_overlapping
        mosaic_errormap.wcs = wcs

        # SAVE THE MOSAIC IN NANOMAGGIES
        mosaic_path = fs.join(mosaics_path, "mosaic.fits")
        mosaic_frame.saveto(mosaic_path)

        # SAVE THE MOSAIC ERROR MAP IN NANOMAGGIES
        mosaic_error_path = fs.join(mosaics_path, "mosaic_errors.fits")
        mosaic_errormap.saveto(mosaic_error_path)

    # -----------------------------------------------------------------

    def convert_sdss_mosaic_and_error_map_to_jansky(self, galaxy_name, band, mosaics_path, results_path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param mosaics_path:
        :param results_path:
        :return:
        """

        # DETERMINE ID STRING TO SAVE THE RESULT
        id_string = galaxy_name + '_SDSS_' + band

        # Inform the user
        log.info("Converting the SDSS mosaic and error map to Jansky per pixel ...")

        # 1 nanomaggy = approximately 3.613e-6 Jy

        # Open the mosaic
        mosaic_path = fs.join(mosaics_path, "mosaic.fits")
        mosaic = Frame.from_file(mosaic_path)

        # Open the mosaic error map
        mosaic_error_path = fs.join(mosaics_path, "mosaic_errors.fits")
        mosaic_errors = Frame.from_file(mosaic_error_path)

        # DO THE CONVERSION FOR THE MOSAIC, SET NEW UNIT
        mosaic *= 3.613e-6
        mosaic.unit = "Jy / pix"

        # DO THE CONVERSION FOR THE MOSAIC ERROR MAP, SET NEW UNIT
        mosaic_errors *= 3.613e-6
        mosaic_errors.unit = "Jy / pix"


        # CREATE IMAGE
        image = Image()
        image.add_frame(mosaic, "primary")
        image.add_frame(mosaic_errors, "errors")

        # SAVE THE MOSAIC IMAGE (with error map) IN JANSKY
        result_path = fs.join(results_path, id_string + ".fits")
        image.saveto(result_path)

    # -----------------------------------------------------------------

    def sdss_mosaic_old_function_things(self, galaxy_name, band, output_path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Making SDSS mosaic and map of relative poisson errors ...")

        # Make rebinned frames in counts (and footprints)
        self.make_sdss_rebinned_frames_in_counts(galaxy_name, band, output_path)

        # Initialize a list to contain the frames to be summed
        ab_frames = []
        b_frames = []

        # Loop over the files in the output directory
        for path, name in fs.files_in_path(output_path, extension="fits", returns=["path", "name"]):

            # Open the image
            image = Image.from_file(path)

            # Get the a and b frame
            a = image.frames.primary    # in counts
            b = image.frames.footprint  # between 0 and 1

            # Add product of primary and footprint and footprint to the appropriate list
            ab = a * b
            ab_frames.append(ab)
            b_frames.append(b)

        # Take the sums
        ab_sum = sum_frames(*ab_frames)
        b_sum = sum_frames(*b_frames)

        # Calculate the relative poisson errors
        rel_poisson_frame = ab_sum**(-0.5)

        # Calculate the mosaic frame
        mosaic_frame = ab_sum / b_sum

        # Make image
        image = Image()

        # Add frames
        image.add_frame(mosaic_frame, "primary")
        image.add_frame(rel_poisson_frame, "rel_poisson") # has no unit, but Image will be saved with unit. Problem?

        # Save the image
        filename = galaxy_name + "_SDSS_" + band + ".fits"
        path = fs.join(output_path, filename)
        image.saveto(path)

    # -----------------------------------------------------------------

    def make_sdss_frames_in_counts(self, galaxy_name, band, raw_path, fields_path, counts_path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param raw_path:
        :param fields_path:
        :param counts_path:
        :return:
        """

        # Inform the user
        log.info("Making SDSS frames in counts for " + galaxy_name + " for SDSS " + band + " band ...")

        # ----

        # Download the FITS files to be used for mosaicing
        urls = self.download_sdss_primary_fields_for_galaxy_for_mosaic(galaxy_name, band, raw_path)

        #### NEW: DOWNLOAD THE FIELD TABLES

        #dr10 / boss / photoObj / frames

        field_url_start = "http://data.sdss3.org/sas/dr12/env/BOSS_PHOTOOBJ" # = $BOSS_PHOTOOBJ/RERUN/RUN
        # NOT WORKING WITH DR13 : BECAUSE OF PERMISSION ISSUES (PASSWORD IS ASKED WHEN ENTERED IN BROWSER !!)

        # FIELD URL = field_url_start / RERUN / RUN    + / photoField-6digits-CAMCOL.fits
        # example: http://data.sdss3.org/sas/dr12/env/BOSS_PHOTOOBJ/301/4294/photoField-004294-5.fits

        # FRAME URL = $BOSS_PHOTOOBJ / frames / RERUN / RUN / CAMCOL    +   /frame-[ugriz]-6digits-CAMCOL-FRAMESEQ.fits.bz2
        # example: http://data.sdss3.org/sas/dr10/boss/photoObj/frames/301/4294/5/frame-i-004294-5-0229.fits.bz2

        field_urls = []
        for url in urls:

            rerun_run_camcol = url.split("frames/")[1].split("/frame-")[0]
            rerun, run, camcol = rerun_run_camcol.split("/")
            rerun_run = rerun + "/" + run

            band, digits, camcol, frameseq = url.split("frame-")[1].split(".fits")[0].split("-")

            field_url = field_url_start + "/" + rerun_run + "/photoField-" + digits + "-" + camcol + ".fits"
            field_urls.append(field_url)

        # Remove duplicates (one field table contains the info for all 5 SDSS bands!)
        field_urls = list(set(field_urls))

        # DOWNLOAD THE FIELD FITS FILES

        # Debugging
        log.debug("Downloading the photoField data files ...")

        # Download the files
        field_paths = network.download_files(field_urls, fields_path)


        ## CONVERT THE IMAGES TO COUNT (DECALIBRATE)

        log.warning("The frames now have to be converted to counts (DN) manually!")

        log.warning("Execute the following IDL commands to convert a frame:")

        # 1. read in the FITS image from HDU0; the resulting image will be
        # sky-subtracted as well as calibrated in nanomaggies/pixel
        log.info("IDL> img= mrdfits(framename,0,hdr)")
        log.info("IDL> nrowc= (size(img,/dim))[1]")
        log.info("")

        # 2. read in sky, and interpolate to full image size; this returns a
        # sky image the same size as the frame image, in units of counts
        log.info("IDL> sky= mrdfits(framename,2)")
        log.info("IDL> simg= interpolate(sky.allsky, sky.xinterp, sky.yinterp, /grid)")
        log.info("")

        # 3. read in calibration, and expand to full image size; this returns
        # a calibration image the same size as the frame image, in units of
        # nanomaggies per count
        log.info("IDL> calib= mrdfits(framename,1)")
        log.info("IDL> cimg= calib#replicate(1.,nrowc)")
        log.info("")

        ##

        # If you have performed the above calculations, you can return the image to very close to the state it
        # was in when input into the photometric pipeline, as follows:

        # 4. Convert to

        log.info("IDL> dn= img/cimg+simg")
        log.info("")

        log.warning("Then save the DN image as a FITS file!")

        log.info("IDL> writefits,'path',dn")


        # Loop over the downloaded "frame" files
        for path in fs.files_in_path(raw_path, extension="fits"): # extension must be specified because there is also the meta.dat and overlap.dat files!!

            # Open the HDUList
            hdulist = open_fits(path)

            # Get calibrated image in nanomaggies/pixel
            img = hdulist[0].data
            nrowc = img.shape[0] # ysize    (x = columns, y = rows)

            # Get the sky HDU (header)
            sky = hdulist[2].data
            allsky = sky.field("allsky")[0] # for some reason, "allsky" has 3 axes (1, 192, 256), therefore [0]
            xinterp = sky.field("xinterp")[0] # columns (xsize = 2048)    # for some reason, "xinterp" has 2 axes (1, 2048)
            yinterp = sky.field("yinterp")[0] # rows    (ysize = 1489)    # for some reason, "yinterp" has 2 axes (1, 1489)

            # Split x, y and z values that are not masked
            #x_values, y_values, z_values = split_xyz(allsky, arrays=True)

            allsky_xrange = np.arange(allsky.shape[1], dtype=float)
            allsky_yrange = np.arange(allsky.shape[0], dtype=float)

            # INTERPOLATE SKY
            skyf = interp2d(allsky_xrange, allsky_yrange, allsky)
            simg = skyf(xinterp, yinterp)

            # Get the calibration HDU (header)
            calib = hdulist[1].data
            replicate = np.ones(nrowc)
            cimg = np.outer(replicate, calib)

            # DECALIBRATE
            dn = img / cimg + simg

            # CREATE FRAME AND SAVE

            dn_frame = Frame(dn)
            dn_frame.unit = "count"

            # Set the WCS of THE DN FRAME
            dn_frame.wcs = CoordinateSystem(hdulist[0].header)

            name = fs.name(path)
            dn_path = fs.join(counts_path, name)
            dn_frame.saveto(dn_path)

    # -----------------------------------------------------------------

    def get_gain_and_dark_variance_from_photofield(self, path, band):

        """
        This function ...
        :param path:
        :param band:
        :return:
        """

        # Get the hdulist
        hdulist = open_fits(path)

        # HDU0: Empty Header
        # HDU1: photoField Table

        # Data Model: photoField: https://data.sdss.org/datamodel/files/BOSS_PHOTOOBJ/RERUN/RUN/photoField.html

        column_names = hdulist[1].columns

        tbdata = hdulist[1].data

        gain_2d = tbdata.field("gain")
        dark_variance_2d = tbdata.field("dark_variance")

        gain_1d = gain_2d[:, sdss_filter_number[band]]
        dark_variance_1d = dark_variance_2d[:, sdss_filter_number[band]]

        if not is_constant_array(gain_1d): raise ValueError("Gain 1D not constant: " + str(gain_1d))
        if not is_constant_array(dark_variance_1d): raise ValueError("Dark variance 1D not constant: " + str(dark_variance_1d))

        gain = gain_1d[0]
        dark_variance = dark_variance_1d[0]

        # Return the values
        return gain, dark_variance

    # -----------------------------------------------------------------

    def make_sdss_rebinned_frames_in_counts(self, galaxy_name, band, output_path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Making SDSS rebinned frames in counts for " + galaxy_name + " for SDSS " + band + " band ...")

        # Determine the path to the temporary directory for downloading the images
        temp_path = fs.join(fs.home(), time.unique_name("SDSS_primary_fields_" + galaxy_name + "_" + band))

        # Create the temporary directory
        fs.create_directory(temp_path)

        # Download the FITS files to be used for mosaicing
        self.download_sdss_primary_fields_for_galaxy_for_mosaic(galaxy_name, band, temp_path)

        # ----

        # Get the target header
        header = self.get_header_for_galaxy(galaxy_name, "SDSS")

        # To coordinate system
        rebin_wcs = CoordinateSystem(header)

        # Loop over the images, convert to count and rebin
        for path, name in fs.files_in_path(temp_path, extension="fits", returns=["path", "name"]):

            # Debugging
            log.debug("Loading the " + name + " frame ...")

            # Load the frame
            frame = Frame.from_file(path, add_meta=True)

            # Print unit
            log.debug("Unit:" + str(frame.unit))

            # Get the 'NMGY' parameter
            nanomaggy_per_count = frame.meta["nmgy"]

            # Debugging
            log.debug("Converting unit to count / sr...")

            # Convert the frame to count
            frame /= nanomaggy_per_count
            frame.unit = "count"

            # Convert the frame to count/sr
            frame /= frame.pixelarea.to("sr").value
            frame.unit = "count/sr"

            # Debugging
            log.debug("Rebinning to the target coordinate system ...")

            # Rebin the frame
            footprint = frame.rebin(rebin_wcs, exact=False)

            # Debugging
            log.debug("Converting the unit to count ...")

            # Convert the frame from count/sr to count
            frame *= frame.pixelarea.to("sr").value
            frame.unit = "count"

            # Debugging
            log.debug("Adding the rebinned field and footprint to the same image ...")

            # Create image, add image frame and footprint
            image = Image()
            image.add_frame(frame, "primary")
            image.add_frame(footprint, "footprint")
            # padded = Mask(footprint < 0.9).disk_dilation(radius=10)
            # image.add_mask(padded, "padded")

            # Determine new path
            new_path = fs.join(output_path, name + ".fits")

            # Debugging
            log.debug("Saving the image to " + new_path + " ...")

            # Save
            image.saveto(new_path)

        # Remove the temporary directory
        #fs.remove_directory(temp_path)

    # -----------------------------------------------------------------

    def make_sdss_mosaic(self, galaxy_name, band, output_path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Making mosaic for " + galaxy_name + " for SDSS " + band + " band ...")

        # Determine the path to the temporary directory for downloading the images
        temp_path = fs.join(fs.home(), time.unique_name("SDSS_" + galaxy_name + "_" + band))

        # Create the temporary directory
        fs.create_directory(temp_path)

        # Download the FITS files to be used for mosaicing
        self.download_sdss_primary_fields_for_galaxy_for_mosaic(galaxy_name, band, temp_path)

        # Determine the path to the result file
        mosaic_path = fs.join(self.temp_path, "mosaic.fits")

        # Get the target header
        header = self.get_header_for_galaxy(galaxy_name, "SDSS")

        # Save the header to the temp location
        header_path = fs.join(self.temp_path, "target.hdr")
        header.totextfile(header_path)

        # Inform the user
        log.info("Performing the mosaicing ...")

        # Set debug level
        debug_level = 4 if log.is_debug() else 0

        # Perform the mosaicing
        montage.commands.mExec("SDSS", band, raw_dir=temp_path, level_only=False, corners=False, debug_level=debug_level,
                               output_image=mosaic_path, region_header=header_path, workspace_dir=output_path)

        # Load the mosaic image
        #return Image.from_file(mosaic_path)

    # -----------------------------------------------------------------

    def download_sdss_primary_fields_for_galaxy_for_mosaic(self, galaxy_name, band, temp_path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param temp_path:
        :return:
        """

        # Inform the user
        log.info("Getting the overlap of the SDSS observations for " + galaxy_name + " in the " + band + " band ...")

        # To use Montage to find out what images cover a given part of the sky, you first want to run mImgTable on the folder
        # containing your image files, like so:
        # montage_wrapper.commands.mImgtbl('/folder/with/maps/', 'Image_Metadata_Table.dat', corners=True)

        # The 'Image_Metadata_Table.dat' output file is a table recording the geometry of all the files in the folder in
        # question. You then use this file as an input to the function mCoverageCheck, which is used like this:

        # montage_wrapper.commands_extra.mCoverageCheck('Image_Metadata_Table.dat', 'Overlap_Table.dat', mode='box', ra=ra,
        # dec=dec, width=width)

        # The 'Overlap_Table.dat' contains the subset of the rows from 'Image_Metadata_Table.dat' that intersect the defined
        # region (not just those that cover the centre coordinate). I then read in this in with:

        # overlap_files = np.genfromtxt('Overlap_Table.dat', skip_header=3, usecols=[31], dtype=('S500'))

        # Download the SDSS primary fields
        urls = self.download_sdss_primary_fields_for_galaxy(galaxy_name, band, temp_path)

        # Get the coordinate range first for this galaxy
        ra, dec, width = self.get_cutout_range_for_galaxy(galaxy_name)
        ra = ra.to("deg").value
        dec = dec.to("deg").value
        width = width.to("deg").value

        # Get the names of the overlapping image files
        overlap_files = self.get_overlapping_files(temp_path, ra, dec, width)

        # Loop over the FITS files in the temp directory, remove non-overlapping
        for path in fs.files_in_path(temp_path, extension="fits"):

            # If the path is not in the overlap_files list, remove the FITS file
            if path not in overlap_files:

                log.debug("Removing the '" + fs.name(path) + "' image since it does not overlap with the target area ...")
                fs.remove_file(path)

        # Return the original urls from which the images were downloaded
        return urls

    # -----------------------------------------------------------------

    def get_overlapping_files(self, path, ra, dec, width):

        """
        This function ...
        :param path to the directory with the images
        :param ra:
        :param dec:
        :param width:
        :return:
        """

        # Generate the meta and then overlap file
        meta_path, overlap_path = self.generate_meta_and_overlap_file(path, ra, dec, width)

        # Load the overlap table
        overlap_files = np.genfromtxt(overlap_path, skip_header=3, usecols=[32], dtype="S500")

        # Return the names of the overlapping images
        return overlap_files

    # -----------------------------------------------------------------

    def generate_meta_and_overlap_file(self, path, ra, dec, width):

        """
        This function ...
        :param path:
        :param ra:
        :param dec:
        :param width:
        :return:
        """

        # Path of the meta file
        meta_path = fs.join(path, "meta.dat")

        # Get the image table of which images cover a given part of the sky
        montage.commands.mImgtbl(path, meta_path, corners=True)

        # Path of the overlap file
        overlap_path = fs.join(path, "overlap.dat")

        # Check the coverage for our galaxy
        montage.commands_extra.mCoverageCheck(meta_path, overlap_path, mode='box', ra=ra, dec=dec, width=width)

        # Return the paths to the created files
        return meta_path, overlap_path

    # -----------------------------------------------------------------

    def download_sdss_primary_fields_for_galaxy(self, galaxy_name, band, path):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :param path:
        :return:
        """

        # Inform the user
        log.info("Downloading the SDSS primary fields for " + galaxy_name + " in the " + band + " band to '" + path + "' ...")

        # Get the urls
        urls = self.get_sdss_primary_field_urls_for_galaxy(galaxy_name, band)

        # Debugging
        log.debug("Number of primary fields that will be downloaded: " + str(len(urls)))

        ## NEW: USE DR13 INSTEAD OF DR10  (BETTER CALIBRATION?)
        #new_urls = [url.replace("dr10/boss", "dr13/eboss") for url in urls] # DOESN'T WORK: RESTRICTED ACCESS: ASKS FOR USERNAME AND PASSWORD!!
        ###

        # USE DR12
        new_urls = [url.replace("dr10/boss", "dr12/boss") for url in urls]

        # Download the files
        paths = network.download_files(new_urls, path)

        # Debugging
        log.debug("Decompressing the files ...")

        # Decompress the files and remove the originals
        archive.decompress_files(paths, remove=True)

        # Return the URLS
        return new_urls

    # -----------------------------------------------------------------

    def get_sdss_primary_field_urls_for_galaxy(self, galaxy_name, band):

        """
        This function ...
        :param galaxy_name:
        :param band:
        :return:
        """

        # Get the coordinate range first for this galaxy
        ra, dec, width = self.get_cutout_range_for_galaxy(galaxy_name)

        # Get the SDSS fields that cover this coordinate range (from Montage) (URLS)
        table = self.get_sdss_fields_for_coordinate_range(band, ra, dec, width)

        # Get the set of unique URLs
        sdss_urls = list(set(table["url"]))

        # ...
        dr12_pri = []
        for i in range(len(self.sdss_field_table)):
            entry = str(int(self.sdss_field_table['RUN'][i])) + ' ' + str(int(self.sdss_field_table['CAMCOL'][i])) + ' ' + str(int(self.sdss_field_table['FIELD'][i]))
            dr12_pri.append(entry)

        # Get the URLS of the primary fields
        sdss_urls_pri = filter_sdss_urls_for_primary_frames(sdss_urls, dr12_pri)

        # Return the list of URLS
        return sdss_urls_pri

    # -----------------------------------------------------------------

    def get_sdss_fields_for_coordinate_range(self, band, ra, dec, width):

        """
        This function ...
        :return:
        """

        # Determine the path to the temporary table file
        path = fs.join(self.temp_path, "fields.tbl")

        ra = ra.to("deg").value
        dec = dec.to("deg").value
        width = width.to("deg").value

        # Get the info
        montage.mArchiveList("SDSS", band, str(ra) + " " + str(dec), width, width, path)

        # Load the table
        table = Table.read(path, format="ascii")

        # Return the table
        return table

    # -----------------------------------------------------------------

    def get_header_for_galaxy(self, galaxy_name, instrument):

        """
        This function ...
        :param galaxy_name:
        :param instrument
        :return:
        """

        # Get cutout range
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

        # Determine the path to the temporary header file
        header_path = fs.join(self.temp_path, "header.hdr")

        # Create the header
        montage.commands.mHdr(str(ra) + ' ' + str(dec), width, header_path, pix_size=pix_size)

        # Load the header
        return Header.fromtextfile(header_path)

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

def filter_sdss_urls_for_primary_frames(urls, index):

    urls_pri = []
    for url in urls:
        run = url.split('/')[9].lstrip('0')
        camcol = url.split('/')[10].lstrip('0')
        field = url.split('/')[11].split('.')[0].split('-')[4].lstrip('0')
        check_string = run + ' ' + camcol + ' ' + field
        if check_string in index:
            urls_pri.append(url)
    return urls_pri

# -----------------------------------------------------------------

def is_constant_array(array):

    """
    This function ...
    :param array:
    :return:
    """

    return np.min(array) == np.max(array)

# -----------------------------------------------------------------
