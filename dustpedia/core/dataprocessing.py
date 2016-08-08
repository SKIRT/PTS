#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.dataprocessing Contains the DustPediaDataProcessing class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import shutil
import tempfile
import numpy as np

# Import astronomical modules
import montage_wrapper as montage
from astropy.table import Table
from astropy.utils import lazyproperty
from astropy.units import Unit
from astropy.io.fits import Header, getheader
#from astroquery.sdss import SDSS

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
from ...magic.core.frame import Frame, sum_frames
from ...magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

dustpedia_dat_path = fs.join(introspection.pts_dat_dir("dustpedia"))

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

        # Inform the user
        log.info("Creating directories ...")

        # Determine the path to the temporary directory for downloading the images
        working_path = fs.join(fs.home(), time.unique_name("GALEX_" + galaxy_name))

        #working_path = fs.join(fs.home(), "GALEX_NGC3031_2016-08-04--10-59-27-047")

        # Create the working directory
        fs.create_directory(working_path)

        # DOWNLOAD PATH
        download_path = fs.join(working_path, "download")
        # Create download directory
        fs.create_directory(download_path)

        # RESPONSE AND BACKGROUND PATH
        response_path = fs.join(working_path, "response")
        background_path = fs.join(working_path, "background")
        fs.create_directories(response_path, background_path)

        # COUNT PATH
        counts_path = fs.join(working_path, "counts")
        fs.create_directory(counts_path)

        # RAW PATH
        raw_path = fs.join(working_path, "raw")
        # Create raw directory
        fs.create_directory(raw_path)

        # TEMP PATH
        temp_path = fs.join(working_path, "temp")
        temp_fuv_path = fs.join(temp_path, "FUV")
        temp_nuv_path = fs.join(temp_path, "NUV")
        # Create temp directory
        fs.create_directory(temp_path)
        fs.create_directories(temp_fuv_path, temp_nuv_path)

        # 1 and 2 RAW directories
        raw_fuv_path = fs.join(raw_path, "FUV")
        raw_nuv_path = fs.join(raw_path, "NUV")
        fs.create_directories(raw_fuv_path, raw_nuv_path)

        # download/images, download/response and download/background
        download_images_path = fs.join(download_path, "images")
        download_response_path = fs.join(download_path, "reponse")
        download_background_path = fs.join(download_path, "background")

        download_counts_path = fs.join(download_path, "counts")
        fs.create_directory(download_counts_path)

        fs.create_directories(download_images_path, download_response_path, download_background_path)


        #

        # Download the GALEX observations to the temporary directory  # they are decompressed here also
        self.download_galex_observations_for_galaxy(galaxy_name, download_images_path, download_response_path, download_background_path, download_counts_path)


        # FUV and NUV response directories
        response_fuv_path = fs.join(response_path, "FUV")
        response_nuv_path = fs.join(response_path, "NUV")
        fs.create_directories(response_fuv_path, response_nuv_path)

        # FUV and NUV background directories
        background_fuv_path = fs.join(background_path, "FUV")
        background_nuv_path = fs.join(background_path, "NUV")
        fs.create_directories(background_fuv_path, background_nuv_path)

        # FUV AND NUV counts directories
        counts_fuv_path = fs.join(counts_path, "FUV")
        counts_nuv_path = fs.join(counts_path, "NUV")
        fs.create_directories(counts_fuv_path, counts_nuv_path)

        ####

        # Inform the user
        log.info("Splitting observations into FUV and NUV ...")

        # Split downloaded images into FUV and NUV
        self.split_galex_observations(download_images_path, raw_fuv_path, raw_nuv_path)

        # Split response maps into FUV and NUV
        self.split_galex_observations(download_response_path, response_fuv_path, response_nuv_path)

        # Split background maps into FUV and NUV
        self.split_galex_observations(download_background_path, background_fuv_path, background_nuv_path)

        # Split count maps into FUV and NUV
        self.split_galex_observations(download_counts_path, counts_fuv_path, counts_nuv_path)

        ###

        # AFTER SPLIT, THIS CAN BE DONE:

        fs.remove_directory(download_images_path)
        fs.remove_directory(download_response_path)
        fs.remove_directory(download_background_path)
        fs.remove_directory(download_counts_path)
        fs.remove_directory(download_path)

        ###

        # Inform the user
        log.info("Getting cutout range ...")

        # Get coordinate range for target image
        ra, dec, width = self.get_cutout_range_for_galaxy(galaxy_name)

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
            if sum(1 for line in open(overlap_path)) <= 3: log.warning("No GALEX " + band + " coverage for " + galaxy_name)
            else: bands_in_dict[band] = bands_dict[band]

            # Remove overlap file
            fs.remove_file(overlap_path)

        # Check if coverage in any band
        if len(bands_in_dict) == 0: raise RuntimeError("No coverage in any GALEX band")

        # Loop over bands, conducting SWarping function
        for band in bands_in_dict.keys():

            # Mosaicing
            mosaic_galex(galaxy_name, ra, dec, width, bands_dict[band], working_path, temp_path, metadata_paths[band], output_path)

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
        image.save(path)

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
            image.save(new_path)

        # Remove the temporary directory
        fs.remove_directory(temp_path)

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
        self.download_sdss_primary_fields_for_galaxy(galaxy_name, band, temp_path)

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

        meta_path = fs.join(path, "meta.dat")

        # Get the image table of which images cover a given part of the sky
        montage.commands.mImgtbl(path, meta_path, corners=True)

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

        # Download the files
        paths = network.download_files(urls, path)

        # Debugging
        log.debug("Decompressing the files ...")

        # Decompress the files and remove the originals
        archive.decompress_files(paths, remove=True)

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
        sdss_urls_pri = SDSS_Primary_Check(sdss_urls, dr12_pri)

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

def SDSS_Primary_Check(urls, index):

    urls_pri = []
    for url in urls:
        run = url.split('/')[9].lstrip('0')
        camcol = url.split('/')[10].lstrip('0')
        field = url.split('/')[11].split('.')[0].split('-')[4].lstrip('0')
        check_string = run+' '+camcol+' '+field
        if check_string in index:
            urls_pri.append(url)
    return urls_pri

# -----------------------------------------------------------------
