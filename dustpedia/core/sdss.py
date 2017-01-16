#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.sdss Contains the SDSSMosaicMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from .dataprocessing import DustPediaDataProcessing
from ...core.tools import filesystem as fs
from ...core.tools import time
from ...magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

sdss_bands = ["u", "g", "r", "i", "z"]

# -----------------------------------------------------------------

class SDSSMosaicMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(SDSSMosaicMaker, self).__init__(config)

        # The DustPedia data processing object
        self.dpdp = None

        self.root_path = None

        self.band_paths = dict()
        self.counts_paths = dict()
        self.footprints_paths = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Create directories
        self.create_directories()

        # 2. Download
        self.get_frames()

        # Convert to counts
        self.convert_to_counts()

        # Make poisson frames
        self.make_poisson_frames()

        # Rebin
        self.rebin()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SDSSMosaicMaker, self).setup()

        # Create the DustPedia data processing instance
        self.dpdp = DustPediaDataProcessing()

        # Get target header
        self.get_target_header()

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating directories ...")

        # Root path
        self.root_path = fs.join(fs.home(), time.unique_name("SDSS_" + self.config.galaxy_name))
        fs.create_directory(self.root_path)

        # Loop over the bands
        for band in self.config.bands:

            path = fs.create_directory_in(self.root_path, band)

            counts_path = fs.create_directory_in(path, "counts")

            footprints_path = fs.create_directory_in(path, "footprints")

            self.band_paths[band] = path
            self.counts_paths[band] = counts_path
            self.footprints_paths[band] = footprints_path

        # temp_path = fs.join(fs.home(), "SDSS_NGC3031_u_2016-08-10--16-48-33-775")

        ###

        # RAW PATH
        #raw_path = fs.join(temp_path, "raw")  # frames with HDU0, HDU1, HDU2 .. (primary, calib, sky, ...) IN NANOMAGGIES PER PIXEL
        #fs.create_directory(raw_path)

        # FIELDS PATH
        #fields_path = fs.join(temp_path, "fields")  # photoField files
        #fs.create_directory(fields_path)

        # COUNTS PATH
        #counts_path = fs.join(temp_path, "counts")  # frames IN COUNTS (DN)
        #fs.create_directory(counts_path)

        # POISSON PATH
        #poisson_path = fs.join(temp_path, "poisson_nmgy")  # error maps in NANOMAGGIES PER PIXEL
        #fs.create_directory(poisson_path)

        # REBINNED PATH
        #rebinned_path = fs.join(temp_path, "rebinned")  # images with frame and corresponding error map in NANOMAGGIES, REBINNED
        #fs.create_directory(rebinned_path)

        # FOOTPRINTS PATH
        #footprints_path = fs.join(temp_path, "footprints")
        #fs.create_directory(footprints_path)

        # MOSAIC PATH
        #mosaics_path = fs.join(temp_path, "mosaics")
        #fs.create_directory(mosaics_path)

        # RESULT PATH
        #results_path = fs.join(temp_path, "results")
        #fs.create_directory(results_path)

    # -----------------------------------------------------------------

    def get_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading ...")

        # Download the FITS files to be used for mosaicing
        urls = self.download_sdss_primary_fields_for_galaxy_for_mosaic(galaxy_name, band, raw_path)

        #### NEW: DOWNLOAD THE FIELD TABLES

        # dr10 / boss / photoObj / frames

        field_url_start = "http://data.sdss3.org/sas/dr12/env/BOSS_PHOTOOBJ"  # = $BOSS_PHOTOOBJ/RERUN/RUN
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

    # -----------------------------------------------------------------

    def convert_to_counts(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def make_poisson_frames(self):

        """
        This function ...
        :return:
        """

        # Make (Poisson) noise maps
        self.make_sdss_noise_maps_in_nanomaggies(raw_path, counts_path, fields_path, poisson_path)

    # -----------------------------------------------------------------

    def get_target_header(self):

        """
        This function ...
        :return:
        """

        # Get the target header
        header = self.dpdp.get_header_for_galaxy(galaxy_name, "SDSS")

        # To coordinate system
        self.rebin_wcs = CoordinateSystem(header)

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # DO REBINNING, CREATE IMAGES WITH REBINNED PRIMARY AND ERROR FRAME IN NANOMAGGIES PER PIXEL, AND FOOTPRINT FILES
        self.rebin_sdss_frames_and_error_maps(rebin_wcs, raw_path, poisson_path, rebinned_path, footprints_path)

        ## Make the footprints: WAS ONLY NECESSARY FOR WHEN FOOTPRINTS WERE NOT CREATED DURING FUNCTION ABOVE
        ##self.make_sdss_footprints(rebinned_path, footprints_path)

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # CONVERT TO JY/PIX
        self.convert_sdss_mosaic_and_error_map_to_jansky(galaxy_name, band, mosaics_path, results_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

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
        relerrors[relerrors < 0.] = 0.0  # set negative values for relative error map to zero
        relerrors.replace_nans(0.0)  # set NaN values (because mosaic was zero) to zero

        # Save mosaic as FITS file
        mosaic_output_path = fs.join(output_path, id_string + ".fits")
        mosaic.saveto(mosaic_output_path)

        # Save error map as FITS file
        errors_output_path = fs.join(output_path, id_string + "_errors.fits")
        mosaic_errors.saveto(errors_output_path)

        # Save relative error map as FITS file
        relerrors_output_path = fs.join(output_path, id_string + "_relerrors.fits")
        relerrors.saveto(relerrors_output_path)

# -----------------------------------------------------------------

def download_sdss_primary_fields_for_galaxy_for_mosaic(galaxy_name, band, temp_path):

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
    urls = download_sdss_primary_fields_for_galaxy(galaxy_name, band, temp_path)

    # Get the coordinate range first for this galaxy
    ra, dec, width = get_cutout_range_for_galaxy(galaxy_name)
    ra = ra.to("deg").value
    dec = dec.to("deg").value
    width = width.to("deg").value

    # Get the names of the overlapping image files
    overlap_files = get_overlapping_files(temp_path, ra, dec, width)

    # Loop over the FITS files in the temp directory, remove non-overlapping
    for path in fs.files_in_path(temp_path, extension="fits"):

        # If the path is not in the overlap_files list, remove the FITS file
        if path not in overlap_files:

            log.debug("Removing the '" + fs.name(path) + "' image since it does not overlap with the target area ...")
            fs.remove_file(path)

    # Return the original urls from which the images were downloaded
    return urls

# -----------------------------------------------------------------

def download_sdss_primary_fields_for_galaxy(galaxy_name, band, path):

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
    urls = get_sdss_primary_field_urls_for_galaxy(galaxy_name, band)

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
    ra, dec, width = get_cutout_range_for_galaxy(galaxy_name)

    # Get the SDSS fields that cover this coordinate range (from Montage) (URLS)
    table = get_sdss_fields_for_coordinate_range(band, ra, dec, width)

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

def get_sdss_fields_for_coordinate_range(band, ra, dec, width):

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
