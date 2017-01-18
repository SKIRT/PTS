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

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from .dataprocessing import DustPediaDataProcessing
from .sample import DustPediaSample
from ...core.tools import filesystem as fs
from ...core.tools import time
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...core.tools import network
from ...magic.basics.coordinate import SkyCoordinate
from ...core.tools import archive

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

        # The DustPedia sample object
        self.sample = None

        # The DustPedia data processing object
        self.dpdp = None

        # The NGC name
        self.ngc_name = None

        # The cutout properties
        self.cutout_center = None
        self.cutout_width = None

        # The WCS
        self.rebin_wcs = None

        # Temporary root path
        self.root_path = None

        # Paths for different bands
        self.fields_paths = dict()
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

        # 2. Create directories
        self.create_directories()

        # 3. Get the cutout range
        self.get_range()

        # Get target header
        self.get_target_header()

        # 3. Download fields
        self.download_fields()

        # Convert to counts (decalibrate)
        self.convert_to_counts()

        # 5. Convert to counts
        self.convert_to_nanomaggies()

        # 7. Rebin
        self.rebin()

        # 8. Show
        self.show()

        # 9. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SDSSMosaicMaker, self).setup()

        # Create the DustPedia sample object
        self.sample = DustPediaSample()

        # Create the DustPedia data processing instance
        self.dpdp = DustPediaDataProcessing()

        # Get the NGC name
        self.ngc_name = self.sample.get_name(self.config.galaxy_name)

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

            fields_path = fs.create_directory_in(path, "fields")
            counts_path = fs.create_directory_in(path, "counts")
            footprints_path = fs.create_directory_in(path, "footprints")

            self.band_paths[band] = path
            self.fields_paths[band] = fields_path
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

    def get_range(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

        # Get the coordinate range first for this galaxy
        ra, dec, width = self.dpdp.get_cutout_range_for_galaxy(self.ngc_name)
        #ra = ra.to("deg").value
        #dec = dec.to("deg").value
        #width = width.to("deg").value

        self.cutout_center = SkyCoordinate(ra, dec, unit="deg", frame="fk5")
        self.cutout_width = width

    # -----------------------------------------------------------------

    def get_target_header(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

        # Get the target header
        header = self.dpdp.get_header_for_galaxy(self.ngc_name, "SDSS")

        # To coordinate system
        self.rebin_wcs = CoordinateSystem(header)

    # -----------------------------------------------------------------

    def download_fields(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading the fields ...")

        # Download the FITS files to be used for mosaicing
        urls = download_sdss_primary_fields_for_galaxy_for_mosaic(galaxy_name, band, raw_path)

        #### NEW: DOWNLOAD THE FIELD TABLES

        # dr10 / boss / photoObj / frames

        # Inform the user
        #log.info("Downloading the field tables ...")

        # Start of URLs
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
        self.field_paths = network.download_files(field_urls, fields_path)

    # -----------------------------------------------------------------

    def convert_to_counts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting to counts ...")

        # Loop over the downloaded "frame" files
        for path in fs.files_in_path(raw_path, extension="fits"):  # extension must be specified because there is also the meta.dat and overlap.dat files!!

            # Open the HDUList
            hdulist = open_fits(path)

            # Get calibrated image in nanomaggies/pixel
            img = hdulist[0].data
            nrowc = img.shape[0]  # ysize    (x = columns, y = rows)

            # Get the sky HDU (header)
            sky = hdulist[2].data
            allsky = sky.field("allsky")[0]  # for some reason, "allsky" has 3 axes (1, 192, 256), therefore [0]
            xinterp = sky.field("xinterp")[
                0]  # columns (xsize = 2048)    # for some reason, "xinterp" has 2 axes (1, 2048)
            yinterp = sky.field("yinterp")[
                0]  # rows    (ysize = 1489)    # for some reason, "yinterp" has 2 axes (1, 1489)

            # Split x, y and z values that are not masked
            # x_values, y_values, z_values = split_xyz(allsky, arrays=True)

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

    def convert_to_nanomaggies(self):

        """
        This function ...
        :return:
        """

        # Make (Poisson) noise maps
        make_sdss_noise_maps_in_nanomaggies(raw_path, counts_path, fields_path, poisson_path)

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning ...")

        # DO REBINNING, CREATE IMAGES WITH REBINNED PRIMARY AND ERROR FRAME IN NANOMAGGIES PER PIXEL, AND FOOTPRINT FILES
        rebin_sdss_frames_and_error_maps(rebin_wcs, raw_path, poisson_path, rebinned_path, footprints_path)

        ## Make the footprints: WAS ONLY NECESSARY FOR WHEN FOOTPRINTS WERE NOT CREATED DURING FUNCTION ABOVE
        ##self.make_sdss_footprints(rebinned_path, footprints_path)

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # CONVERT TO JY/PIX
        convert_sdss_mosaic_and_error_map_to_jansky(galaxy_name, band, mosaics_path, results_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        ## WRITE RESULT TO OUTPUT DIRECTORY

        # Load the image
        id_string = self.ngc_name + '_SDSS_' + band
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
    #ra, dec, width = get_cutout_range_for_galaxy(galaxy_name)
    #ra = ra.to("deg").value
    #dec = dec.to("deg").value
    #width = width.to("deg").value

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

def make_sdss_noise_maps_in_nanomaggies(raw_path, counts_path, fields_path, poisson_path):

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
        gain, dark_variance = get_gain_and_dark_variance_from_photofield(field_file_path, band)

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

def get_gain_and_dark_variance_from_photofield(path, band):

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
