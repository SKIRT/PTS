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
from scipy.interpolate import interp2d

# Import astronomical modules
from astropy.io.fits import open as open_fits

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
from ...magic.tools import mosaic
from ...core.basics.filter import Filter
from ...core.tools import formatting as fmt
from ...magic.core.frame import Frame, sum_frames, sum_frames_quadratically

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
        self.band_paths = dict()
        self.fields_paths = dict()
        self.photofield_paths = dict()
        self.counts_paths = dict()
        self.poisson_paths = dict()
        self.rebinned_paths = dict()
        self.mosaics_paths = dict()

        # The URLs from where the images were downloaded
        self.urls = dict()

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

        # 4. Get target header
        self.get_target_header()

        # 5. Get the fields
        self.get_fields()

        # 6. Convert to counts (decalibrate)
        self.convert_to_counts()

        # 7. Convert to counts
        self.convert_to_nanomaggies()

        # 8. Rebin
        self.rebin()

        # 9. Make the mosaics
        self.mosaic()

        # 10. Convert to desired unit
        self.convert_units()

        # 11. Show
        self.show()

        # 12. Write
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

            # Create root directory
            path = fs.create_directory_in(self.root_path, band)

            # Create directories
            fields_path = fs.create_directory_in(path, "fields") # frames with HDU0, HDU1, HDU2 .. (primary, calib, sky, ...) IN NANOMAGGIES PER PIXEL
            photofield_path = fs.create_directory_in(path, "photofield") # photoField files
            counts_path = fs.create_directory_in(path, "counts") # frames IN COUNTS (DN)
            poisson_path = fs.create_directory_in(path, "poisson") # error maps in NANOMAGGIES PER PIXEL
            rebinned_path = fs.create_directory_in(path, "rebinned") # images with frame and corresponding error map in NANOMAGGIES, REBINNED
            mosaics_path = fs.create_directory_in(path, "mosaics") # the mosaic images

            # Set paths
            self.band_paths[band] = path
            self.fields_paths[band] = fields_path
            self.photofield_paths[band] = photofield_path
            self.counts_paths[band] = counts_path
            self.poisson_paths[band] = poisson_path
            self.rebinned_paths[band] = rebinned_path
            self.mosaics_paths[band] = mosaics_path

    # -----------------------------------------------------------------

    def get_range(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting cutout range ...")

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

    def get_fields(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the fields ...")

        # Download SDSS fields
        self.download_fields()

        # Filter those not overlapping the target position
        self.filter_non_overlapping()

        # Download photofield data
        self.download_photofield()

    # -----------------------------------------------------------------

    def download_fields(self):

        """
        This function ...
        :return:
        """

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

        # Get the coordinate range first for this galaxy
        # ra, dec, width = get_cutout_range_for_galaxy(galaxy_name)
        # ra = ra.to("deg").value
        # dec = dec.to("deg").value
        # width = width.to("deg").value

        # Inform the user
        #log.info("Downloading fields ...")

        # Loop over the bands
        for band in self.config.bands:

            download_to = self.fields_paths[band]

            # Inform the user
            log.info("Downloading the SDSS primary fields for " + self.ngc_name + " in the " + band + " band to '" + download_to + "' ...")

            # Get the urls
            urls = self.get_sdss_primary_field_urls_for_galaxy(band)

            # Debugging
            log.debug("Number of primary fields that will be downloaded: " + str(len(urls)))

            ## NEW: USE DR13 INSTEAD OF DR10  (BETTER CALIBRATION?)
            # new_urls = [url.replace("dr10/boss", "dr13/eboss") for url in urls] # DOESN'T WORK: RESTRICTED ACCESS: ASKS FOR USERNAME AND PASSWORD!!
            ###

            # USE DR12
            new_urls = [url.replace("dr10/boss", "dr12/boss") for url in urls]

            # Download the files
            paths = network.download_files(new_urls, download_to)

            # Debugging
            log.debug("Decompressing the files ...")

            # Decompress the files and remove the originals
            archive.decompress_files(paths, remove=True)

            # Return the URLS
            #return new_urls

            # Save the URLS for this band
            # The original urls from which the images were downloaded
            self.urls[band] = new_urls

    # -----------------------------------------------------------------

    def get_sdss_primary_field_urls_for_galaxy(self, band):

        """
        This function ...
        :param band:
        :return:
        """

        # Inform the user
        log.info("Getting the URLs of the primary field for galaxy " + self.ngc_name + " for the " + band + " band ...")

        # Get the SDSS fields that cover this coordinate range (from Montage) (URLS)
        table = mosaic.get_field_table(self.cutout_center, self.cutout_width, band)

        # Get the set of unique URLs
        sdss_urls = list(set(table["url"]))

        # ...
        dr12_pri = []
        for i in range(len(self.dpdp.sdss_field_table)):
            entry = str(int(self.dpdp.sdss_field_table['RUN'][i])) + ' ' + str(int(self.dpdp.sdss_field_table['CAMCOL'][i])) + ' ' + str(int(self.dpdp.sdss_field_table['FIELD'][i]))
            dr12_pri.append(entry)

        # Get the URLS of the primary fields
        sdss_urls_pri = filter_sdss_urls_for_primary_frames(sdss_urls, dr12_pri)

        # Return the list of URLS
        return sdss_urls_pri

    # -----------------------------------------------------------------

    def filter_non_overlapping(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking overlap of the SDSS primary fields ...")

        # Loop over the bands
        for band in self.config.bands:

            # Path where the fields are
            fields_path = self.fields_paths[band]

            # Inform the user
            log.info("Getting the overlap of the SDSS observations for " + self.ngc_name + " in the " + band + " band ...")

            # Generate meta file
            meta_path = mosaic.generate_meta_file(fields_path)

            # Generate overlap file
            overlap_path = mosaic.generate_overlap_file(fields_path, self.cutout_center.ra, self.cutout_center.dec, self.cutout_width, meta_path)

            # Get the names of the overlapping image files
            overlap_files = np.genfromtxt(overlap_path, skip_header=3, usecols=[32], dtype="S500")

            # Loop over the FITS files in the temp directory, remove non-overlapping
            for path in fs.files_in_path(fields_path, extension="fits"):

                # If the path is not in the overlap_files list, remove the FITS file
                if path not in overlap_files:
                    log.debug("Removing the '" + fs.name(path) + "' image since it does not overlap with the target area ...")
                    fs.remove_file(path)

    # -----------------------------------------------------------------

    def download_photofield(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading the photofield tables ...")

        # Start of URLs
        field_url_start = "http://data.sdss3.org/sas/dr12/env/BOSS_PHOTOOBJ"  # = $BOSS_PHOTOOBJ/RERUN/RUN
        # NOT WORKING WITH DR13 : BECAUSE OF PERMISSION ISSUES (PASSWORD IS ASKED WHEN ENTERED IN BROWSER !!)

        # FIELD URL = field_url_start / RERUN / RUN    + / photoField-6digits-CAMCOL.fits
        # example: http://data.sdss3.org/sas/dr12/env/BOSS_PHOTOOBJ/301/4294/photoField-004294-5.fits

        # FRAME URL = $BOSS_PHOTOOBJ / frames / RERUN / RUN / CAMCOL    +   /frame-[ugriz]-6digits-CAMCOL-FRAMESEQ.fits.bz2
        # example: http://data.sdss3.org/sas/dr10/boss/photoObj/frames/301/4294/5/frame-i-004294-5-0229.fits.bz2

        # Loop over the bands
        for band in self.config.bands:

            # Inform the user
            log.info("Downloading the photofield tables for the " + band + " band ...")

            # Initialize a list to contain the URLs of the field tables
            field_urls = []
            for url in self.urls[band]:

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
            field_paths = network.download_files(field_urls, self.photofield_paths[band])

    # -----------------------------------------------------------------

    def convert_to_counts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting the fields to counts (decalibrating) ...")

        # Loop over the bands
        for band in self.config.bands:

            fields_path = self.fields_paths[band]
            counts_path = self.counts_paths[band]

            # Loop over the downloaded "frame" files
            for path in fs.files_in_path(fields_path, extension="fits"):  # extension must be specified because there is also the meta.dat and overlap.dat files!!

                # Open the HDUList
                hdulist = open_fits(path)

                # Get calibrated image in nanomaggies/pixel
                img = hdulist[0].data
                nrowc = img.shape[0]  # ysize    (x = columns, y = rows)

                # Get the sky HDU (header)
                sky = hdulist[2].data
                allsky = sky.field("allsky")[0]  # for some reason, "allsky" has 3 axes (1, 192, 256), therefore [0]
                xinterp = sky.field("xinterp")[0]  # columns (xsize = 2048)    # for some reason, "xinterp" has 2 axes (1, 2048)
                yinterp = sky.field("yinterp")[0]  # rows    (ysize = 1489)    # for some reason, "yinterp" has 2 axes (1, 1489)

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

        # Inform the user
        log.info("Converting fields to nanomaggies ...")

        # Loop over the bands
        for band in self.config.bands:

            # Inform the user
            log.info("Computing error maps for each SDSS frame in nanomaggies per pixel ...")

            fields_path = self.fields_paths[band]
            photofields_path = self.photofield_paths[band]
            counts_path = self.counts_paths[band]
            poisson_path = self.poisson_paths[band]

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
                field_file_path = fs.join(photofields_path, "photoField-" + digits + "-" + camcol + ".fits")

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
                raw_frame_path = fs.join(fields_path, name + ".fits")
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

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning frames ...")

        # Loop over the bands
        for band in self.config.bands:

            # Get paths
            fields_path = self.fields_paths[band]
            poisson_path = self.poisson_paths[band]
            rebinned_path = self.rebinned_paths[band]

            # DO REBINNING, CREATE IMAGES WITH REBINNED PRIMARY AND ERROR FRAME IN NANOMAGGIES PER PIXEL, AND FOOTPRINT FILES

            # Inform the user
            log.info("Rebinning primary SDSS frames and corresponding error maps to the target coordinate system ...")

            # Loop over the files in the raw directory
            for path, name in fs.files_in_path(fields_path, extension="fits", returns=["path", "name"]):

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
                footprint = frame.rebin(self.rebin_wcs, exact=False)

                # REBIN THE ERROR MAP
                errormap.rebin(self.rebin_wcs, exact=False)

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

                # Add frames
                image.add_frame(frame, "primary")
                image.add_frame(errormap, "noise")
                image.add_frame(footprint, "footprint")

                # Determine path for the image
                image_path = fs.join(rebinned_path, name + ".fits")
                image.saveto(image_path)

    # -----------------------------------------------------------------

    def mosaic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the mosaicing ...")

        # Loop over the bands
        for band in self.config.bands:

            rebinned_path = self.rebinned_paths[band]
            mosaics_path = self.mosaics_paths[band]

            # Inform the user
            log.info("Creating mosaic image for the " + band + " band ...")

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
                #footprint_path = fs.join(footprints_path, name + ".fits")

                # Open the footprint
                #footprint = Frame.from_file(footprint_path)

                # Open the image
                image = Image.from_file(path)

                # Get the footprint
                footprint = image.frames.footprint

                # Get the a and b frame
                a = image.frames.primary  # IN NANOMAGGIES PER PIXEL
                b = image.frames.noise  # IN NANOMAGGIES PER PIXEL

                # SET NANS TO ZERO
                a.replace_nans(0.0)
                b.replace_nans(0.0)

                # Add to the list
                primary_frames.append(a)
                error_frames.append(b)
                n_overlapping_frames.append(footprint)

            # Calculate n_overlapping_frames array
            n_overlapping = sum_frames(*n_overlapping_frames)  # = n_overlapping_frames[x,y] = 2D ARRAY !

            # CALCULATE THE MOSAIC FRAME IN NANOMAGGIES
            mosaic_frame = sum_frames(*primary_frames) / n_overlapping
            mosaic_frame.wcs = self.rebin_wcs

            # CALCULATE THE MOSAIC ERROR MAP IN NANOMAGGIES
            mosaic_errormap = sum_frames_quadratically(*error_frames) / n_overlapping
            mosaic_errormap.wcs = self.rebin_wcs

            # SAVE THE MOSAIC IN NANOMAGGIES
            mosaic_path = fs.join(mosaics_path, "mosaic_nanomaggy.fits")
            #mosaic_frame.saveto(mosaic_path)

            # SAVE THE MOSAIC ERROR MAP IN NANOMAGGIES
            #mosaic_error_path = fs.join(mosaics_path, "mosaic_errors.fits")
            #mosaic_errormap.saveto(mosaic_error_path)

            # Create image
            mosaic = Image()
            mosaic.add_frame(mosaic_frame, "primary")
            mosaic.add_frame(mosaic_errormap, "errors")

            # Save
            mosaic.saveto(mosaic_path)

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting to Jansky ...")

        # Loop over the bands
        for band in self.config.bands:

            mosaics_path = self.mosaics_paths[band]

            # DETERMINE ID STRING TO SAVE THE RESULT
            id_string = self.ngc_name + '_SDSS_' + band

            # Inform the user
            log.info("Converting the SDSS mosaic and error map to Jansky per pixel ...")

            # 1 nanomaggy = approximately 3.613e-6 Jy

            # Open the mosaic
            mosaic_path = fs.join(mosaics_path, "mosaic.fits")
            #mosaic_frame = Frame.from_file(mosaic_path)

            mosaic = Image.from_file(mosaic_path)

            # Open the mosaic error map
            #mosaic_error_path = fs.join(mosaics_path, "mosaic_errors.fits")
            #mosaic_errors = Frame.from_file(mosaic_error_path)

            # DO THE CONVERSION FOR THE MOSAIC, SET NEW UNIT
            #mosaic *= 3.613e-6
            #mosaic.unit = "Jy / pix"

            mosaic.convert_to("Jy")

            # DO THE CONVERSION FOR THE MOSAIC ERROR MAP, SET NEW UNIT
            #mosaic_errors *= 3.613e-6
            #mosaic_errors.unit = "Jy / pix"

            # CREATE IMAGE
            #image = Image()
            #image.add_frame(mosaic, "primary")
            #image.add_frame(mosaic_errors, "errors")

            # SAVE THE MOSAIC IMAGE (with error map) IN JANSKY
            #result_path = fs.join(results_path, id_string + ".fits")
            #image.saveto(result_path)

            # Determine new path and save
            new_path = fs.join(mosaics_path, "mosaic_jansky.fits")
            mosaic.saveto(new_path)

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finishing ...")

        # Loop over the bands
        for band in self.config.bands:

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
            #mosaic_output_path = fs.join(output_path, id_string + ".fits")
            mosaic_output
            mosaic.saveto(mosaic_output_path)

            # Save error map as FITS file
            errors_output_path = fs.join(output_path, id_string + "_errors.fits")
            mosaic_errors.saveto(errors_output_path)

            # Save relative error map as FITS file
            relerrors_output_path = fs.join(output_path, id_string + "_relerrors.fits")
            relerrors.saveto(relerrors_output_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        print("")

        # Loop over the bands
        for band in self.config.bands:

            # Construct filter
            fltr = Filter.sdss(band)

            print(fmt.green + fmt.bold + str(fltr.description()))
            print("")

            print(" - Number of fields: " + str(len(self.urls[band])))
            #print(" - Average poisson error: " + )

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

def filter_sdss_urls_for_primary_frames(urls, index):

    """
    This function ...
    :param urls:
    :param index:
    :return:
    """

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

def convert_sdss_mosaic_and_error_map_to_jansky(self, galaxy_name, band, mosaics_path, results_path):

    """
    This function ...
    :param galaxy_name:
    :param band:
    :param mosaics_path:
    :param results_path:
    :return:
    """



# -----------------------------------------------------------------
