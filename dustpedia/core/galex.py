#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.sdss Contains the GALEXMosaicMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import lmfit
import numpy as np
from scipy.ndimage import interpolation

# Import astronomical modules
from astropy.io.fits import getheader, PrimaryHDU, HDUList, writeto, Header
from astropy.io.fits import open as open_fits
from astropy.wcs import WCS
from astropy.convolution.kernels import Tophat2DKernel
from astropy.convolution import convolve_fft

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from .dataprocessing import DustPediaDataProcessing
from .sample import DustPediaSample
from ...core.tools import filesystem as fs
from ...core.tools import time
from ...core.tools import tables, network
from ...magic.basics.coordinate import SkyCoordinate
from ...magic.tools import mosaicing
from ...core.tools.parallelization import ParallelTarget
from ...magic.misc import chrisfuncs
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.frame import Frame, sum_frames, sum_frames_quadratically

# -----------------------------------------------------------------

galex_bands = ["NUV", "FUV"]

# -----------------------------------------------------------------

class GALEXMosaicMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(GALEXMosaicMaker, self).__init__(config)

        # The DustPedia sample object
        self.sample = None

        # The DustPedia data processing instance
        self.dpdp = None

        # The NGC name of the galaxy
        self.ngc_name = None

        # The cutout center and width
        self.cutout_center = None
        self.cutout_width = None

        # Temporary root path
        self.root_path = None

        # Paths of subdirectories
        self.band_paths = dict()
        self.download_paths = dict()
        self.download_observations_paths = dict()
        self.download_response_paths = dict()
        self.download_background_paths = dict()
        self.download_counts_paths = dict()
        self.response_paths = dict()
        self.background_paths = dict()
        self.counts_paths = dict()

        # The URLs from where the images were downloaded
        self.observation_urls = dict()
        self.response_urls = dict()
        self.background_urls = dict()
        self.counts_urls = dict()

        # The rebin header and WCS
        self.rebin_header = None
        self.rebin_wcs = None

        # The mosaic frames and error maps
        self.mosaics = dict()
        self.mosaics_swarp = dict()
        self.error_maps = dict()
        self.relative_error_maps = dict()

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

        # 4. Get the data
        self.get_data()

        # 5. Filter
        self.filter()

        # 6. Clean
        self.clean()

        # 7. Level
        self.level()

        # 8. Get exposure times
        self.get_exposure_times()

        # 9. Mosaic
        self.mosaic()

        # 10. Finish
        self.finish()

        # 11. Show
        self.show()

        # 12. Write
        self.write()

    # -----------------------------------------------------------------

    def nurls_for_band(self, band):

        """
        This function ...
        :param band:
        :return:
        """

        return len(self.observation_urls[band])

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GALEXMosaicMaker, self).setup()

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
        self.root_path = fs.join(fs.home(), time.unique_name("GALEX_" + self.config.galaxy_name))
        fs.create_directory(self.root_path)

        # Loop over the bands
        for band in self.config.bands:

            # Create directory for this band
            path = fs.create_directory_in(self.root_path, band)

            # Create the download path
            download_path = fs.create_directory_in(path, "download")

            # Create subdirectories of the download directory
            download_observations_path = fs.create_directory_in(download_path, "observations")
            download_response_path = fs.create_directory_in(download_path, "response")
            download_background_path = fs.create_directory_in(download_path, "background")
            download_counts_path = fs.create_directory_in(download_path, "counts")

            # Create directories
            response_path = fs.create_directory_in(path, "response")
            background_path = fs.create_directory_in(path, "background")
            counts_path = fs.create_directory_in(path, "counts")

            # Set paths
            self.download_paths[band] = download_path
            self.download_observations_paths[band] = download_observations_path
            self.download_response_paths[band] = download_response_path
            self.download_background_paths[band] = download_background_path
            self.download_counts_paths[band] = download_counts_path
            self.response_paths[band] = response_path
            self.background_paths[band] = background_path
            self.counts_paths[band] = counts_path

    # -----------------------------------------------------------------

    def get_range(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting cutout range ...")

        # Get coordinate range for target image
        ra, dec, width = self.dpdp.get_cutout_range_for_galaxy(self.ngc_name)

        # Set the cutout center and cutout width
        self.cutout_center = SkyCoordinate(ra, dec, unit="deg", frame="fk5")
        self.cutout_width = width

    # -----------------------------------------------------------------

    def get_target_header(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the header and WCS for the mosaic ...")

        # Get the target header
        self.rebin_header = self.dpdp.get_header_for_galaxy(self.ngc_name, "GALEX")

        # Create WCS
        self.rebin_wcs = CoordinateSystem(self.rebin_header)

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the data ...")

        # 1. Get the URLs
        self.get_urls()

        # 2. Download
        self.download()

    # -----------------------------------------------------------------

    def get_urls(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the observation URLs ...")

        # Get the urls
        nuv_urls, fuv_urls = self.get_galex_observation_urls_for_galaxy()
        urls = dict()
        urls["NUV"] = nuv_urls
        urls["FUV"] = fuv_urls

        # Loop over the bands
        for band in self.config.bands:

            # Set observation urls
            self.observation_urls[band] = urls[band]

            # Determine the URLS for the response maps
            self.response_urls[band] = [url.replace("-int.fits.gz", "-rr.fits.gz") for url in urls[band]]

            # Determine the URLs for the background maps
            self.background_urls[band] = [url.replace("-int.fits.gz", "-skybg.fits.gz") for url in urls[band]]

            # Determine the URLS for the maps in original detector counts
            self.counts_urls[band] = [url.replace("-int.fits.gz", "-cnt.fits.gz") for url in urls[band]]

            # Debugging
            log.debug("Number of observations for the " + band + " band: " + str(self.nurls_for_band(band)))

    # -----------------------------------------------------------------

    def get_galex_observation_urls_for_galaxy(self):

        """
        This function ...
        :return:
        """

        # Get the tilenames
        tilenames = self.get_galex_tilenames_for_galaxy()

        urls_nuv = []
        urls_fuv = []

        # Search through the URL table to get all the URLS that contain one of the tilenames
        for i in range(len(self.dpdp.galex_url_table)):

            # Get the url
            url = self.dpdp.galex_url_table["URL"][i]

            # Check whether NUV of FUV observation
            nuv_or_fuv = "nuv" if fs.strip_extension(fs.name(url), double=True).split("-int")[0].split("-")[1] == "nd" else "fuv"

            for tilename in tilenames:

                if tilename in url:
                    if nuv_or_fuv == "nuv": urls_nuv.append(url)
                    else: urls_fuv.append(url)
                    break  # URL added, so no need to look at the other tilenames for this URL

        # Return the list of URLS
        return urls_nuv, urls_fuv

    # -----------------------------------------------------------------

    def get_galex_tilenames_for_galaxy(self):

        """
        This function ...
        :return:
        """

        # Find the indices in the table
        indices = tables.find_indices(self.dpdp.galex_observations_table, self.ngc_name, column_name="uploadID")

        # Initialize set of tilenames
        tilenames = set()

        # Loop over the indices corresponding to the specified galaxy
        for index in indices:

            tilename = self.dpdp.galex_observations_table["tilename"][index]
            tilenames.add(tilename)

        # Return the tilenames as a list
        return list(tilenames)

    # -----------------------------------------------------------------

    def download(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading the data (in parallel) ...")

        # Loop over the bands
        for band in self.config.bands:

            # Debugging
            log.debug("Downloading the data for the " + band + " band ...")

            # Parallel execution
            with ParallelTarget(network.download_files, self.config.nprocesses) as target:

                # Call the target function
                target(self.observation_urls[band], self.download_observations_paths[band])
                target(self.response_urls[band], self.download_response_paths[band])
                target(self.background_urls[band], self.download_background_paths[band])
                target(self.counts_urls[band], self.download_counts_paths[band])

    # -----------------------------------------------------------------

    def filter(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Filtering the observations on different criteria ...")

        # 1. Filter non overlapping
        self.filter_non_overlapping()

        # 2. Filter observations based on the coverage
        self.filter_coverage()

    # -----------------------------------------------------------------

    def filter_non_overlapping(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Filtering non-overlapping observations ...")

        # Parallel execution
        with ParallelTarget(mosaicing.filter_non_overlapping, self.config.nprocesses) as target:

            # Loop over the bands
            for band in self.config.bands: target(self.ngc_name, band, self.download_observations_paths[band], self.cutout_center, self.cutout_width, mode='point')

    # -----------------------------------------------------------------

    def filter_coverage(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Filtering observations based on their (flux) coverage ...")

        # Parallel execution
        with ParallelTarget(filter_galex_tiles, self.config.nprocesses) as target:

            # Loop over the bands
            # galaxy_name, tiles_path, ra, dec, width_deg, temp_raw_path, band_dict
            for band in self.config.bands: target(self.ngc_name, self.download_observations_paths[band], self.cutout_center, self.cutout_width, band)

    # -----------------------------------------------------------------

    def clean(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cleaning GALEX tiles ...")

        # CLEAN GALEX TILES
        #clean_galex_tiles(raw_files, id_string, nprocesses, working_path, temp_path_band, temp_reproject_path, band_dict)

        # Loop over raw tiles, creating exposure maps, and cleaning images to remove null pixels (also, creating convolved maps for later background fitting)
        #log.info('Cleaning ' + str(len(raw_files)) + ' raw maps for ' + id_string + " ...")

        # CLEAN
        #for raw_file in raw_files: clean_galex_tile(raw_file, working_path, temp_path_band, temp_reproject_path,
        #                                            band_dict)

        # Loop over the bands
        for band in self.config.bands:

            # Debugging
            log.debug("Cleaning GALEX " + band + " tiles ...")

            # Parallel execution
            with ParallelTarget(clean_galex_tile, self.config.nprocesses) as target:

                # raw_file, working_path, temp_path_band, temp_reproject_path, band_dict
                #for raw_file in raw_files: target(raw_file, working_path, temp_path_band, temp_reproject_path, band_dict)
                for path in fs.files_in_path(self.download_observations_paths[band]): target(path)

    # -----------------------------------------------------------------

    def level(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Peforming background-leveling on the observations ...")

        # Count image files, and move to reprojection directory
        mosaic_count = 0
        for filename in os.listdir(temp_raw_path):
            if '.fits' in filename:
                mosaic_count += 1
                new_path = fs.join(temp_reproject_path, filename)  ## NEW
                if fs.is_file(new_path): fs.remove_file(new_path)  ## NEW
                fs.move_file(filename, temp_reproject_path)

        # If more than one image file, commence background-matching
        if mosaic_count > 1:

            # Inform the user
            log.info("Matching background of " + id_string + " maps ...")

            # ...
            level_galex_maps(temp_reproject_path, temp_convolve_path, 'int.fits')

    # -----------------------------------------------------------------

    def get_exposure_times(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting exposure times ...")

        # Create a dictionary for the exposure times for each image
        exposure_times = dict()

        # Determine how the files are named
        filename_ends = "-" + band_dict['band_short'] + "-int"  # -fd-int for FUV, -nd-int for NUV
        filename_ends_no_int = "-" + band_dict['band_short']  # -fd for FUV, -nd for NUV

        # Get the exposure time for each image
        for path, name in fs.files_in_path(temp_reproject_path, extension="fits", contains=filename_ends,
                                           returns=["path", "name"]):
            # Open the header
            header = getheader(path)

            # Search for the exposure time
            exp_time = get_total_exposure_time(header)

            # Set the exposure time
            image_name = name.split(filename_ends)[0]
            exposure_times[image_name] = exp_time

    # -----------------------------------------------------------------

    def reproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Reprojecting ...")

        # Reproject image and weight prior to coaddition
        montage.commands.mImgtbl(temp_reproject_path, metatable_path, corners=True)
        proj_stats_path = fs.join(temp_path_band, id_string + "_Proj_Stats.txt")
        montage.commands.mProjExec(metatable_path, header_path, temp_swarp_path, proj_stats_path,
                                   raw_dir=temp_reproject_path, debug=False, exact=True, whole=False)
        # WHOLE IS IMPORTANT HERE
        # WE ACTUALLY DON'T WANT TO REPROJECT TO THE EXACT PIXELGRID DEFINED BY THE HEADER HERE,
        # BUT RATHER CUT OFF THE ORIGINAL MAPS WHERE THEY ARE OUTSIDE OF THE FIELD OF VIEW OF THE HEADER DEFINED AREA
        # SO NO ADDING OF NEW PIXELS IS DONE TO HAVE THE EXACT PIXELGRID DEFINED BY THE HEADER
        # THIS IS PRESUMABLY DONE JUST TO MAKE THE SWARPING MORE EFFICIENT ??

        # Rename reprojected files for SWarp
        for listfile in os.listdir(temp_swarp_path):
            if '_area.fits' in listfile:
                os.remove(fs.join(temp_swarp_path, listfile))
            elif 'hdu0_' in listfile:
                os.rename(fs.join(temp_swarp_path, listfile), fs.join(temp_swarp_path, listfile.replace('hdu0_', '')))

    # -----------------------------------------------------------------

    def mosaic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making GALEX mosaic for " + self.ngc_name + " and map of relative poisson errors ...")

        # Mosaic using Swarp
        self.mosaic_swarp()

        # Mosaic PTS method
        self.mosaic_pts()

    # -----------------------------------------------------------------

    def mosaic_swarp(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making mosaics with SWARP ...")

        # Parallel execution
        output = dict()
        with ParallelTarget(mosaic_with_swarp, self.config.nprocesses) as target:

            # Loop over the bands
            for band in self.config.bands: output[band] = mosaic_with_swarp(id_string, width_deg, pix_size, temp_swarp_path, ra_deg, dec_deg)

        # Loop over bands, LOAD SWARP RESULT, PRESUMABLY IN COUNTS/S
        for band in self.config.bands:

            swarp_result_path = output[band][0]

            # Load the resulting frame
            out_image = Frame.from_file(swarp_result_path)
            out_image.unit = "count/s"

            out_image[out_image == 0] = np.NaN
            out_image[out_image < -1E3] = np.NaN
            out_image[out_image <= 1E-8] = 0

            # Set the mosaic
            self.mosaics_swarp[band] = out_image

            # Write the swarped image
            #swarp_output_path = fs.join(output_path, id_string + "_swarp.fits")
            #out_image.saveto(swarp_output_path)

    # -----------------------------------------------------------------

    def mosaic_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making mosaic with PTS ...")

        # Loop over the bands
        for band in self.config.bands:

            # Get list of images
            self.get_images_for_mosaic(band)

            # Create directories
            temp_noise_path, temp_converted_path, temp_rebinned_path, temp_mosaic_path = self.create_directories_for_mosaic(band)

            # The path to the directory with all the tiles in counts for the current band (FUV or NUV)
            counts_path_band = fs.join(working_path, "counts", band_dict["band_long"])

        self.make_noise()

        self.convert_to_count_s_sr()

        self.rebin_for_mosaic()

    # -----------------------------------------------------------------

    def get_images_for_mosaic(self, band):

        """
        This function ...
        :return:
        """

        # The list of image names to be used for the mosaic
        image_names_for_mosaic = []

        # Loop over the files in the temp_swarp_path directory, where the tiles from temp_reproject_path are _peakand saved to
        for filename in os.listdir(temp_swarp_path):

            # SKIP WEIGHT FILES
            if not filename.endswith("-int.fits"): continue

            # Get the image name
            image_name = filename.split(filename_ends)[0]

            # Add the image name to the list
            image_names_for_mosaic.append(image_name)

    # -----------------------------------------------------------------

    def create_directories_for_mosaic(self, band):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating temporary directories ...")

        # Create a directory for the noise maps in counts per second
        temp_noise_path = fs.create_directory_in(temp_path_band, "noise")

        # Temp directory for the images and poisson noise maps in count/s/sr
        temp_converted_path = fs.create_directory_in(temp_path_band, "converted")

        # REBINNED in counts per second PER SR
        temp_rebinned_path = fs.create_directory_in(temp_path_band, "rebinned")

        # MOSAICING
        temp_mosaic_path = fs.create_directory_in(temp_path_band, "mosaic")

        # RESULT AFTER MOSAICING AND BACK TO COUNTS/S
        # temp_result_path = fs.create_directory_in(temp_path_band, "result")

        return temp_noise_path, temp_converted_path, temp_rebinned_path, temp_mosaic_path

    # -----------------------------------------------------------------

    def make_noise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making noise maps ...")

        # Make noise maps in count/s
        make_noise_maps_in_cps(band_dict, image_names_for_mosaic, counts_path_band, temp_noise_path, exposure_times)

    # -----------------------------------------------------------------

    def convert_to_count_s_sr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting to counts/s/sr ...")

        # Convert to counts/s/sr
        convert_frames_and_error_maps_to_per_solid_angle(band_dict, image_names_for_mosaic, temp_reproject_path,
                                                         temp_noise_path, temp_converted_path)

    # -----------------------------------------------------------------

    def rebin_for_mosaic(self):

        """
        This function ...
        :return:
        """

        # Rebin frames and error maps
        rebin_frames_and_error_maps(temp_converted_path, temp_rebinned_path, header_path, metatable_path,
                                    proj_stats_path)

        # Rebin the weight maps
        rebin_weight_maps(band_dict, image_names_for_mosaic, temp_reproject_path, temp_rebinned_path, header_path)

        for band in self.config.bands:

            with ParallelTarget() as target:



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

            # SOME THINGS
            mosaic[mosaic.data == 0] = np.NaN
            mosaic[mosaic.data < -1E3] = np.NaN
            mosaic[mosaic.data <= 1E-8] = 0

            # CALCULATE RELATIVE POISSON ERROR MAP
            relerrors = errors / mosaic
            relerrors[relerrors < 0.] = 0.0  # set negative values for relative error map to zero
            relerrors.replace_nans(0.0)  # set NaN values (because mosaic was zero) to zero

            # Set
            self.relative_error_maps[band] = relerrors

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

        # Write the mosaics
        self.write_mosaics()

        # Write the error maps
        self.write_error_maps()

        # Write the relative error maps
        self.write_relative_error_maps()

    # -----------------------------------------------------------------

    def write_mosaics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mosaics ...")

        # Loop over the bands
        for band in self.config.bands:

            # Determine path
            id_string = self.ngc_name + "_SDSS_" + band
            path = fs.join(output_path, id_string + ".fits")

            # Save mosaic frame as FITS file
            self.mosaics[band].saveto(path)

    # -----------------------------------------------------------------

    def write_error_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the error maps ...")

        # Loop over the bands
        for band in self.config.bands:

            # Determine path
            id_string = self.ngc_name + "_SDSS_" + band
            path = fs.join(output_path, id_string + "_errors.fits")

            # Save error map as FITS file
            self.error_maps[band].saveto(path)

    # -----------------------------------------------------------------

    def write_relative_error_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the relative error maps ...")

        # Loop over the bands
        for band in self.config.bands:

            # Determine path
            id_string = self.ngc_name + "_SDSS_" + band
            path = fs.join(output_path, id_string + "_relerrors.fits")

            # Save relative error map as FITS file
            self.relative_error_maps[band].saveto(path)

# -----------------------------------------------------------------

def nuv_or_fuv(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Get header
    header = getheader(path)

    # Check band
    band = header["BAND"]

    # Return
    if band == 1: return "NUV"
    elif band == 2: return "FUV"
    else: raise RuntimeError("Invalid band: " + str(band))

# -----------------------------------------------------------------

def filter_galex_tiles(galaxy_name, tiles_path, center, width, band):

    """
    This function ...
    :return:
    """

    meta_path = fs.join(tiles_path, "meta.dat")
    new_overlap_path = fs.join(tiles_path, "overlap_circle.dat")

    # Get overlapping file paths
    ra = center.ra.to("deg").value
    dec = center.dec.to("deg").value
    width = width.to("deg").value
    overlapping_file_paths = mosaicing.generate_overlapping_file_paths(meta_path, ra, dec, meta_path, mode="circle", radius=(0.5 * width) * (2.0 ** 0.5))

    # Check
    #if len(overlapping_file_paths.shape) == 0: overlapping_file_paths = [overlapping_file_paths.tolist()]
    #for overlapping_file_path in overlapping_file_paths: fs.copy_file(overlapping_file_path, temp_raw_path)

    # NEW: FILTER IN-PLACE
    for path in fs.files_in_path(tiles_path):
        if path in overlapping_file_paths: continue
        else: fs.remove_file(path)

    temp_raw_path = tiles_path

    # Uncompress .fits.gz files
    # [os.system('gunzip '+ listfile) for listfile in os.listdir(raw_in_temp_dir)]

    # Ensure that at least one of the raw GALEX tiles has actual flux coverage at location of source
    # raw_files = os.listdir(temp_raw_path) ## THERE WAS A 'deg' FILE IN THE DIRECTORY AS WELL WHICH WAS NOT A FITS FILE BUT A PLAIN TEXT HEADER FILE, AND SO THIS WASN'T WORKING
    raw_files = fs.files_in_path(temp_raw_path, extension="fits", returns="name")
    raw_files = [name + ".fits" for name in raw_files]
    coverage = False
    for raw_file in raw_files:

        # Read in map
        in_fitsdata = open_fits(fs.join(temp_raw_path, raw_file))
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_fitsdata.close()

        # Locate pixel coords
        in_wcs = WCS(in_header)
        location_pix = in_wcs.wcs_world2pix(np.array([[np.float(ra), np.float(dec)]]), 0)[0]
        pix_i, pix_j = location_pix[1], location_pix[0]

        # Evalulate coverage at location, and proceed accordingly
        if True in [coord <= 0 for coord in [pix_i - 10, pix_i + 11, pix_j - 10, pix_j + 11]]:
            continue
        try:
            image_slice = in_image[pix_i - 10:pix_i + 11, pix_j - 10:pix_j + 11]
        except: continue
        if np.where(image_slice > 0)[0].shape[0] > 0:
            coverage = True

    # No coverage
    if not coverage: raise RuntimeError('No GALEX ' + band + ' coverage for ' + galaxy_name)

    # Return ...
    #return raw_files

# -----------------------------------------------------------------

def clean_galex_tile(raw_file_path, working_path, temp_path_band, temp_reproject_path, band_dict):

    """
    Function to clean GALEX tiles and create exposure maps
    :param path: raw file path
    :param working_path:
    :param temp_path_band:
    :param temp_reproject_path:
    :param band_dict:
    :return:
    """

    raw_file_name = fs.name(raw_file_path)

    # Inform the user ...
    log.info("Cleaning map " + raw_file_path + " ...")
    #log.info(" - raw_file = " + raw_file)
    #log.info(" - working_path = " + working_path)
    #log.info(" - temp_path_band = " + temp_path_band)
    #log.info(" - temp_reproject_path = " + temp_reproject_path)
    #log.info(" - band_dict = " + str(band_dict))

    # Response and background paths for this band
    #response_path = fs.join(working_path, "response", band_dict['band_long'])
    #background_path = fs.join(working_path, "background", band_dict['band_long'])

    temp_raw_path = fs.join(temp_path_band, "raw")

    # Read in image
    #raw_file_path = fs.join(temp_raw_path, raw_file)
    in_fitsdata = open_fits(raw_file_path)
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    out_image = in_image.copy()

    # Load and align response map
    rr_path = fs.join(response_path, raw_file_name.replace('-int.fits','-rr.fits'))
    rr_fitsdata = open_fits(rr_path)
    rr_image = rr_fitsdata[0].data
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean image using response map
    out_image[ np.where( rr_image <= 1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_path = fs.join(background_path, raw_file_name.replace('-int.fits','-skybg.fits'))
    bg_fitsdata = open_fits(bg_path)
    bg_image = bg_fitsdata[0].data
    bg_zoom = np.float(out_image.shape[0]) / np.float(bg_image.shape[0])
    bg_image = interpolation.zoom(bg_image, bg_zoom, order=0)

    # Clean image using sky background map
    out_image[ np.where( bg_image <= 1E-10 ) ] = np.NaN

    """
    # Load and align flag map
    flag_path = root_dir+'Flags/'+band_dict['band_long']+'/'+raw_file.replace('-int.fits','-flags.fits.gz')
    flag_fitsdata = fits.open(flag_path)
    flag_image = flag_fitsdata[0].data
    flag_zoom = np.float(out_image.shape[0]) / np.float(flag_image.shape[0])
    flag_image = scipy.ndimage.interpolation.zoom(flag_image, flag_zoom, order=0)

    # Nullify pixels where the bitwise flags indicate dichoric reflections, window refections, variable-pixel masks or hot-pixel masks
    out_image[ np.where( flag_image.astype(int) & (1<<1) > 0 ) ] = np.NaN # Dichoric reflection
    out_image[ np.where( flag_image.astype(int) & (1<<2) > 0 ) ] = np.NaN # Window reflection
    out_image[ np.where( flag_image.astype(int) & (1<<7) > 0 ) ] = np.NaN # Variable-pixel mask
    out_image[ np.where( flag_image.astype(int) & (1<<8) > 0 ) ] = np.NaN # Hot-pixel mask
    """

    # Set all remaining, and hence "true", zero pixels to be ever-so-slighly non-zero
    out_image += 1E-8

    # Find centre of coverage area
    cov_i = ((np.where( np.isnan(out_image)==False ))[0])
    cov_j = ((np.where( np.isnan(out_image)==False ))[1])
    cov_ellipse = chrisfuncs.EllipseFit(cov_i, cov_j)
    cov_centre = cov_ellipse[0]
    cov_centre_i, cov_centre_j = cov_centre[0], cov_centre[1]

    # Set all pixels more than 35 arcmin (1400 pizels) from centre to be NaN, as these are typically low-quality
    cov_trim_mask = chrisfuncs.EllipseMask(out_image, 1400, 1.0, 0.0, cov_centre_i, cov_centre_j)
    out_image[ np.where(cov_trim_mask==0) ] = np.NaN

    # Save cleaned image
    out_hdu = PrimaryHDU(data=out_image, header=in_header)
    out_hdulist = HDUList([out_hdu])
    out_hdulist.writeto(fs.join(temp_raw_path, raw_file_name), clobber=True)

    # Create convolved version of map, for later use in background-matching
    """
    if np.isnan(out_image).sum()==0:
        conv_image = scipy.ndimage.filters.gaussian_filter(out_image, 20)
    else:
    """

    temp_convolve_path = fs.join(temp_path_band, "convolve")

    kernel = Tophat2DKernel(10)
    conv_image = convolve_fft(out_image, kernel, interpolate_nan=False, normalize_kernel=True, ignore_edge_zeros=False, allow_huge=True) #, interpolate_nan=True, normalize_kernel=True)

    # Write
    temp_convolve_image_path = fs.join(temp_convolve_path, raw_file_name)                  ## NEW
    if fs.is_file(temp_convolve_image_path): fs.remove_file(temp_convolve_image_path) ## NEW
    writeto(temp_convolve_image_path, conv_image, in_header)

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[ np.where( np.isnan(out_image)==False ) ] = (float(in_header['exptime']))**0.5 # SQUARE ROOT OF THE EXPOSURE TIME
    exp_hdu = PrimaryHDU(data=exp_image, header=in_header)
    exp_hdulist = HDUList([exp_hdu])

    # Write
    temp_reproject_image_path = fs.join(temp_reproject_path, raw_file_name.replace('.fits','.wgt.fits'))  ## NEW
    if fs.is_file(temp_reproject_image_path): fs.remove_file(temp_reproject_image_path)              ## NEW
    exp_hdulist.writeto(temp_reproject_image_path)

# -----------------------------------------------------------------

def level_chi_squared(level_params, image):

    """
    Fit flat plane to the image to find level
    :param level_params:
    :param image:
    :return:
    """

    level = level_params['level'].value
    chi = image - level
    chisq = chi**2.0

    return chisq

# -----------------------------------------------------------------

def level_galex_maps(fitsfile_dir, convfile_dir, target_suffix):

    """
    original name: GALEX_Zero
    Set a set of maps to the same level
    :param fitsfile_dir:
    :param convfile_dir:
    :param target_suffix:
    :return:
    """

    # Inform the user
    log.info("Inside the 'level_galex_maps' function with:")
    log.info(" - fitsfile_dir = " + fitsfile_dir)
    log.info(" - convfile_dir = " + convfile_dir)
    log.info(" - target_suffix = " + target_suffix)

    # Make list of files in target directory that have target suffix
    allfile_list = os.listdir(fitsfile_dir)
    fitsfile_list = []
    for allfile in allfile_list:
        if target_suffix in allfile:
            fitsfile_list.append(allfile)

    # Loop over each file
    for i in range(0, len(fitsfile_list)):

        # Inform the user
        log.info('Matching background of map ' + fitsfile_list[i] + " ...")

        # Read in corresponding map from directory containing convolved images
        fitsdata_conv = open_fits(convfile_dir + '/' + fitsfile_list[i])
        image_conv = fitsdata_conv[0].data
        fitsdata_conv.close()

        # Fit to level of image; save if first image, otherwise calculate appropriate offset
        level_params = lmfit.Parameters()

        #level_params.add('level', value=np.nanmedian(image_conv), vary=True) ## FOR NUMPY VERSION 1.9.0 AND ABOVE
        image_conv_nonans = image_conv[np.logical_not(np.isnan(image_conv))]
        level_params.add('level', value=np.median(image_conv_nonans), vary=True) # BELOW NUMPY VERSION 1.9.0
        #

        image_conv_clipped = chrisfuncs.SigmaClip(image_conv, tolerance=0.005, median=False, sigma_thresh=3.0)[2]
        level_result = lmfit.minimize(level_chi_squared, level_params, args=(image_conv_clipped.flatten(),))
        level = level_result.params['level'].value
        if i==0:
            level_ref = level
            continue
        average_offset = level_ref - level
        #print 'Applying offset of '+str(average_offset)+' to '+fitsfile_list[i]

        """
        # Save floor and peak values
        floor_value = np.nanmin(image_conv)
        peak_value = chrisfuncs.SigmaClip( image_conv, tolerance=0.00025, median=False, sigma_thresh=3.0)[1]
        floor_value_list.append(floor_value)
        peak_value_list.append(peak_value)
        if i==0:
            floor_value_ref = floor_value
            peak_value_ref = peak_value
            continue

        # Calculate offsets
        floor_offset = floor_value_ref - floor_value
        peak_offset = peak_value_ref - peak_value
        average_offset = peak_offset#np.mean([ floor_offset, peak_offset ])
        """

        # Read in unconvolved file, and apply offset
        fitsdata_in = open_fits(fitsfile_dir+'/'+fitsfile_list[i])
        image_in = fitsdata_in[0].data
        header_in = fitsdata_in[0].header
        fitsdata_in.close()
        image_out = image_in + average_offset
        #print 'Map mean of '+fitsfile_list[i]+' changed from '+str(np.nanmean(image_in))+' to '+str(np.nanmean(image_out))

        # Save corrected file
        image_out_hdu = PrimaryHDU(data=image_out, header=header_in)
        image_out_hdulist = HDUList([image_out_hdu])
        image_out_hdulist.writeto(fitsfile_dir+'/'+fitsfile_list[i], clobber=True)

# -----------------------------------------------------------------

def make_noise_maps_in_cps(band_dict, image_names_for_mosaic, counts_path_band, temp_noise_path, exposure_times):

    """
    This function ...
    :param band_dict:
    :param image_names_for_mosaic:
    :param counts_path_band:
    :param temp_noise_path:
    :param exposure_times:
    :return:
    """

    # Inform the user
    log.info("Creating maps of the poisson noise for each GALEX tile in counts per second ...")

    # Loop over the file names
    for image_name in image_names_for_mosaic:

        # Get the path to the file in counts
        counts_filepath = fs.join(counts_path_band, image_name + "-" + band_dict["band_short"] + "-cnt.fits")

        # Load the counts map
        counts_frame = Frame.from_file(counts_filepath)
        counts_frame.unit = "ct"

        # Calculate the poisson frame
        poisson = Frame(np.sqrt(counts_frame.data))   # calculate the poisson error (in counts) in every pixel
        poisson.wcs = counts_frame.wcs
        poisson.unit = "ct"

        # Get the exposure time in seconds
        exposure_time = exposure_times[image_name]

        # Convert the poisson frame to counts/second
        poisson /= exposure_time
        poisson.unit = "ct/s"

        # Determine the path to the poisson noise map in counts per second
        new_path = fs.join(temp_noise_path, image_name + ".fits")

        # Save the poisson noise map in counts per second
        poisson.saveto(new_path)

# -----------------------------------------------------------------

def convert_frames_and_error_maps_to_per_solid_angle(band_dict, image_names_for_mosaic, temp_reproject_path,
                                                     temp_noise_path, temp_converted_path):

    """
    This function ...
    :param band_dict:
    :param image_names_for_mosaic:
    :param temp_reproject_path:
    :param temp_noise_path:
    :param temp_converted_path:
    :return:
    """

    # Inform the user
    log.info("Converting the units to luminosity per solid angle ...")

    # Determine how the files are named
    filename_ends = "-" + band_dict['band_short'] + "-int"  # -fd-int for FUV, -nd-int for NUV
    filename_ends_no_int = "-" + band_dict['band_short']    # -fd for FUV, -nd for NUV

    # DO THE UNIT CONVERSION
    for image_name in image_names_for_mosaic:

        # Debugging
        log.debug("Converting " + image_name + " frame and error map to count / s / sr...")

        # Determine the full path to the tile in temp_reproject_path
        filepath = fs.join(temp_reproject_path, image_name + filename_ends + ".fits")

        # Load the frame
        frame = Frame.from_file(filepath)
        frame.unit = "count/s"

        # NUMBER OF SR PER PIXEL
        pixelsr = frame.pixelarea.to("sr").value

        # CONVERT FRAME
        frame /= pixelsr  # IN COUNTS PER SECOND PER SR NOW
        frame.unit = "count/(s*sr)" # set new unit

        # Determine the full path to the poisson noise map for this tile
        poisson_filepath = fs.join(temp_noise_path, image_name + ".fits")

        # Load the poisson noise frame
        poisson = Frame.from_file(poisson_filepath)
        # the unit should be set already to count/s

        # normally, this shouldn't be necessary
        poisson_pixelsr = poisson.pixelarea.to("sr").value
        assert np.isclose(pixelsr, poisson_pixelsr)

        # CONVERT ERROR MAP TO COUNTS PER SECOND PER SR
        poisson /= pixelsr  # IN COUNTS PER SECOND PER SR NOW
        poisson.unit = "count/(s*sr)" # set new unit

        # SAVE THE TILE IN COUNTS/S/SR
        frame_path = fs.join(temp_converted_path, image_name + ".fits")
        frame.saveto(frame_path)

        # SAVE THE ERROR MAP OF THE TILE IN COUNTS/S/SR
        errors_path = fs.join(temp_converted_path, image_name + "_error.fits")
        poisson.saveto(errors_path)

# -----------------------------------------------------------------

def rebin_frames_and_error_maps(temp_converted_path, temp_rebinned_path, header_path, metatable_path, proj_stats_path):

    """
    This function ...
    :param temp_converted_path:
    :param temp_rebinned_path:
    :param header_path:
    :param metatable_path:
    :param proj_stats_path:
    :return:
    """

    # Inform the user
    log.info("Rebinning the frames and poisson noise maps in counts / second to the target coordinate system ...")

    # Rebin the batch of images (tiles and noise maps) in counts/s/sr
    #rebin_batch_with_montage(temp_converted_path, temp_rebinned_path, metatable_path, header_path, proj_stats_path)

    # GET REBIN WCS
    rebin_header = Header.fromtextfile(header_path)
    rebin_wcs = CoordinateSystem(rebin_header)

    # Loop over all FITS files
    for path, name in fs.files_in_path(temp_converted_path, extension="fits", returns=["path", "name"]):

        # Load the image
        frame = Frame.from_file(path)

        # Rebin
        frame.rebin(rebin_wcs)

        # Determine the new path
        new_path = fs.join(temp_rebinned_path, name + ".fits")

        # Save
        frame.saveto(new_path)

# -----------------------------------------------------------------

def rebin_weight_maps(band_dict, image_names_for_mosaic, temp_reproject_path, temp_rebinned_path, header_path):

    """
    This function ...
    :param band_dict:
    :param image_names_for_mosaic:
    :param temp_reproject_path:
    :param temp_rebinned_path:
    :param header_path:
    :return:
    """

    # Inform the user
    log.info("Rebinning weight maps ...")

    # Determine how the files are named
    filename_ends = "-" + band_dict['band_short'] + "-int"  # -fd-int for FUV, -nd-int for NUV
    filename_ends_no_int = "-" + band_dict['band_short']  # -fd for FUV, -nd for NUV

    # GET REBIN WCS
    rebin_header = Header.fromtextfile(header_path)
    rebin_wcs = CoordinateSystem(rebin_header)

    # Loop over the image_names
    for image_name in image_names_for_mosaic:

        # Determine the path to the weight map
        weight_path = fs.join(temp_reproject_path, image_name + filename_ends + ".wgt.fits")

        # Load the weight map
        weights = Frame.from_file(weight_path)

        # Rebin
        weights.rebin(rebin_wcs)

        # Determine the new path
        new_path = fs.join(temp_rebinned_path, image_name + "_weight.fits")

        # Save
        weights.saveto(new_path)

# -----------------------------------------------------------------

def combine_frames_and_error_maps(image_names_for_mosaic, temp_rebinned_path, temp_mosaic_path, wcs):

    """
    This function ...
    :param image_names_for_mosaic:
    :param temp_rebinned_path:
    :param temp_mosaic_path:
    :param wcs:
    :return:
    """

    primary_frames = []
    error_frames = []
    weight_frames = []

    for image_name in image_names_for_mosaic:

        # Determine paths
        filepath = fs.join(temp_rebinned_path, image_name + ".fits")
        errorpath = fs.join(temp_rebinned_path, image_name + "_error.fits")
        weightpath = fs.join(temp_rebinned_path, image_name + "_weight.fits")

        # Load the frames
        frame = Frame.from_file(filepath)
        errors = Frame.from_file(errorpath)
        weights = Frame.from_file(weightpath)

        # Calculate the weighted frame
        frame_weighted = frame * weights
        errors_weighted = errors / weights

        # Convert the units from counts/s/sr to counts/s
        pixelsr = frame_weighted.pixelarea.to("sr").value
        frame_weighted *= pixelsr
        errors_weighted *= pixelsr

        # Create mask where the weights are nans
        mask = weights.nans()

        # Set zero
        frame_weighted[mask] = 0.0
        errors_weighted[mask] = 0.0
        weights[mask] = 0.0

        # Add to the list
        primary_frames.append(frame_weighted)
        error_frames.append(errors_weighted)
        weight_frames.append(weights)

    #  Calculate denominator for weighted average (mosaicing)
    normalization = sum_frames(*weight_frames)

    # CALCULATE THE MOSAIC FRAME IN COUNTS/S
    mosaic_frame = sum_frames(*primary_frames) / normalization
    mosaic_frame.wcs = wcs

    # CALCULATE THE MOSAIC ERROR MAP IN COUNTS/S
    mosaic_errormap = sum_frames_quadratically(*error_frames) / normalization
    mosaic_errormap.wcs = wcs

    ## DONE

    # SAVE THE MOSAIC IN COUNTS/S
    mosaic_path = fs.join(temp_mosaic_path, "mosaic.fits")
    mosaic_frame.saveto(mosaic_path)

    # SAVE THE MOSAIC ERROR MAP IN COUNTS/S
    mosaic_error_path = fs.join(temp_mosaic_path, "mosaic_errors.fits")
    mosaic_errormap.saveto(mosaic_error_path)

    # Return the mosaic and the mosaic error map in counts/s
    return mosaic_frame, mosaic_errormap

# -----------------------------------------------------------------

def mosaic_with_swarp(id_string, width_deg, pix_size, temp_swarp_path, ra_deg, dec_deg):

    """
    This function ...
    :return:
    """

    # Use SWarp to co-add images weighted by their error maps
    log.info("Co-adding " + id_string + " maps ...")

    # Determine image width in pixels
    image_width_pixels = str(int((float(width_deg) * 3600.) / pix_size))

    # Change directories
    os.chdir(temp_swarp_path)

    # EXECUTE SWARP
    swarp_command_string = 'swarp *int.fits -IMAGEOUT_NAME ' + id_string + '_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER ' + str(ra_deg) + ',' + str(dec_deg) + ' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE ' + image_width_pixels + ',' + image_width_pixels + ' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N  -RESAMPLE N -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT'
    os.system(swarp_command_string)

    # Swarp result path
    swarp_result_path = fs.join(temp_swarp_path, id_string + "_SWarp.fits")

    # Return the path to the resulting mosaic
    return swarp_result_path

# -----------------------------------------------------------------

def get_total_exposure_time(header):

    """
    This function ....
    :param header:
    :return:
    """

    # If the 'EXPTIME' keyword is defined in the header, return the value
    if "EXPTIME" in header: return float(header["EXPTIME"])

    # Otherwise, look at individual 'EXPTIME' keywords
    else:

        total = 0.0
        for key in header:

            # Add the values of all keys that state an exposure time (EXPT0001, EXPT0002, ...)
            if key.startswith("EXPT"): total += float(header[key])

        # Return the total exposure time
        return total

# -----------------------------------------------------------------

def rebin_batch_with_montage(in_path, out_path, metatable_path, header_path, proj_stats_path, remove_area=True):

    """
    This function ...
    :param in_path:
    :param out_path:
    :param metatable_path:
    :param header_path:
    :param proj_stats_path:
    :param remove_area:
    :return:
    """

    # Inform the user
    log.info("Rebinning a batch of images in " + in_path + " to " + out_path + " ...")

    # REPROJECT ALL INPUT MAPS TO TARGET HEADER
    montage.commands.mProjExec(metatable_path, header_path, out_path, proj_stats_path, raw_dir=in_path, debug=False, exact=True, whole=True)

    # Rename reprojected files
    for listfile in os.listdir(out_path):

        if '_area.fits' in listfile:
            if remove_area: os.remove(fs.join(out_path, listfile)) # remove area if requested
        elif 'hdu0_' in listfile:
            os.rename(fs.join(out_path, listfile), fs.join(out_path, listfile.replace('hdu0_', '')))

# -----------------------------------------------------------------

def mosaic_pts_for_band():

    """
    This function ...
    :return:
    """



    # DO THE COMBINING
    mosaic, errors = combine_frames_and_error_maps(image_names_for_mosaic, temp_rebinned_path, temp_mosaic_path,
                                                   self.rebin_wcs)

# -----------------------------------------------------------------
