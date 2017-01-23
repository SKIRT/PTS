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
import numpy as np

# Import astronomical modules
from astropy.io.fits import getheader
from astropy.io.fits import open as open_fits
from astropy.wcs import WCS

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from .dataprocessing import DustPediaDataProcessing
from .sample import DustPediaSample
from ...core.tools import filesystem as fs
from ...core.tools import time
from ...core.tools import tables, network
from ...magic.basics.coordinate import SkyCoordinate
from ...magic.tools import mosaic
from ...core.tools.parallelization import ParallelTarget

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

        # The mosaic frames and error maps
        self.mosaics = dict()
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

        # Create Montage FITS header
        location_string = str(ra_deg) + ' ' + str(dec_deg)
        pix_size = 3.2

        # Make header
        header_path = fs.join(temp_path_band, id_string + ".hdr")
        montage.commands.mHdr(location_string, width_deg, header_path, pix_size=pix_size)

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
        with ParallelTarget(mosaic.filter_non_overlapping, self.config.nprocesses) as target:

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

        # # Path to the overlap table
        #overlap_path = fs.join(temp_path_band, "overlap_table.dat")

        # FILTER GALEX tiles for this band
        #raw_files = filter_galex_tiles(galaxy_name, meta_path, overlap_path, ra, dec, width_deg, temp_raw_path, band_dict)

        # Parallel execution
        with ParallelTarget(filter_galex_tiles, self.config.nprocesses) as target:

            # Loop over the bands
            for band in self.config.bands: target(self.ngc_name, meta_path, overlap_path, ra, dec, width_deg, temp_raw_path, band_dict)

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

        # Parallel execution
        with ParallelTarget(clean_galex_tile, self.config.nprocesses):

            for raw_file in raw_files: target(raw_file, working_path, temp_path_band, temp_reproject_path, band_dict)

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
                shutil.move(filename, temp_reproject_path)

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

        # Create a dictionary for the exposure times for each image
        exposure_times = dict()

        # Determine how the files are named
        filename_ends = "-" + band_dict['band_short'] + "-int"  # -fd-int for FUV, -nd-int for NUV
        filename_ends_no_int = "-" + band_dict['band_short']  # -fd for FUV, -nd for NUV

        # Get the exposure time for each image
        for path, name in fs.files_in_path(temp_reproject_path, extension="fits", contains=filename_ends,
                                           returns=["path", "name"]):
            # Open the header
            header = fits.getheader(path)

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

        # MOSAIC WITH SWARP
        swarp_result_path = mosaic_with_swarp(id_string, width_deg, pix_size, temp_swarp_path, ra_deg, dec_deg)

        #### LOAD SWARP RESULT, PRESUMABLY IN COUNTS/S

        # Load the resulting frame
        out_image = Frame.from_file(swarp_result_path)
        out_image.unit = "count/s"

        out_image[out_image == 0] = np.NaN
        out_image[out_image < -1E3] = np.NaN
        out_image[out_image <= 1E-8] = 0

        # Write the swarped image
        swarp_output_path = fs.join(output_path, id_string + "_swarp.fits")
        out_image.saveto(swarp_output_path)

    # -----------------------------------------------------------------

    def mosaic_pts(self):

        """
        This function ...
        :return:
        """

        ####################### IMAGE NAMES

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

        ####################### DIRECTORIES

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

        ######################

        # The path to the directory with all the tiles in counts for the current band (FUV or NUV)
        counts_path_band = fs.join(working_path, "counts", band_dict["band_long"])

        # Make noise maps in count/s
        make_noise_maps_in_cps(band_dict, image_names_for_mosaic, counts_path_band, temp_noise_path, exposure_times)

        # Convert to counts/s/sr
        convert_frames_and_error_maps_to_per_solid_angle(band_dict, image_names_for_mosaic, temp_reproject_path,
                                                         temp_noise_path, temp_converted_path)

        # Rebin frames and error maps
        rebin_frames_and_error_maps(temp_converted_path, temp_rebinned_path, header_path, metatable_path,
                                    proj_stats_path)

        # Rebin the weight maps
        rebin_weight_maps(band_dict, image_names_for_mosaic, temp_reproject_path, temp_rebinned_path, header_path)

        # DO THE COMBINING
        # rebinned_path, footprints_path, mosaics_path, wcs
        wcs = CoordinateSystem(Header.fromtextfile(header_path))
        mosaic, errors = combine_frames_and_error_maps(image_names_for_mosaic, temp_rebinned_path, temp_mosaic_path,
                                                       wcs)

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return:
        """

        # SOME THINGS
        mosaic[mosaic.data == 0] = np.NaN
        mosaic[mosaic.data < -1E3] = np.NaN
        mosaic[mosaic.data <= 1E-8] = 0

        # CALCULATE RELATIVE POISSON ERROR MAP
        relerrors = errors / mosaic
        relerrors[relerrors < 0.] = 0.0  # set negative values for relative error map to zero
        relerrors.replace_nans(0.0)  # set NaN values (because mosaic was zero) to zero

        # Save mosaic as FITS file
        mosaic_output_path = fs.join(output_path, id_string + ".fits")
        mosaic.saveto(mosaic_output_path)

        # Save error map as FITS file
        errors_output_path = fs.join(output_path, id_string + "_errors.fits")
        errors.saveto(errors_output_path)

        # Save relative error map as FITS file
        relerrors_output_path = fs.join(output_path, id_string + "_relerrors.fits")
        relerrors.saveto(relerrors_output_path)

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

def filter_galex_tiles(galaxy_name, meta_path, overlap_path, ra, dec, width_deg, temp_raw_path, band_dict):

    """
    This fucntion ...
    :return:
    """

    # Use Montage image metadata table to identify and retrieve which raw GALEX tiles overlap with entire region of interest (handling the case of only a single file)
    montage.commands_extra.mCoverageCheck(meta_path, overlap_path, mode='circle', ra=ra.to("deg").value, dec=dec.to("deg").value, radius=(0.5 * width_deg) * (2.0 ** 0.5))

    # Get file paths of overlapping observations
    overlapping_file_paths = np.genfromtxt(overlap_path, skip_header=3, usecols=[32], dtype=str)

    if len(overlapping_file_paths.shape) == 0: overlapping_file_paths = [overlapping_file_paths.tolist()]
    for overlapping_file_path in overlapping_file_paths: fs.copy_file(overlapping_file_path, temp_raw_path)

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
        location_pix = in_wcs.wcs_world2pix(np.array([[np.float(ra.to("deg").value), np.float(dec.to("deg").value)]]), 0)[0]
        pix_i, pix_j = location_pix[1], location_pix[0]

        # Evalulate coverage at location, and proceed accordingly
        if True in [coord <= 0 for coord in [pix_i - 10, pix_i + 11, pix_j - 10, pix_j + 11]]:
            continue
        try:
            image_slice = in_image[pix_i - 10:pix_i + 11, pix_j - 10:pix_j + 11]
        except:
            continue
        if np.where(image_slice > 0)[0].shape[0] > 0:
            coverage = True

    # No coverage
    if not coverage: raise RuntimeError('No GALEX ' + band_dict['band_long'] + ' coverage for ' + galaxy_name)

    # Return ...
    return raw_files

# -----------------------------------------------------------------

def clean_galex_tile(raw_file, working_path, temp_path_band, temp_reproject_path, band_dict):

    """
    Function to clean GALEX tiles and create exposure maps
    :return:
    """

    # Inform the user ...
    log.info("Cleaning map " + raw_file + " ...")
    log.info(" - raw_file = " + raw_file)
    log.info(" - working_path = " + working_path)
    log.info(" - temp_path_band = " + temp_path_band)
    log.info(" - temp_reproject_path = " + temp_reproject_path)
    log.info(" - band_dict = " + str(band_dict))

    # Response and background paths for this band
    response_path = fs.join(working_path, "response", band_dict['band_long'])
    background_path = fs.join(working_path, "background", band_dict['band_long'])

    temp_raw_path = fs.join(temp_path_band, "raw")

    # Read in image
    in_fitsdata = fits.open(fs.join(temp_raw_path, raw_file))
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    out_image = in_image.copy()

    # Load and align response map
    rr_path = fs.join(response_path, raw_file.replace('-int.fits','-rr.fits'))
    rr_fitsdata = fits.open(rr_path)
    rr_image = rr_fitsdata[0].data
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = scipy.ndimage.interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean image using response map
    out_image[ np.where( rr_image <= 1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_path = fs.join(background_path, raw_file.replace('-int.fits','-skybg.fits'))
    bg_fitsdata = fits.open(bg_path)
    bg_image = bg_fitsdata[0].data
    bg_zoom = np.float(out_image.shape[0]) / np.float(bg_image.shape[0])
    bg_image = scipy.ndimage.interpolation.zoom(bg_image, bg_zoom, order=0)

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
    out_hdu = fits.PrimaryHDU(data=out_image, header=in_header)
    out_hdulist = fits.HDUList([out_hdu])
    out_hdulist.writeto(fs.join(temp_raw_path, raw_file), clobber=True)

    # Create convolved version of map, for later use in background-matching
    """
    if np.isnan(out_image).sum()==0:
        conv_image = scipy.ndimage.filters.gaussian_filter(out_image, 20)
    else:
    """

    temp_convolve_path = fs.join(temp_path_band, "convolve")

    kernel = astropy.convolution.kernels.Tophat2DKernel(10)
    conv_image = astropy.convolution.convolve_fft(out_image, kernel, interpolate_nan=False, normalize_kernel=True, ignore_edge_zeros=False, allow_huge=True) #, interpolate_nan=True, normalize_kernel=True)

    # Write
    temp_convolve_image_path = fs.join(temp_convolve_path, raw_file)                  ## NEW
    if fs.is_file(temp_convolve_image_path): fs.remove_file(temp_convolve_image_path) ## NEW
    fits.writeto(temp_convolve_image_path, conv_image, in_header)

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[ np.where( np.isnan(out_image)==False ) ] = (float(in_header['exptime']))**0.5 # SQUARE ROOT OF THE EXPOSURE TIME
    exp_hdu = fits.PrimaryHDU(data=exp_image, header=in_header)
    exp_hdulist = fits.HDUList([exp_hdu])

    # Write
    temp_reproject_image_path = fs.join(temp_reproject_path, raw_file.replace('.fits','.wgt.fits'))  ## NEW
    if fs.is_file(temp_reproject_image_path): fs.remove_file(temp_reproject_image_path)              ## NEW
    exp_hdulist.writeto(temp_reproject_image_path)

# -----------------------------------------------------------------
