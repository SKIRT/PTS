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
import traceback
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
from ...core.basics.log import log
from .dataprocessing import DustPediaDataProcessing
from .sample import DustPediaSample
from ...core.tools import filesystem as fs
from ...core.tools import tables, network
from ...magic.basics.coordinate import SkyCoordinate
from ...magic.tools import mosaicing
from ...core.tools.parallelization import ParallelTarget
from ...magic.misc import chrisfuncs
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.frame import Frame, sum_frames, sum_frames_quadratically
from ...core.basics.configuration import print_mapping
from ...core.tools import stringify
from ...core.tools.formatting import print_files_in_list, print_files_in_path, print_directories_in_path
from ...core.tools import formatting as fmt
from ...core.tools import archive
from ...core.tools import types
from ...modeling.preparation import unitconversion
# -----------------------------------------------------------------

galex_bands = ["NUV", "FUV"]
band_short = {"NUV": "nd", "FUV": "fd"}
galex_fuv_zero_point = 18.82
galex_nuv_zero_point = 20.08

# -----------------------------------------------------------------

class GALEXMosaicMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(GALEXMosaicMaker, self).__init__(*args, **kwargs)

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
        self.convolve_paths = dict()
        self.reproject_paths = dict()
        self.swarp_paths = dict()

        self.noise_paths = dict()
        self.converted_paths = dict()
        self.rebinned_paths = dict()
        self.mosaic_paths = dict()

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
        self.rebin_header_path = None
        self.rebin_wcs = None

        # The mosaic frames and error maps
        self.mosaics = dict()
        self.mosaics_swarp = dict()
        self.error_maps = dict()
        self.relative_error_maps = dict()

        # Create a dictionary for the exposure times for each image
        self.exposure_times = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

        # Reproject
        self.reproject()

        # 9. Mosaic
        self.mosaic()

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

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(GALEXMosaicMaker, self).setup(**kwargs)

        # Show the configuration
        if log.is_debug:
            log.debug("Configuration:")
            print_mapping(self.config, empty_lines=False, indent="  ")
            print("")

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
        #self.root_path = fs.join(fs.home, time.unique_name("GALEX_" + self.config.galaxy_name))
        #fs.create_directory(self.root_path)

        self.root_path = self.config.path

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
            convolve_path = fs.create_directory_in(path, "convolve")
            reproject_path = fs.create_directory_in(path, "reproject")
            swarp_path = fs.create_directory_in(path, "swarp")

            noise_path = fs.create_directory_in(path, "noise")
            converted_path = fs.create_directory_in(path, "converted")
            rebinned_path = fs.create_directory_in(path, "rebinned")
            mosaic_path = fs.create_directory_in(path, "mosaic")

            response_path = fs.create_directory_in(path, "response")
            background_path = fs.create_directory_in(path, "background")
            counts_path = fs.create_directory_in(path, "counts")

            # Set paths
            self.download_paths[band] = download_path
            self.download_observations_paths[band] = download_observations_path
            self.download_response_paths[band] = download_response_path
            self.download_background_paths[band] = download_background_path
            self.download_counts_paths[band] = download_counts_path
            self.convolve_paths[band] = convolve_path
            self.reproject_paths[band] = reproject_path
            self.swarp_paths[band] = swarp_path

            self.noise_paths[band] = noise_path
            self.converted_paths[band] = converted_path
            self.rebinned_paths[band] = rebinned_path
            self.mosaic_paths[band] = mosaic_path

            self.response_paths[band] = response_path
            self.background_paths[band] = background_path
            self.counts_paths[band] = counts_path

            # Debugging
            print_directories_in_path(path)

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

        # Debugging
        log.debug("Cutout center: " + stringify.stringify(self.cutout_center)[1])
        log.debug("Cutout width: " + stringify.stringify(self.cutout_width)[1])

    # -----------------------------------------------------------------

    def get_target_header(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the header and WCS for the mosaic ...")

        # Get the target header
        self.rebin_header, self.rebin_header_path = self.dpdp.get_header_for_galaxy(self.ngc_name, "GALEX", returns=["header", "path"])

        # Create WCS
        self.rebin_wcs = CoordinateSystem(self.rebin_header)

        # Debugging
        log.debug("Rebin WCS:")
        if log.is_debug:
            print("")
            print(self.rebin_wcs)
            print("")

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

            # Inform the user
            log.info("Getting the URLs for the '" + band + "' band ...")

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

        # Inform the user
        log.info("Getting the observation urls in both GALEX bands ...")

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

            # Stop adding more observations if the maximum has been achieved
            #if self.config.max_nobservations_fuv is not None: log.debug(str(len(urls_fuv)) + " observations of " + str(self.config.max_nobservations_fuv) + " (max) reached")
            #if self.config.max_nobservations_nuv is not None: log.debug(str(len(urls_nuv)) + " observations of " + str(self.config.max_nobservations_nuv) + " (max) reached")
            if self.config.max_nobservations_fuv is not None and len(urls_fuv) == self.config.max_nobservations_fuv and nuv_or_fuv == "fuv": continue
            if self.config.max_nobservations_nuv is not None and len(urls_nuv) == self.config.max_nobservations_nuv and nuv_or_fuv == "nuv": continue

            # Loop over all tilenames for this galaxy, find whether the tilename matches the current URL
            for tilename in tilenames:
                if tilename in url:
                    if nuv_or_fuv == "nuv": urls_nuv.append(url)
                    else: urls_fuv.append(url)
                    break  # URL added, so no need to look at the other tilenames for this URL
            #else: raise RuntimeError("Tilename was not found in the DustPedia GALEX tile urls table for the url '" + url + "'")
            #else: log.warning("Tilename was not found in the DustPedia GALEX tile urls table for the url '" + url + "'")

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

    def has_existing_downloads(self, band):

        """
        This fucntion ...
        :param band: 
        :return: 
        """

        if self.config.download_directories is not None and band in self.config.download_directories: return True
        else: return False

    # -----------------------------------------------------------------

    def existing_downloads_path(self, band):

        """
        This function ...
        :param band: 
        :return: 
        """

        if self.config.download_directories is not None and band in self.config.download_directories: return fs.absolute_path(self.config.download_directories[band])
        else: return None

    # -----------------------------------------------------------------

    @property
    def has_manual_selection(self):

        """
        This function ...
        :return: 
        """

        return self.config.manual_selection is not None

    # -----------------------------------------------------------------

    def has_manual_selection_for_band(self, band):

        """
        This function ...
        :param band: 
        :return: 
        """

        return self.config.manual_selection is not None and band in self.config.manual_selection

    # -----------------------------------------------------------------

    def selected_names_for_band(self, band):

        """
        This function ...
        :param band: 
        :return: 
        """

        return self.config.manual_selection[band]

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

            short = band_short[band]

            # Debugging
            log.debug("Downloading the data for the " + band + " band ...")

            # Copy from existing downloads
            if self.has_existing_downloads(band): self.copy_downloads_for_band(band)

            # Otherwise
            else:

                # Parallel execution
                with ParallelTarget(network.download_and_decompress_files, self.config.nprocesses) as target:

                    # Create selection
                    if self.has_manual_selection_for_band(band):

                        # Get selected names
                        selected = self.selected_names_for_band(band)

                        # Loop over the observation URLS
                        observation_urls = []
                        for url in self.observation_urls[band]:

                            # Determine the name
                            name = fs.strip_extension(fs.name(url), double=True).split("-" + short)[0] # .split("-int")[0].split("-")[1] == "nd" else "fuv"

                            if name not in selected: continue
                            observation_urls.append(url)

                        # Loop over the response URLS
                        response_urls = []
                        for url in self.response_urls[band]:

                            # Determine the name
                            name = fs.strip_extension(fs.name(url), double=True).split("-" + short)[0] #

                            if name not in selected: continue
                            response_urls.append(url)

                        # Loop over the background URLs
                        background_urls = []
                        for url in self.background_urls[band]:

                            # Determine the name
                            name = fs.strip_extension(fs.name(url), double=True).split("-" + short)[0]

                            if name not in selected: continue
                            background_urls.append(url)

                        # Loop over the counts URLs
                        counts_urls = []
                        for url in self.counts_urls[band]:

                            # Determine the name
                            name = fs.strip_extension(fs.name(url), double=True).split("-" + short)[0]

                            if name not in selected: continue
                            counts_urls.append(url)

                    else:

                        observation_urls = self.observation_urls[band]
                        response_urls = self.response_urls[band]
                        background_urls = self.background_urls[band]
                        counts_urls = self.counts_urls[band]

                    # Call the target function
                    target(observation_urls, self.download_observations_paths[band], info="observations")
                    target(response_urls, self.download_response_paths[band], info="response maps")
                    target(background_urls, self.download_background_paths[band], info="background maps")
                    target(counts_urls, self.download_counts_paths[band], info="count maps")

            # Debugging
            print_files_in_path(self.download_observations_paths[band])
            print_files_in_path(self.download_response_paths[band])
            print_files_in_path(self.download_background_paths[band])
            print_files_in_path(self.download_counts_paths[band])

    # -----------------------------------------------------------------

    def copy_downloads_for_band(self, band):

        """
        This function ...
        :return: 
        """

        # Debugging
        log.debug("Existing data for band in: " + self.existing_downloads_path(band))

        extensions = ["fits"]
        extensions.extend(archive.extensions)

        not_contains = ["meta", "overlap"]

        # OBSERVATIONS

        existing_observations_path = fs.join(self.existing_downloads_path(band), "observations")

        if self.has_manual_selection_for_band(band): selected_names = self.selected_names_for_band(band)
        else: selected_names = None

        print("OBSERVATIONS:")
        fmt.print_files_in_path(existing_observations_path)

        paths = fs.files_in_path(existing_observations_path, not_contains=not_contains, extension=extensions, contains=selected_names, contains_operator="OR")
        fmt.print_files_in_list(paths, "existing observations", only_name=True)

        # Copy and decompress
        fs.copy_and_decompress_files(paths, self.download_observations_paths[band])



        # RESPONSE

        existing_response_path = fs.join(self.existing_downloads_path(band), "response")

        print("RESPONSE:")
        fmt.print_files_in_path(existing_response_path)

        paths = fs.files_in_path(existing_response_path, not_contains=not_contains, extension=extensions, contains=selected_names, contains_operator="OR")
        fmt.print_files_in_list(paths, "existing response maps", only_name=True)

        # Copy and decompress
        fs.copy_and_decompress_files(paths, self.download_response_paths[band])



        # BACKGROUND

        existing_background_path = fs.join(self.existing_downloads_path(band), "background")

        print("BACKGROUND:")
        fmt.print_files_in_path(existing_background_path)

        paths = fs.files_in_path(existing_background_path, not_contains=not_contains, extension=extensions, contains=selected_names, contains_operator="OR")
        fmt.print_files_in_list(paths, "existing background maps", only_name=True)

        # Copy and decompress
        fs.copy_and_decompress_files(paths, self.download_background_paths[band])



        # COUNTS

        existing_counts_path = fs.join(self.existing_downloads_path(band), "counts")

        print("COUNTS:")
        fmt.print_files_in_path(existing_counts_path)

        paths = fs.files_in_path(existing_counts_path, not_contains=not_contains, extension=extensions, contains=selected_names, contains_operator="OR")
        fmt.print_files_in_list(paths, "existing count maps")

        # Copy and decompress
        fs.copy_and_decompress_files(paths, self.download_counts_paths[band])

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
            for band in self.config.bands:

                # ALREADY MANUAL SELECTION
                if self.has_manual_selection_for_band(band): continue

                target(self.ngc_name, band, self.download_observations_paths[band], self.cutout_center, self.cutout_width, mode='point')

        # Debugging
        for band in self.config.bands: print_files_in_path(self.download_observations_paths[band])

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
            for band in self.config.bands:

                # ALREADY MANUAL SELECTION
                if self.has_manual_selection_for_band(band): continue

                target(self.ngc_name, self.download_observations_paths[band], self.cutout_center, self.cutout_width, band)

        # Debugging
        for band in self.config.bands: print_files_in_path(self.download_observations_paths[band])

    # -----------------------------------------------------------------

    def clean(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cleaning GALEX tiles ...")

        # Loop over the bands
        for band in self.config.bands:

            # Debugging
            log.debug("Cleaning GALEX " + band + " tiles ...")

            # Determine paths for this band
            response_path = self.download_response_paths[band]
            convolve_path = self.convolve_paths[band]
            background_path = self.download_background_paths[band]
            reproject_path = self.reproject_paths[band]

            # Debugging
            log.debug("Response path: " + response_path)
            log.debug("Convolve path: " + convolve_path)
            log.debug("Background path: " + background_path)
            log.debug("Reproject path: " + reproject_path)

            # Parallel execution
            with ParallelTarget(clean_galex_tile, self.config.nprocesses) as target:

                # Loop over the raw files
                # raw_file_path, response_path, convolve_path, background_path, reproject_path
                for path in fs.files_in_path(self.download_observations_paths[band]): target(path, response_path, convolve_path, background_path, reproject_path)

            # Debugging
            print_files_in_path(convolve_path)
            print_files_in_path(reproject_path)

    # -----------------------------------------------------------------

    def level(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Peforming background-leveling on the observations ...")

        #print("Performing background-leveling on the observations ...")

        target_suffix = "int.fits"

        # Loop over the bands
        for band in self.config.bands:

            reproject_path = self.reproject_paths[band]
            convolve_path = self.convolve_paths[band]

            # Make list of files in target directory that have target suffix
            allfile_list = os.listdir(reproject_path)
            fitsfile_list = []
            for allfile in allfile_list:
                if target_suffix in allfile:
                    fitsfile_list.append(allfile)

            #
            offsets = [None] * len(fitsfile_list)
            results = [None] * len(fitsfile_list)

            # Parallel execution
            with ParallelTarget(determine_background_level, self.config.nprocesses) as target:

                # Loop over each file
                for i in range(len(fitsfile_list)):

                    filename = fitsfile_list[i]

                    # Call the target function
                    result = target(filename, convolve_path)

                    # Set the result handle
                    results[i] = result

            level_ref = None

            # For
            for index in range(len(results)):

                results[index].request()
                output = results[index].output

                if types.is_real_type(output) or types.is_integer_type(output): level = float(output)
                elif types.is_sequence(output): level = output[0]
                else: raise RuntimeError("Output type not recognized: " + str(output))
                #level = output[0]

                if index == 0:
                    level_ref = level
                    continue

                average_offset = level_ref - level
                offsets[index] = average_offset

            # Parallel execution
            with ParallelTarget(level_galex_map, self.config.nprocesses) as target:

                # Loop over the files
                for i in range(1, len(fitsfile_list)):

                    filename = fitsfile_list[i]
                    offset = offsets[i]

                    # Call the target function
                    target(filename, offset, reproject_path)

            # Debugging
            print_files_in_path(reproject_path)

    # -----------------------------------------------------------------

    def get_exposure_times(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting exposure times ...")

        # Loop over the bands
        for band in self.config.bands:

            # Debugging
            log.debug("Getting exposure times for the " + band + " band ...")

            # Initialize dictionary
            self.exposure_times[band] = dict()

            reproject_path = self.reproject_paths[band]

            # Determine how the files are named
            #filename_ends = "-" + band_short[band] + "-int"  # -fd-int for FUV, -nd-int for NUV
            filename_ends = "-" + band_short[band] + "-int.wgt"

            # Debugging
            log.debug("Looking for files with a name ending on '" + filename_ends + "' ...")

            # Get the exposure time for each image
            for path, name in fs.files_in_path(reproject_path, extension="fits", contains=filename_ends, returns=["path", "name"]):

                # Debugging
                log.debug("Getting the exposure time of the " + name + " image ...")

                # Open the header
                header = getheader(path)

                # Search for the exposure time
                exp_time = get_total_exposure_time(header)

                # Debugging
                log.debug("Found an exposure time of " + str(exp_time) + " ...")

                # Set the exposure time
                image_name = name.split(filename_ends)[0]
                self.exposure_times[band][image_name] = exp_time

    # -----------------------------------------------------------------

    def reproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Reprojecting ...")

        # Loop over the bands
        for band in self.config.bands:

            reproject_path = self.reproject_paths[band]
            swarp_path = self.swarp_paths[band]

            # Generate metatable
            metatable_path = mosaicing.generate_meta_file(reproject_path)

            # Reproject image and weight prior to coaddition
            mosaicing.reproject(reproject_path, swarp_path, metatable_path, self.rebin_header_path)

            # Rename reprojected files for SWarp
            for listfile in os.listdir(swarp_path):
                if '_area.fits' in listfile: os.remove(fs.join(swarp_path, listfile))
                elif 'hdu0_' in listfile: os.rename(fs.join(swarp_path, listfile), fs.join(swarp_path, listfile.replace('hdu0_', '')))

            # Debugging
            print_files_in_path(swarp_path)

    # -----------------------------------------------------------------

    def mosaic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making GALEX mosaic for " + self.ngc_name + " and map of relative poisson errors ...")

        # Mosaic using Swarp
        try: self.mosaic_swarp()
        except Exception as e:
            traceback.print_exc()
            log.warning("Something went wrong during the SWARPing, but that's no problem ...")

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

        # Get the pixelscale
        pixelscale_arcsec = self.dpdp.get_pixelscale_for_instrument("GALEX").to("arcsec").value

        # Parallel execution
        results = dict()
        with ParallelTarget(mosaic_with_swarp, self.config.nprocesses) as target:

            # Loop over the bands and execute
            # band, width, pix_size, swarp_path, center
            for band in self.config.bands: results[band] = target(band, self.cutout_width, pixelscale_arcsec, self.swarp_paths[band], self.cutout_center)

        # Loop over bands, LOAD SWARP RESULT, PRESUMABLY IN COUNTS/S
        for band in self.config.bands:

            results[band].request()
            output = results[band].output
            swarp_result_path = output

            # Debugging
            print_files_in_path(self.swarp_paths[band])

            # Load the resulting frame
            out_image = Frame.from_file(swarp_result_path)
            out_image.unit = "count/s"

            # Overwrite invalid data
            out_image[out_image == 0] = np.NaN
            out_image[out_image < -1E3] = np.NaN
            out_image[out_image <= 1E-8] = 0

            # Set the mosaic
            self.mosaics_swarp[band] = out_image

    # -----------------------------------------------------------------

    def mosaic_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making mosaic with PTS ...")

        image_names_bands = dict()

        # Loop over the bands
        for band in self.config.bands:

            # Get list of images
            image_names = self.get_images_for_mosaic(band)

            # Debugging
            print_files_in_list(image_names, "image_names (for mosaic pts)")

            # Make noise map
            self.make_noise(band, image_names)

            # Convert
            self.convert_to_count_s_sr(band, image_names)

            # Rebin
            self.rebin_for_mosaic(band, image_names)

            # Set the image names
            image_names_bands[band] = image_names

        # Combine
        self.combine(image_names_bands)

    # -----------------------------------------------------------------

    def get_images_for_mosaic(self, band):

        """
        This function ...
        :return:
        """

        # The list of image names to be used for the mosaic
        image_names_for_mosaic = []

        swarp_path = self.swarp_paths[band]

        filename_ends = "-" + band_short[band] + "-int"

        # Loop over the files in the temp_swarp_path directory, where the tiles from temp_reproject_path are _peakand saved to
        for filename in os.listdir(swarp_path):

            # SKIP WEIGHT FILES !
            if not filename.endswith("-int.fits"): continue

            # Get the image name
            image_name = filename.split(filename_ends)[0]

            # Add the image name to the list
            image_names_for_mosaic.append(image_name)

        # Return the list of image names
        return image_names_for_mosaic

    # -----------------------------------------------------------------

    def make_noise(self, band, image_names):

        """
        This function ...
        :param band:
        :param image_names:
        :return:
        """

        # Inform the user
        log.info("Making noise maps for the " + band + " ...")

        # The path to the directory with all the tiles in counts for the current band (FUV or NUV)
        counts_path = self.download_counts_paths[band]
        noise_path = self.noise_paths[band]

        # Execute in parallel
        with ParallelTarget(make_noise_map_in_cps, self.config.nprocesses) as target:

            # Inform the user
            log.info("Creating maps of the poisson noise for each GALEX tile in counts per second ...")

            # Loop over the file names
            for image_name in image_names:

                # Get the exposure time
                exposure_time = self.exposure_times[band][image_name]

                # Make noise maps in count/s
                # band, image_name, counts_path_band, temp_noise_path, exposure_time
                target(band, image_name, counts_path, noise_path, exposure_time)

        # Debugging
        print_files_in_path(noise_path)

    # -----------------------------------------------------------------

    def convert_to_count_s_sr(self, band, image_names):

        """
        This function ...
        :param band:
        :param image_names:
        :return:
        """

        # Inform the user
        log.info("Converting to counts/s/sr for the " + band + "...")

        # Paths
        reproject_path = self.reproject_paths[band]
        noise_path = self.noise_paths[band]
        converted_path = self.converted_paths[band]

        # Parallel execution
        with ParallelTarget(convert_frame_and_error_map_to_per_solid_angle, self.config.nprocesses) as target:

            # DO THE UNIT CONVERSION
            for image_name in image_names:

                # Convert to counts/s/sr
                # band, image_name, temp_reproject_path, temp_noise_path, temp_converted_path
                target(band, image_name, reproject_path, noise_path, converted_path)

        # Debugging
        print_files_in_path(converted_path)

    # -----------------------------------------------------------------

    def rebin_for_mosaic(self, band, image_names):

        """
        This function ...
        :param band:
        :param image_names:
        :return:
        """

        # Inform the user
        log.info("Rebinning frames for mosaic for the " + band + " ...")

        # Get paths
        reproject_path = self.reproject_paths[band]
        converted_path = self.converted_paths[band]
        rebinned_path = self.rebinned_paths[band]

        # Parallel execution
        with ParallelTarget(rebin_frame_and_error_map, self.config.nprocesses) as target:

            # Loop over all FITS files
            for path, name in fs.files_in_path(converted_path, extension="fits", returns=["path", "name"]):

                # Rebin frames and error maps
                target(path, rebinned_path, self.rebin_header_path)

        # Inform the user
        log.info("Rebinning weight maps ...")

        # Parallel execution
        with ParallelTarget(rebin_weight_map, self.config.nprocesses) as target:

            # Loop over the image_names
            for image_name in image_names:

                # Rebin the weight maps
                target(band, image_name, reproject_path, rebinned_path, self.rebin_header_path)

        # Debugging
        print_files_in_path(rebinned_path)

    # -----------------------------------------------------------------

    def combine(self, image_names_bands):

        """
        This function ...
        :param image_names_bands:
        :return:
        """

        # Inform the user
        log.info("Combining the maps ...")

        results = dict()

        # Parallel execution
        with ParallelTarget(combine_frames_and_error_maps, self.config.nprocesses) as target:

            # Loop over the bands
            for band in self.config.bands:

                rebinned_path = self.rebinned_paths[band]
                mosaic_path = self.mosaic_paths[band]

                # image_names_for_mosaic, temp_rebinned_path, temp_mosaic_path, wcs
                # RETURNS: mosaic_frame_path, mosaic_errors_path
                result = target(image_names_bands[band], rebinned_path, mosaic_path, self.rebin_header_path)
                results[band] = result

        # Inform the user
        log.info("Finishing ...")

        # Load the mosaics
        for band in self.config.bands:

            results[band].request()
            output = results[band].output

            mosaic_path, error_path = output[0], output[1]
            mosaic = Frame.from_file(mosaic_path)
            errors = Frame.from_file(error_path)

            # SOME THINGS
            mosaic[mosaic.data == 0] = np.NaN
            mosaic[mosaic.data < -1E3] = np.NaN
            mosaic[mosaic.data <= 1E-8] = 0

            # CALCULATE RELATIVE POISSON ERROR MAP
            relerrors = errors / mosaic
            relerrors[relerrors < 0.] = 0.0  # set negative values for relative error map to zero
            relerrors.replace_nans(0.0)  # set NaN values (because mosaic was zero) to zero

            # CONVERT TO JANSKY / PIX

            # FROM COUNT / S TO AB MAG:
            if band == "FUV":
                mosaic = Frame(galex_fuv_zero_point - (2.5 * np.log10(mosaic)), wcs=mosaic.wcs)
                errors = Frame(galex_fuv_zero_point - (2.5 * np.log10(errors)), wcs=errors.wcs)
            elif band == "NUV":
                mosaic = Frame(galex_nuv_zero_point - (2.5 * np.log10(mosaic)), wcs=mosaic.wcs)
                errors = Frame(galex_nuv_zero_point - (2.5 * np.log10(errors)), wcs=errors.wcs)

            # FROM AB MAG TO FLUX (JANSKY):
            mosaic = Frame(unitconversion.ab_to_jansky(mosaic), wcs=mosaic.wcs)
            errors = Frame(unitconversion.ab_to_jansky(errors), wcs=errors.wcs)

            # SET IMAGE UNITS
            mosaic.unit = "Jy/pix"
            errors.unit = "Jy/pix"

            # SET IMAGE FILTER
            mosaic.filter = band
            errors.filter = band

            # Set
            self.mosaics[band] = mosaic
            self.error_maps[band] = errors
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
            id_string = self.ngc_name + "_GALEX_" + band
            #path = fs.join(output_path, id_string + ".fits")
            path = self.output_path_file(id_string + ".fits")

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
            id_string = self.ngc_name + "_GALEX_" + band
            #path = fs.join(output_path, id_string + "_errors.fits")
            path = self.output_path_file(id_string + "_errors.fits")

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
            id_string = self.ngc_name + "_GALEX_" + band
            #path = fs.join(output_path, id_string + "_relerrors.fits")
            path = self.output_path_file(id_string + "_relerrors.fits")

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
    circle_overlap_path = fs.join(tiles_path, "overlap_circle.dat")

    # Debuggging
    log.debug("Circle overlap path: " + circle_overlap_path)

    # Get overlapping file paths
    #ra = center.ra.to("deg").value
    #dec = center.dec.to("deg").value
    #width = width.to("deg").value
    ra = center.ra
    dec = center.dec

    #overlapping_file_paths = mosaicing.generate_overlapping_file_paths(tiles_path, ra, dec, meta_path, mode="circle", radius=(0.5 * width) * (2.0 ** 0.5))

    # Generate circle overlap file
    mosaicing.generate_overlap_file(tiles_path, ra, dec, meta_path, mode="circle", radius=(0.5 * width) * (2.0 ** 0.5), overlap_path=circle_overlap_path)

    # Get overlapping file paths
    overlapping_file_paths = mosaicing.get_overlapping_file_paths(circle_overlap_path)

    fmt.print_files_in_list(overlapping_file_paths, "overlapping files", only_name=True)

    ra_deg = ra.to("deg").value
    dec_deg = dec.to("deg").value

    # Check
    #if len(overlapping_file_paths.shape) == 0: overlapping_file_paths = [overlapping_file_paths.tolist()]
    #for overlapping_file_path in overlapping_file_paths: fs.copy_file(overlapping_file_path, temp_raw_path)

    # NEW: FILTER IN-PLACE
    removed = []
    for path in fs.files_in_path(tiles_path, extension="fits"):

        # Debugging
        log.debug("Checking image " + path + " ...")

        if path in overlapping_file_paths: continue
        else:
            log.debug("Removing image " + path + " ...")
            fs.remove_file(path)
            removed.append(path)

    print("REMOVED:")
    fmt.print_files_in_list(removed, "removed files", only_name=True)

    temp_raw_path = tiles_path

    # Uncompress .fits.gz files
    # [os.system('gunzip '+ listfile) for listfile in os.listdir(raw_in_temp_dir)]

    # Ensure that at least one of the raw GALEX tiles has actual flux coverage at location of source
    # raw_files = os.listdir(temp_raw_path) ## THERE WAS A 'deg' FILE IN THE DIRECTORY AS WELL WHICH WAS NOT A FITS FILE BUT A PLAIN TEXT HEADER FILE, AND SO THIS WASN'T WORKING
    raw_files = fs.files_in_path(temp_raw_path, extension="fits", returns="name", extensions=True)
    #raw_files = [name + ".fits" for name in raw_files]

    fmt.print_files_in_list(raw_files, "raw files", only_name=True)

    coverage = False
    covering = []

    # Loop over the files
    for raw_file in raw_files:

        # Debugging
        log.debug("Checking coverage at galaxy center for image " + raw_file + " ...")

        # Determine path
        path = fs.join(temp_raw_path, raw_file)

        # Read in map
        in_fitsdata = open_fits(path)
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_fitsdata.close()

        # Locate pixel coords
        in_wcs = WCS(in_header)

        print("WCS:")
        print(in_wcs)

        location_pix = in_wcs.wcs_world2pix(np.array([[np.float(ra_deg), np.float(dec_deg)]]), 0)[0]
        pix_i, pix_j = location_pix[1], location_pix[0]

        log.debug("Pixel coordinate: " + str(pix_i) + ", " + str(pix_j))

        # Evalulate coverage at location, and proceed accordingly
        check_coordinates = [pix_i - 10, pix_i + 11, pix_j - 10, pix_j + 11]
        print("check coordinates: " + str(check_coordinates))
        if True in [coord <= 0 for coord in check_coordinates]:
            print("FIRST CONDITION WAS EVALUATED AS TRUE, SO IMAGE IS NOT COUNTED AS COVERING")
            continue
        try:
            print(in_image.shape)
            print(pix_i - 10, pix_i + 11, pix_j - 10, pix_j + 11)
            image_slice = in_image[pix_i - 10:pix_i + 11, pix_j - 10:pix_j + 11]
        except:
            print("AN ERROR OCCURED TRYING TO GET THE IMAGE SLICE, SO THE IMAGE IS NOT COUNTED AS COVERING")
            continue

        pixels_y, pixels_x = np.where(image_slice > 0)

        #if pixels_y.shape[0] > 0:
        #if pixels_y.size > 0: # should be the same
        if len(pixels_y) > 0: # should be the same

            coverage = True
            covering.append(path)
            #break

    # Debugging
    print("COVERING:")
    fmt.print_files_in_list(covering, "covering", only_name=True)

    # No coverage
    if not coverage:
        log.warning("It seems that there is no flux coverage for the GALEX " + band + " band for " + galaxy_name)
        #raise RuntimeError('No GALEX ' + band + ' coverage for ' + galaxy_name)

    # Return ...
    #return raw_files

# -----------------------------------------------------------------

def clean_galex_tile(raw_file_path, response_path, convolve_path, background_path, reproject_path):

    """
    Function to clean GALEX tiles and create exposure maps
    :param raw_file_path: raw file path
    :param response_path:
    :param convolve_path:
    :param background_path:
    :param reproject_path:
    :return:
    """

    raw_file_name = fs.name(raw_file_path)

    # Inform the user ...
    log.info("Cleaning map " + raw_file_path + " ...")

    #print("Cleaning map " + raw_file_path + " ...")

    # Create subdirectories
    #temp_raw_path = fs.create_directory_in(temp_path_band, "raw")
    #temp_convolve_path = fs.create_directory_in(temp_path_band, "convolve")

    # Read in image
    #raw_file_path = fs.join(temp_raw_path, raw_file)
    if raw_file_path.endswith('.dat'): return
    in_fitsdata = open_fits(raw_file_path, ignore_missing_end=True)
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    out_image = in_image.copy()

    # Load and align response map
    rr_path = fs.join(response_path, raw_file_name.replace('-int.fits','-rr.fits'))
    print("rr_path: " + rr_path)
    if not fs.is_file(rr_path): print("NOT PRESENT")
    rr_fitsdata = open_fits(rr_path)
    rr_image = rr_fitsdata[0].data
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean image using response map
    out_image[ np.where( rr_image <= 1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_path = fs.join(background_path, raw_file_name.replace('-int.fits','-skybg.fits'))
    print("bg_path: " + bg_path)
    if not fs.is_file(bg_path): print("NOT PRESENT")
    bg_fitsdata = open_fits(bg_path)
    bg_image = bg_fitsdata[0].data
    bg_zoom = np.float(out_image.shape[0]) / np.float(bg_image.shape[0])
    bg_image = interpolation.zoom(bg_image, bg_zoom, order=0)

    # Clean image using sky background map
    out_image[np.where(bg_image <= 1E-10)] = np.NaN

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
    #out_hdulist.writeto(fs.join(temp_raw_path, raw_file_name), clobber=True)

    out_file_path = fs.join(reproject_path, raw_file_name)
    #out_hdulist.writeto(raw_file_path, clobber=True)
    out_hdulist.writeto(out_file_path, clobber=True)

    # Create convolved version of map, for later use in background-matching
    """
    if np.isnan(out_image).sum()==0:
        conv_image = scipy.ndimage.filters.gaussian_filter(out_image, 20)
    else:
    """

    # Do the convolution
    kernel = Tophat2DKernel(10)
    conv_image = convolve_fft(out_image, kernel, nan_treatment='fill', normalize_kernel=True, allow_huge=True) #, ignore_edge_zeros=False, interpolate_nan=True, normalize_kernel=True)

    # Write the convolved map
    temp_convolve_image_path = fs.join(convolve_path, raw_file_name)                  ## NEW
    print("Writing convolved map to " + temp_convolve_image_path)
    if fs.is_file(temp_convolve_image_path): fs.remove_file(temp_convolve_image_path) ## NEW
    writeto(temp_convolve_image_path, conv_image, in_header)

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[ np.where( np.isnan(out_image)==False ) ] = (float(in_header['exptime']))**0.5 # SQUARE ROOT OF THE EXPOSURE TIME
    exp_hdu = PrimaryHDU(data=exp_image, header=in_header)
    exp_hdulist = HDUList([exp_hdu])

    # Write the weight map
    temp_reproject_weight_path = fs.join(reproject_path, raw_file_name.replace('.fits','.wgt.fits'))  ## NEW
    print("Writing weight map to " + temp_reproject_weight_path)
    if fs.is_file(temp_reproject_weight_path): fs.remove_file(temp_reproject_weight_path)              ## NEW
    exp_hdulist.writeto(temp_reproject_weight_path)

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

def determine_background_level(filename, convolve_path):

    """
    original name: GALEX_Zero
    Set a set of maps to the same level
    :return:
    """

    convfile_dir = convolve_path

    # Read in corresponding map from directory containing convolved images
    conv_filepath = fs.join(convfile_dir, filename)
    fitsdata_conv = open_fits(conv_filepath)
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

    # Return the level
    return level

# -----------------------------------------------------------------

def level_galex_map(filename, average_offset, reproject_path):

    """
    This fcuntion ...
    :param filename:
    :param average_offset:
    :param reproject_path
    :return:
    """

    # Inform the user
    log.info("Matching background of map " + filename + " ...")

    fitsfile_dir = reproject_path

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

    filepath = fs.join(fitsfile_dir, filename)

    # Read in unconvolved file, and apply offset
    fitsdata_in = open_fits(filepath)
    image_in = fitsdata_in[0].data
    header_in = fitsdata_in[0].header
    fitsdata_in.close()
    image_out = image_in + average_offset
    #print 'Map mean of '+fitsfile_list[i]+' changed from '+str(np.nanmean(image_in))+' to '+str(np.nanmean(image_out))

    # Save corrected file
    image_out_hdu = PrimaryHDU(data=image_out, header=header_in)
    image_out_hdulist = HDUList([image_out_hdu])
    image_out_hdulist.writeto(filepath, clobber=True)

# -----------------------------------------------------------------

def make_noise_map_in_cps(band, image_name, counts_path_band, temp_noise_path, exposure_time):

    """
    This function ...
    :param band:
    :param image_name:
    :param counts_path_band:
    :param temp_noise_path:
    :param exposure_time:
    :return:
    """

    # Get the path to the file in counts
    counts_filepath = fs.join(counts_path_band, image_name + "-" + band_short[band] + "-cnt.fits")

    # Load the counts map
    counts_frame = Frame.from_file(counts_filepath)
    counts_frame.unit = "ct"

    # Calculate the poisson frame
    poisson = Frame(np.sqrt(counts_frame.data))   # calculate the poisson error (in counts) in every pixel
    poisson.wcs = counts_frame.wcs
    poisson.unit = "ct"

    # Get the exposure time in seconds
    #exposure_time = exposure_times[image_name]

    # Convert the poisson frame to counts/second
    poisson /= exposure_time
    poisson.unit = "ct/s"

    # Determine the path to the poisson noise map in counts per second
    new_path = fs.join(temp_noise_path, image_name + ".fits")

    # Save the poisson noise map in counts per second
    poisson.saveto(new_path)

# -----------------------------------------------------------------

def convert_frame_and_error_map_to_per_solid_angle(band, image_name, temp_reproject_path, temp_noise_path, temp_converted_path):

    """
    This function ...
    :param band:
    :param image_name:
    :param temp_reproject_path:
    :param temp_noise_path:
    :param temp_converted_path:
    :return:
    """

    # Inform the user
    log.info("Converting the units to luminosity per solid angle ...")

    # Debugging
    log.debug("Converting " + image_name + " frame and error map to count / s / sr...")

    # Determine ...
    filename_ends = "-" + band_short[band] + "-int"  # -fd-int for FUV, -nd-int for NUV

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

def rebin_frame_and_error_map(filepath, temp_rebinned_path, header_path):

    """
    This function ...
    :param filepath:
    :param temp_rebinned_path:
    :param header_path:
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
    #for path, name in fs.files_in_path(temp_converted_path, extension="fits", returns=["path", "name"]):

    filename = fs.strip_extension(fs.name(filepath))

    # Load the image
    frame = Frame.from_file(filepath)

    # Rebin
    frame.rebin(rebin_wcs)

    # Determine the new path
    new_path = fs.join(temp_rebinned_path, filename + ".fits")

    # Save
    frame.saveto(new_path)

# -----------------------------------------------------------------

def rebin_weight_map(band, image_name, temp_reproject_path, temp_rebinned_path, header_path):

    """
    This function ...
    :param band:
    :param image_name:
    :param temp_reproject_path:
    :param temp_rebinned_path:
    :param header_path:
    :return:
    """

    # Determine how the files are named
    filename_ends = "-" + band_short[band] + "-int"  # -fd-int for FUV, -nd-int for NUV

    # GET REBIN WCS
    rebin_header = Header.fromtextfile(header_path)
    rebin_wcs = CoordinateSystem(rebin_header)

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

def combine_frames_and_error_maps(image_names_for_mosaic, temp_rebinned_path, temp_mosaic_path, header_path):

    """
    This function ...
    :param image_names_for_mosaic:
    :param temp_rebinned_path:
    :param temp_mosaic_path:
    :param header_path:
    :return:
    """

    # GET REBIN WCS
    rebin_header = Header.fromtextfile(header_path)
    wcs = CoordinateSystem(rebin_header)

    #primary_frames = []
    #error_frames = []
    weight_frames = []
    weighted_frames = []
    frames = []

    # Loop over the images to be used for the mosaic
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
        #errors_weighted = errors / weights

        # Convert the units from counts/s/sr to counts/s
        pixelsr = frame_weighted.pixelarea.to("sr").value
        frame_weighted *= pixelsr
        frame *= pixelsr
        #errors_weighted *= pixelsr

        # Create mask where the weights are nans
        mask = weights.nans

        # Set zero
        frame_weighted[mask] = 0.0
        frame[mask] = 0.0
        #errors_weighted[mask] = 0.0
        weights[mask] = 0.0

        # Add to the list
        #primary_frames.append(frame_weighted)
        #error_frames.append(errors_weighted)
        #weight_frames.append(weights)

        ## NEW: FIX MISTAKE
        weight_frames.append(weights)
        weighted_frames.append(frame_weighted)
        frames.append(frame)

    # Calculate denominator for weighted average (mosaicing)
    normalization = sum_frames(*weight_frames)

    # CALCULATE THE MOSAIC FRAME IN COUNTS/S
    #mosaic_frame = sum_frames(*primary_frames) / normalization
    #mosaic_frame.wcs = wcs
    ## NEW
    mosaic_frame = sum_frames(*weighted_frames) / normalization
    mosaic_frame.wcs = wcs
    mosaic_frame.unit = "ct/s"

    # CALCULATE THE MOSAIC ERROR MAP IN COUNTS/S
    #mosaic_errormap = sum_frames_quadratically(*error_frames) / normalization
    #mosaic_errormap.wcs = wcs
    ## NEW
    frames_sum = sum_frames(*frames)
    mosaic_errormap = Frame(np.sqrt(frames_sum.data)) / normalization
    mosaic_errormap.wcs = wcs
    mosaic_errormap.unit = "ct/s"

    ## DONE

    # SAVE THE MOSAIC IN COUNTS/S
    mosaic_path = fs.join(temp_mosaic_path, "mosaic.fits")
    mosaic_frame.saveto(mosaic_path)

    # SAVE THE MOSAIC ERROR MAP IN COUNTS/S
    mosaic_error_path = fs.join(temp_mosaic_path, "mosaic_errors.fits")
    mosaic_errormap.saveto(mosaic_error_path)

    # Return the mosaic path and the mosaic error map path in counts/s
    return mosaic_path, mosaic_error_path

# -----------------------------------------------------------------

def mosaic_with_swarp(band, width, pix_size, swarp_path, center):

    """
    This function ...
    :param band:
    :param width:
    :param pix_size:
    :param swarp_path:
    :param center:
    :return:
    """

    # Use SWarp to co-add images weighted by their error maps
    log.info("Co-adding GALEX " + band + " maps ...")

    # Determine image width in pixels
    image_width_pixels = str(int((float(width.to("deg").value) * 3600.) / pix_size))

    # Do the mosaicing with Swarp
    return mosaicing.mosaic(swarp_path, band, center, image_width_pixels)

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

    import montage_wrapper as montage

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
