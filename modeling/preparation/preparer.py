#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.preparer Contains the DataPreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...magic.core.image import Image
from ...magic.region.list import PixelRegionList
from .component import PreparationComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.dataset import DataSet
from ...core.basics.containers import NamedFileList
from ...magic.services.extinction import GalacticExtinction
from ...core.filter.filter import parse_filter
from ...magic.sources.extractor import SourceExtractor
from ...magic.sky.skysubtractor import SkySubtractor
from ...dustpedia.core.properties import DustPediaProperties
from ...dustpedia.core.photometry import DustPediaPhotometry
from ...magic.basics.mask import Mask as oldMask
from ...magic.core.frame import Frame, sum_frames_quadratically
from . import unitconversion
from ...magic.region.point import PixelPointRegion
from ...magic.basics.coordinate import PixelCoordinate
from ...magic.region.ellipse import PixelEllipseRegion
from ...core.basics.composite import SimplePropertyComposite
from ...core.tools.serialization import write_dict
from ...magic.core.mask import Mask
from ...core.remote.remote import Remote

# -----------------------------------------------------------------

# Define names for the resulting images
initialized_name = "initialized.fits"
extracted_name = "extracted.fits"
corrected_name = "corrected_for_extinction.fits"
subtracted_name = "sky_subtracted.fits"
with_errors_name = "with_errormaps.fits"
result_name = "result.fits"

# Sources and sky directories names
sources_name = "sources"
sky_name = "sky"

# Statistics file
statistics_name = "statistics.txt"

# -----------------------------------------------------------------

# STEPS:

# Source extraction
# Correct for galactic extinction
# Sky SUBTRACTION
# CALCULATE ERROR MAPS
# CONVERT UNITS

steps = ["extraction", "extinction", "subtraction", "errormaps", "units"]

# -----------------------------------------------------------------

status_list = ["initialized", "extracted", "corrected", "subtracted", "with_errors", "result"]

# -----------------------------------------------------------------

filename_for_status = dict()
filename_for_status["initialized"] = initialized_name
filename_for_status["extracted"] = extracted_name
filename_for_status["corrected"] = corrected_name
filename_for_status["subtracted"] = subtracted_name
filename_for_status["with_errors"] = with_errors_name
filename_for_status["result"] = result_name

# -----------------------------------------------------------------

def filenames_before_status(status):

    """
    This function ...
    :param status: 
    :return: 
    """

    filenames = []

    for status_i in status_list:

        if status_i == status: break
        else: filenames.append(filename_for_status[status_i])

    return filenames

# -----------------------------------------------------------------

def status_to_steps(status):

    """
    This function ...
    :param status: 
    :return: 
    """

    if status == "initialized": return []
    elif status == "extracted": return steps_before_and_including("extraction")
    elif status == "corrected": return steps_before_and_including("extinction")
    elif status == "subtracted": return steps_before_and_including("subtraction")
    elif status == "with_errors": return steps_before_and_including("errormaps")
    elif status == "result": return steps_before_and_including("units")
    else: raise ValueError("Invalid status: '" + status + "'")

# -----------------------------------------------------------------

def steps_before(step):

    """
    This function ...
    :param step: 
    :return: 
    """

    # Check
    if step not in steps: raise ValueError("Invalid step: '" + step + "'")

    the_steps = []
    for stepi in steps:
        if stepi == step: break
        else: the_steps.append(stepi)
    return the_steps

# -----------------------------------------------------------------

def steps_before_and_including(step):

    """
    This function ...
    :param step: 
    :return: 
    """

    return steps_before(step) + [step]

# -----------------------------------------------------------------

class PreparationStatistics(SimplePropertyComposite):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PreparationStatistics, self).__init__()

        # Define properties
        self.add_property("nsources", "integer", "total number of sources")
        self.add_property("ngalaxy_sources", "integer", "number of galaxy sources")
        self.add_property("nstar_sources", "integer", "number of star sources")
        self.add_property("nother_sources", "integer", "number of other sources")
        self.add_property("nforeground", "integer", "number of foreground sources")
        self.add_property("nfailed", "integer", "number of failed extractions")
        self.add_property("nsuccess", "integer", "number of succesful extractions")
        self.add_property("nwith_saturation", "integer", "number of stars with a saturation source")
        self.add_property("galaxy_regions", "boolean", "has galaxy regions")
        self.add_property("star_regions", "boolean", "has star regions")
        self.add_property("saturation_regions", "boolean", "has saturation regions")
        self.add_property("other_regions", "boolean", "has other regions")
        self.add_property("galaxy_segments", "boolean", "has galaxy segments")
        self.add_property("star_segments", "boolean", "has star segments")
        self.add_property("other_segments", "boolean", "has other segments")
        self.add_property("attenuation", "real", "attenuation value")

        self.add_property("subtracted", "boolean", "subtraction succesful")
        self.add_property("mean_frame", "real", "mean value")
        self.add_property("median_frame", "real", "median value")
        self.add_property("stddev_frame", "real", "stddev")
        self.add_property("mean_frame_not_clipped", "real", "mean value not clipped")
        self.add_property("median_frame_not_clipped", "real", "median value not clipped")
        self.add_property("stddev_frame_not_clipped", "real", "stddev not clipped")
        self.add_property("mean_sky", "real", "mean sky value")
        self.add_property("median_sky", "real", "median sky value")
        self.add_property("mean_noise", "real", "mean noise")
        self.add_property("mean_subtracted", "real", "mean subtracted value")
        self.add_property("median_subtracted", "real", "median subtracted value")
        self.add_property("stddev_subtracted", "real", "stddev of subtracted frame")

        self.add_property("error_contributions", "string_list", "error contributions")

        self.add_property("original_unit", "string", "original unit")
        self.add_property("conversion_factor", "real", "unit conversion factor")

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class DataPreparer(PreparationComponent):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(DataPreparer, self).__init__(*args, **kwargs)

        # The extinction calculator
        self.extinction = None

        # The DustPedia properties
        self.properties = None

        # The DustPedia photometry object
        self.photometry = None

        # The path lists
        self.initialized_paths = NamedFileList()
        self.extracted_paths = NamedFileList()
        self.corrected_paths = NamedFileList()
        self.subtracted_paths = NamedFileList()
        self.with_errormaps_paths = NamedFileList()
        self.result_paths = NamedFileList()

        # The prepared dataset
        self.prepared_dataset = DataSet()

        # The statistics for each
        self.statistics = dict()

        # The remote host (if needed)
        self.remote = None

        # The remote cache path
        self.remote_preparation_path = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function runs the data preparation ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Extract sources
        self.extract_sources()

        # 3. Correct for galactic extinction
        self.correct_for_extinction()

        # 4. Subtract sky
        self.subtract_sky()

        # 5. Calculate the error maps
        self.create_errormaps()

        # 6. If requested, convert the unit
        self.convert_units()

        # 7. Create the prepared dataset
        self.create_dataset()

        # 8. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataPreparer, self).setup(**kwargs)

        # Setup the remote
        if self.config.cache: self.remote = Remote(host_id=self.environment.cache_host_id)

        # Create the cache directory
        if self.config.cache:
            self.remote_preparation_path = fs.join(self.remote.home_directory, self.galaxy_name + "_preparation")
            if not self.remote.is_directory(self.remote_preparation_path): self.remote.create_directory(self.remote_preparation_path)

        # Get paths
        self.get_paths()

        # Load statistics
        self.load_statistics()

        # Create the galactic extinction calculator
        self.extinction = GalacticExtinction(self.galaxy_center)

        # Create the DustPedia properties
        self.properties = DustPediaProperties()

        # Create the DustPedia photometry object
        self.photometry = DustPediaPhotometry()

    # -----------------------------------------------------------------

    @lazyproperty
    def all_initialized_names(self):

        """
        This function ...
        :return: 
        """

        return self.initial_dataset.names

    # -----------------------------------------------------------------

    @property
    def all_initialized_paths(self):

        """
        This function ...
        :return: 
        """

        paths = dict()
        for name in self.all_initialized_names:
            path = self.initial_dataset.paths[name]
            paths[name] = path
        return paths

    # -----------------------------------------------------------------

    @property
    def all_initialized_directories(self):

        """
        This function ...
        :return: 
        """

        paths = dict()
        for name in self.all_initialized_paths:
            paths[name] = fs.directory_of(self.all_initialized_paths[name])
        return paths

    # -----------------------------------------------------------------

    def get_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Collecting and sorting the image paths ...")

        # Loop over all images of the initial dataset
        for name in self.all_initialized_directories:

            # Determine preparation directory for this image
            path = self.all_initialized_directories[name]

            # Check
            check_initialized(name, path)

            # Sort
            label, filepath = sort_image(name, path, rerun=self.config.rerun)
            if label == "initialized": self.initialized_paths.append(name, filepath)
            elif label == "extracted": self.extracted_paths.append(name, filepath)
            elif label == "corrected": self.corrected_paths.append(name, filepath)
            elif label == "subtracted": self.subtracted_paths.append(name, filepath)
            elif label == "with_errors": self.with_errormaps_paths.append(name, filepath)
            elif label == "result": self.result_paths.append(name, filepath)
            else: raise RuntimeError("Invalid answer for label: " + label)

            # Cache results from all previous steps
            if self.config.cache: self.cache_all_before(name, path, label)

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return: 
        """

        return self.remote.host_id

    # -----------------------------------------------------------------

    def cache_all_before(self, name, directory_path, status):

        """
        This function ...
        :param name:
        :param directory_path: 
        :param status: 
        :return: 
        """

        # Debugging
        log.debug("Caching intermediate results before " + status + " of " + name + " image to remote host '" + self.host_id + "' ...")

        # Get filenames that can be cached
        filenames = filenames_before_status(status)

        # Determine the remote directory for this image
        remote_image_directory_path = fs.join(self.remote_preparation_path, name)
        if not self.remote.is_directory(remote_image_directory_path): self.remote.create_directory(remote_image_directory_path)

        # Upload the files one by one, and delete them afterwards from the local
        for filename in filenames:

            # Debugging
            log.debug("Uploading file '" + filename + "' ...")

            # Determine the local filepath
            filepath = fs.join(directory_path, filename)

            # Upload
            self.remote.upload(filepath, remote_image_directory_path, compress=True)

            # Remove the local file
            fs.remove_file(filepath)

    # -----------------------------------------------------------------

    def load_statistics(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the statistics file ...")

        # Loop over all the names of the initial dataset
        for name in self.all_initialized_directories:

            # Determine prep path
            path = self.all_initialized_directories[name]

            # Determine statistics file path
            statistics_path = fs.join(path, statistics_name)

            # Load or initialize the statistics file
            if fs.is_file(statistics_path): statistics = PreparationStatistics.from_file(statistics_path)
            else: statistics = PreparationStatistics()

            # Set path anyway (so that we can do .save() everytime)
            statistics.path = statistics_path

            # Set
            self.statistics[name] = statistics

    # -----------------------------------------------------------------

    def cache(self, name, path):

        """
        This function ...
        :param name:
        :param path: 
        :return: 
        """

        filename = fs.strip_extension(fs.name(path))

        # Debugging
        log.debug("Caching " + filename + " " + name + " image ...")

        # Determine the remote directory for this image
        remote_image_directory_path = fs.join(self.remote_preparation_path, name)
        if not self.remote.is_directory(remote_image_directory_path): self.remote.create_directory(remote_image_directory_path)

        # Upload
        self.remote.upload(path, remote_image_directory_path)

        # Remove the file
        fs.remove_file(path)

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Extracting the sources ...")

        # Loop over the images
        for name in self.initialized_paths.names:

            # Get the current path
            path = self.initialized_paths.pop(name)

            # Get directry path
            directory_path = fs.directory_of(path)

            # Get sources path
            sources_path = fs.join(directory_path, sources_name)

            # Load the image
            image = Image.from_file(path)

            config = dict()
            config["write"] = False
            config["only_foreground"] = True

            # Extract the sources
            extractor = extract_sources(image, config, sources_path)

            # Determine the new path
            new_path = fs.join(directory_path, extracted_name)

            # Save the image
            image.saveto(new_path)

            # Pop the path and add to extracted
            self.extracted_paths.append(name, new_path)

            # Set statistics
            self.statistics[name].nsources = extractor.nsources
            self.statistics[name].ngalaxy_sources = extractor.ngalaxy_sources
            self.statistics[name].nstar_sources = extractor.nstar_sources
            self.statistics[name].nother_sources = extractor.nother_sources
            self.statistics[name].nforeground = extractor.nforeground
            self.statistics[name].nfailed = extractor.nfailed
            self.statistics[name].nsuccess = extractor.nsuccess
            self.statistics[name].nwith_saturation = extractor.nwith_saturation

            self.statistics[name].galaxy_regions = extractor.galaxy_region is not None
            self.statistics[name].star_regions = extractor.star_region is not None
            self.statistics[name].saturation_regions = extractor.saturation_region is not None
            self.statistics[name].other_regions = extractor.other_region is not None

            self.statistics[name].galaxy_segments = extractor.galaxy_segments is not None
            self.statistics[name].star_segments = extractor.star_segments is not None
            self.statistics[name].other_segments = extractor.other_segments is not None

            # Save the statistics
            self.statistics[name].save()

            # Cache the initialized file on the remote host
            if self.config.cache: self.cache(name, path)

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Correcting for extinction ...")

        # Loop over the images
        for name in self.extracted_paths.names:

            # Get the current path
            path = self.extracted_paths.pop(name)

            # Get directry path
            directory_path = fs.directory_of(path)

            # Get the filter
            fltr = parse_filter(name)

            # Get the extinction value
            attenuation = self.extinction.extinction_for_filter(fltr)

            # Load the image
            image = Image.from_file(path)

            # Correct
            # Do the correction
            # Correct all data frames for galactic extinction (primary, poisson error frame, ...)
            image *= 10 ** (0.4 * attenuation)

            # IMPORTANT: SET FLAG
            image.corrected_for_extinction = True

            # Determine the new path
            new_path = fs.join(directory_path, corrected_name)

            # Save the image
            image.saveto(new_path)

            # Pop the path and add to corrected
            self.corrected_paths.append(name, new_path)

            # Add info to the statistics
            self.statistics[name].attenuation = attenuation
            self.statistics[name].save()

            # Cache the extracted image
            if self.config.cache: self.cache(name, path)

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Subtracting the sky ...")

        # Loop over the corrected images
        for name in self.corrected_paths.names:

            # Get the current path
            path = self.corrected_paths.pop(name)

            # Determine the directory path
            directory_path = fs.directory_of(path)

            # Determine sources path
            sources_path = fs.join(directory_path, sources_name)

            # load the image
            image = Image.from_file(path)

            # If the FWHM of the image is undefined, set it now
            if image.fwhm is None:
                log.warning("The FWHM of the " + name + " image is still undefined. Getting the value from the DustPedia data properties ...")
                image.fwhm = self.properties.get_fwhm(image.filter)
                log.warning("Set the value of the FWHM to " + str(image.fwhm))

            config = dict()

            config["write"] = False
            config["plot"] = False

            # Load regions
            principal_shape = get_principal_shape_sky_from_sources_path(sources_path, image.wcs)
            saturation_regions = get_saturation_regions_sky_from_sources_path(sources_path, image.wcs)
            star_regions = get_star_regions_sky_from_sources_path(sources_path, image.wcs)

            # Determine and create the sky path
            sky_path = fs.join(directory_path, "sky")
            if fs.is_directory(sky_path):
                log.warning("There is alread a sky directory present for the " + name + " image: removing ...")
                fs.remove_directory(sky_path)
            fs.create_directory(sky_path)

            # Use DustPedia aperture
            if self.config.dustpedia_aperture: aperture = self.photometry.get_aperture(self.ngc_name)

            # Use the galaxy shape
            else: aperture = principal_shape * self.config.aperture_galaxy_region_factor

            # Set aperture settings
            config["mask"] = dict()
            config["mask"]["annulus_inner_factor"] = self.config.annulus_inner_factor
            config["mask"]["annulus_outer_factor"] = self.config.annulus_outer_factor
            config["mask"]["saturation_expansion_factor"] = self.config.saturation_expansion_factor
            config["mask"]["stars_expansion_factor"] = self.config.stars_expansion_factor

            # Subtract
            # image, sky_path, config, principal_sky_region, saturation_sky_region=None, star_sky_region=None, visualisation_path=None
            try: subtractor = subtract_sky(image, sky_path, config, aperture, saturation_regions, star_regions)
            except RuntimeError:
                log.warning("The " + name + " image could not be sky subtracted")
                subtractor = None

            # Determine the new path
            new_path = fs.join(directory_path, subtracted_name)

            # Save the image
            image.saveto(new_path)

            # Pop the path and add to subtracted
            self.subtracted_paths.append(name, new_path)

            # Set statistics
            if subtractor is not None:

                self.statistics[name].subtracted = True

                mean_sky = np.nanmean(subtractor.sky)
                median_sky = np.nanmedian(subtractor.sky)
                mean_noise = np.nanmean(subtractor.noise)

                # Set statistics
                self.statistics[name].mean_frame = subtractor.mean_frame
                self.statistics[name].median_frame = subtractor.median_frame
                self.statistics[name].stddev_frame = subtractor.stddev_frame
                self.statistics[name].mean_frame_not_clipped = subtractor.mean_frame_not_clipped
                self.statistics[name].median_frame_not_clipped = subtractor.median_frame_not_clipped
                self.statistics[name].stddev_frame_not_clipped = subtractor.stddev_frame_not_clipped

                self.statistics[name].mean_sky = mean_sky
                self.statistics[name].median_sky = median_sky
                self.statistics[name].mean_noise = mean_noise

                self.statistics[name].mean_subtracted = subtractor.mean_subtracted
                self.statistics[name].median_subtracted = subtractor.median_subtracted
                self.statistics[name].stddev_subtracted = subtractor.stddev_subtracted

            # Subtraction failed
            else: self.statistics[name].subtracted = False

            # Save statistics
            self.statistics[name].save()

            # Cache the extinction corrected image
            if self.config.cache: self.cache(name, path)

    # -----------------------------------------------------------------

    def create_errormaps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the error maps ...")

        # Loop over the subtracted images
        for name in self.subtracted_paths.names:

            # Get the path
            path = self.subtracted_paths.pop(name)

            # Get the directory path
            directory_path = fs.directory_of(path)

            # Get the sky directory path
            sky_path = fs.join(directory_path, sky_name)

            # Get the filter
            fltr = parse_filter(name)

            # Load the image
            image = Image.from_file(path)

            # Calculate calibration error frame
            if self.properties.has_calibration_error_magnitude(fltr):
                magnitude = self.properties.get_calibration_error_magnitude(fltr)
                calibration_frame = get_calibration_uncertainty_frame_from_magnitude(image, magnitude)
            else:

                # Calculate calibration errors with percentage
                fraction = self.properties.get_calibration_error_relative(fltr)
                calibration_frame = image.primary * fraction

            # Add the calibration frame
            image.add_frame(calibration_frame, "calibration_errors")

            # Create a list to contain (the squares of) all the individual error contributions so that we can sum these arrays element-wise later
            error_maps = []
            error_contributions = []

            # Add the Poisson errors
            if "poisson_errors" in image.frames:
                error_maps.append(image.frames["poisson_errors"])
                error_contributions.append("poisson")

            # Load noise frame
            noise_frame = get_noise_frame_from_sky_path(sky_path)

            # Add the sky errors
            error_maps.append(noise_frame)
            error_contributions.append("noise")

            # Add the calibration errors
            error_maps.append(image.frames["calibration_errors"])
            error_contributions.append("calibration")

            # Add additional error frames indicated by the user
            #if self.config.error_frame_names is not None:
            #    for error_frame_name in self.config.error_frame_names: error_maps.append(image.frames[error_frame_name])

            # Calculate the final error map
            errors = sum_frames_quadratically(*error_maps)

            # Add the combined errors frame
            image.add_frame(errors, "errors")

            # Determine the new path
            new_path = fs.join(directory_path, with_errors_name)

            # Save the image
            image.saveto(new_path)

            # Pop the path and add to with_errormaps
            self.with_errormaps_paths.append(name, new_path)

            # Set statistics
            self.statistics[name].error_contributions = error_contributions

            # Save statistics
            self.statistics[name].save()

            # Cache the subtracted image
            if self.config.cache: self.cache(name, path)

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Converting the units ...")

        # Loop over the with_errormaps images
        for name in self.with_errormaps_paths.names:

            # Pop the path
            path = self.with_errormaps_paths[name]

            # Load the image
            image = Image.from_file(path)

            # Remove all frames except for the primary and errors
            image.remove_frames_except("primary", "errors")

            # Get original unit
            original_unit = str(image.unit)

            # Convert to Jansky, except for the Ha image
            if name != "Ha": conversion_factor = image.convert_to("Jy")
            else: conversion_factor = None

            # Pop the path and add to result paths
            self.result_paths.append(name, path)

            # Set statistics
            self.statistics[name].original_unit = original_unit

            if conversion_factor is not None: self.statistics[name].conversion_factor = conversion_factor

            # Save statistics
            self.statistics[name].save()

            # Cache the with_errormaps image
            if self.config.cache: self.cache(name, path)

    # -----------------------------------------------------------------

    def create_dataset(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the prepared dataset ...")

        # Loop over the result paths
        for name in self.result_paths.names:

            # Get the path
            path = self.result_paths[name]

            # Add to the dataset
            self.prepared_dataset.add_path(name, path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the result dataset
        self.write_dataset()

    # -----------------------------------------------------------------

    def write_dataset(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the output dataset ...")

        # Write the dataset
        self.prepared_dataset.saveto(self.prepared_dataset_path)

# -----------------------------------------------------------------

def load_sources(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Load the galaxy region
    galaxy_region_path = fs.join(path, "galaxies.reg")
    galaxy_region = PixelRegionList.from_file(galaxy_region_path)

    # Load the star region (if present)
    star_region_path = fs.join(path, "stars.reg")
    star_region = PixelRegionList.from_file(star_region_path) if fs.is_file(star_region_path) else None

    # load the saturation region (if present)
    saturation_region_path = fs.join(path, "saturation.reg")
    saturation_region = PixelRegionList.from_file(saturation_region_path) if fs.is_file(saturation_region_path) else None

    # Load the region of other sources
    other_region_path = fs.join(path, "other_sources.reg")
    other_region = PixelRegionList.from_file(other_region_path) if fs.is_file(other_region_path) else None

    # Load the image with segmentation maps
    segments_path = fs.join(path, "segments.fits")

    # Get the different segmentation frames
    if fs.is_file(segments_path):
        segments = Image.from_file(segments_path, no_filter=True)
        galaxy_segments = segments.frames["galaxies"] if "galaxies" in segments.frames else None
        star_segments = segments.frames["stars"] if "stars" in segments.frames else None
        other_segments = segments.frames["other_sources"] if "other_sources" in segments.frames else None
    else: galaxy_segments = star_segments = other_segments = None

    # Return the regions and segmentation maps
    return galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments

# -----------------------------------------------------------------

def check_initialized(name, path):

    """
    THis function ...
    :param name:
    :param path: 
    :return: 
    """

    # Debugging
    log.debug("Checking the " + name + " image ...")

    initialized_path = fs.join(path, initialized_name)
    sources_path = fs.join(path, sources_name)

    # Look if an initialized image file is present
    if not fs.is_file(initialized_path):
        # log.warning("Initialized image could not be found for " + path)
        # continue
        raise RuntimeError("Initialized " + name + " image could not be found")

    # Look if the 'sources' directory is present
    if not fs.is_directory(sources_path):
        # log.warning("Sources directory could not be found for " + path)
        # continue
        raise RuntimeError("Sources directory could not be found for the " + name + " image")

# -----------------------------------------------------------------

def sort_image(name, path, rerun=None, read_only=False):

    """
    This function ...
    :param name: 
    :param path: 
    :param rerun:
    :param read_only:
    :return: 
    """

    # Debugging
    log.debug("Sorting the " + name + " image ...")

    initialized_path = fs.join(path, initialized_name)

    # Check if the intermediate results have already been produced for this image and saved to the
    # corresponding preparation subdirectory
    extracted_path = fs.join(path, extracted_name)
    corrected_path = fs.join(path, corrected_name)
    subtracted_path = fs.join(path, subtracted_name)
    with_errors_path = fs.join(path, with_errors_name)
    result_path = fs.join(path, result_name)

    # Sky directory path
    sky_path = fs.join(path, sky_name)

    # Run through the different checks
    if check_result(name, result_path, read_only=read_only): return "result", result_path
    elif check_with_errors(name, with_errors_path, rerun=rerun, read_only=read_only): return "with_errors", with_errors_path
    elif check_subtracted(name, subtracted_path, sky_path, rerun=rerun, read_only=read_only): return "subtracted", subtracted_path
    elif check_extinction_corrected(name, corrected_path, rerun=rerun, read_only=read_only): return "corrected", corrected_path
    elif check_extracted(name, extracted_path, rerun=rerun, read_only=read_only): return "extracted", extracted_path
    else: return "initialized", initialized_path

# -----------------------------------------------------------------

def check_result(name, path, rerun=None, read_only=False):

    """
    This function ...
    :param name: 
    :param path:
    :param rerun: 
    :param read_only:
    :return: 
    """

    if fs.is_file(path):

        if rerun is not None and rerun in steps_before_and_including("units"):

            if read_only: raise RuntimeError("Cannot rerun when read_only=True")
            fs.remove_file(path)
            return False

        else: return True

    else: return False

# -----------------------------------------------------------------

def check_with_errors(name, path, rerun=None, read_only=False):

    """
    This function ...
    :param name: 
    :param rerun: 
    :param read_only:
    :return: 
    """

    if fs.is_file(path):

        if rerun is not None and rerun in steps_before_and_including("errormaps"):

            if read_only: raise RuntimeError("Cannot rerun when read_only=True")
            fs.remove_file(path)
            return False

        else: return True

    else: return False

# -----------------------------------------------------------------

def check_subtracted(name, subtracted_path, sky_path, rerun=None, read_only=False):

    """
    This fucntion ...
    :param name:
    :param subtracted_path: 
    :param sky_path: 
    :param rerun:
    :param read_only:
    :return: 
    """

    # Subtracted file is present
    if fs.is_file(subtracted_path):

        # NO: because it's OK that for some images the sky subtraction failed and therefore they don't have the sky directory
        # Check whether the sky directory is present
        #if not fs.is_directory(sky_path): #raise IOError("The sky subtraction output directory is not present for the '" + name + "' image")
            #log.warning("The sky subtraction output directory is not present for the " + name + " image")
            #log.warning("Removing the subtracted image file ...")
            #fs.remove_file(subtracted_path)

        #elif fs.is_empty(sky_path):
        if fs.is_directory(sky_path) and fs.is_empty(sky_path):

            log.warning("The sky subtraction output directory is empty for the " + name + " image")
            log.warning("Removing the subtracted image file and the empty directory ...")

            if read_only: raise RuntimeError("Cannot remove when read_only=True")

            fs.remove_file(subtracted_path)
            fs.remove_directory(sky_path)

            return False

        elif rerun is not None and rerun in steps_before_and_including("subtraction"):

            if read_only: raise RuntimeError("Cannot remove when read_only=True")

            if fs.is_directory(sky_path): fs.remove_directory(sky_path)
            fs.remove_file(subtracted_path)

            return False

        else: return True

    else: return False

# -----------------------------------------------------------------

def check_extinction_corrected(name, path, rerun=None, read_only=False):

    """
    This function ...
    :param name: 
    :param path: 
    :param rerun:
    :param read_only:
    :return: 
    """

    if fs.is_file(path):

        if rerun is not None and rerun in steps_before_and_including("extinction"):

            if read_only: raise RuntimeError("Cannot remove when read_only=True")
            fs.remove_file(path)
            return False

        else: return True

    else: return False

# -----------------------------------------------------------------

def check_extracted(name, path, rerun=None, read_only=False):

    """
    This function ...
    :param name: 
    :param path: 
    :param rerun: 
    :param read_only:
    :return: 
    """

    if fs.is_file(path):

        if rerun is not None and rerun in steps_before_and_including("extraction"):

            if read_only: raise RuntimeError("Cannot remove when read_only=True")
            fs.remove_file(path)
            return False

        else: return True

    else: return False

# -----------------------------------------------------------------

def extract_sources(image, config, sources_path, visualisation_path=None):

    """
    This function ...
    :param image:
    :param config:
    :param sources_path:
    :param visualisation_path:
    :return:
    """

    # Inform the user
    log.info("Extracting the sources ...")

    # Create an animation to show the result of this step
    #if visualisation_path is not None:
    #    # Create an animation
    #    animation = ImageBlinkAnimation()
    #    animation.add_image(image.primary)
    #    # Create an animation to show the progress of the SourceExtractor
    #    source_extractor_animation = SourceExtractionAnimation(image.primary)
    #else: animation = source_extractor_animation = None
    source_extractor_animation = None

    # Create the source extractor
    extractor = SourceExtractor(config)

    # Set configuration settings
    extractor.config.input = sources_path

    # Run the extraction
    special_region = None
    extractor.run(frame=image.primary, animation=source_extractor_animation, special_region=special_region)

    # Get the saturation region
    #saturation_region = extractor.saturation_region

    # Write the animation
    #if visualisation_path is not None:
        # Determine the path to the animation
        #path = fs.join(visualisation_path, time.unique_name(image.name + "_sourceextraction") + ".gif")
        # Save the animation
        #source_extractor_animation.saveto(path)
        # ...
        #animation.add_image(image.primary)
        # Determine the path to the animation
        #path = fs.join(visualisation_path, time.unique_name(image.name + "_imagepreparation_extractsources") + ".gif")
        # Save the animation
        #animation.saveto(path)

    # Add the sources mask to the image
    image.add_mask(extractor.mask, sources_name)

    # Get the principal shape in sky coordinates
    #principal_shape_sky = extractor.principal_shape.to_sky(image.wcs)

    # Get the saturation region in sky coordinates
    #saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

    # IMPORTANT: SET FLAG
    image.source_extracted = True

    # Return
    #return principal_shape_sky, saturation_region_sky

    return extractor

# -----------------------------------------------------------------

def subtract_sky(image, sky_path, config, principal_sky_region, saturation_sky_region=None, star_sky_region=None, visualisation_path=None):

    """
    This function ...
    :param image:
    :param sky_path:
    :param config:
    :param principal_sky_region:
    :param saturation_sky_region:
    :param star_sky_region:
    :param visualisation_path:
    :return:
    """

    # Inform the user
    log.info("Subtracting the sky ...")

    # SET FLAG TO FALSE
    image.sky_subtracted = False

    # Create an animation to show the result of this step
    #if visualisation_path is not None:
    #    # Create an animation to show the result of the sky subtraction step
    #    animation = ImageBlinkAnimation()
    #    animation.add_image(image.primary)
    #    # Create an animation to show the progress of the SkySubtractor
    #    skysubtractor_animation = Animation()
    #else: animation = skysubtractor_animation = None
    skysubtractor_animation = None

    # Convert the principal ellipse in sky coordinates into pixel coordinates
    principal_shape = principal_sky_region.to_pixel(image.wcs)

    # Convert the saturation region in sky coordinates into pixel coordinates
    if saturation_sky_region is not None: saturation_region = saturation_sky_region.to_pixel(image.wcs)
    else: saturation_region = None

    # Convert the star region in sky coordinate into pixel coordinates
    if star_sky_region is not None: star_region = star_sky_region.to_pixel(image.wcs)
    else: star_region = None

    # Create the 'extra' mask (bad and padded pixels)
    extra_mask = None
    if "bad" in image.masks:
        # Combine with padded mask or just use bad mask
        if "padded" in image.masks:
            extra_mask = image.masks.bad + image.masks.padded
        else: extra_mask = image.masks.bad
    elif "padded" in image.masks: extra_mask = image.masks.padded  # just use padded mask

    # Create the sky subtractor
    sky_subtractor = SkySubtractor(config)

    # Run the sky subtractor
    sky_subtractor.run(frame=image.primary, principal_shape=principal_shape, sources_mask=image.masks.sources,
                       extra_mask=extra_mask, saturation_region=saturation_region, star_region=star_region, animation=skysubtractor_animation)

    # Set the subtracted frame as the primary frame
    image.replace_frame("primary", sky_subtractor.subtracted)

    # Add the sky frame to the image
    image.add_frame(sky_subtractor.sky_frame, "sky")

    # Add the mask that is used for the sky estimation
    image.add_mask(sky_subtractor.mask, "sky")

    # Add the sky noise frame
    image.add_frame(sky_subtractor.noise_frame, "sky_errors")

    # Write sky annuli maps if requested
    #if sky_apertures_path is not None:

    # Write the sky region
    region_path = fs.join(sky_path, "annulus.reg")
    sky_subtractor.region.saveto(region_path)

    # Write the apertures frame
    apertures_frame_path = fs.join(sky_path, "apertures.fits")
    if sky_subtractor.apertures_frame is not None: sky_subtractor.apertures_frame.saveto(apertures_frame_path)

    # Write the apertures mean frame
    apertures_mean_path = fs.join(sky_path, "apertures_values.fits")
    if sky_subtractor.apertures_values_frame is not None: sky_subtractor.apertures_values_frame.saveto(apertures_mean_path)

    # Write the apertures noise frame
    apertures_noise_path = fs.join(sky_path, "apertures_noise.fits")
    if sky_subtractor.apertures_noise_frame is not None: sky_subtractor.apertures_noise_frame.saveto(apertures_noise_path)

    # Write photutils results
    background_mesh_path = fs.join(sky_path, "background_mesh.fits")
    if sky_subtractor.phot_background_mesh is not None: sky_subtractor.phot_background_mesh.saveto(background_mesh_path)

    background_rms_mesh_path = fs.join(sky_path, "background_rms_mesh.fits")
    if sky_subtractor.phot_background_rms_mesh is not None: sky_subtractor.phot_background_rms_mesh.saveto(background_rms_mesh_path)

    estimated_sky_path = fs.join(sky_path, "estimated_sky.fits")
    if sky_subtractor.phot_sky is not None: sky_subtractor.phot_sky.saveto(estimated_sky_path)

    estimated_sky_rms_path = fs.join(sky_path, "estimated_sky_rms.fits")
    if sky_subtractor.phot_rms is not None: sky_subtractor.phot_rms.saveto(estimated_sky_rms_path)

    # WRITE THE MASK
    #mask_path = fs.join(mask_path, )

    # Write boundaries for cutting out a piece of the frame for photutils
    if sky_subtractor.phot_boundaries is not None:
        phot_boundaries_path = fs.join(sky_path, "boundaries.dat")
        write_dict(sky_subtractor.phot_boundaries, phot_boundaries_path)

    # WRITE MASK CONTRIBUTIONS
    if sky_subtractor.sources_mask is not None:
        sources_mask_path = fs.join(sky_path, "sources_mask.fits")
        mask = Mask(sky_subtractor.sources_mask)
        mask.saveto(sources_mask_path)

    if sky_subtractor.outside_mask is not None:
        outside_mask_path = fs.join(sky_path, "outside_mask.fits")
        mask = Mask(sky_subtractor.outside_mask)
        mask.saveto(outside_mask_path)

    if sky_subtractor.principal_mask is not None:
        principal_mask_path = fs.join(sky_path, "principal_mask.fits")
        mask = Mask(sky_subtractor.principal_mask)
        mask.saveto(principal_mask_path)

    if sky_subtractor.saturation_mask is not None:
        saturation_mask_path = fs.join(sky_path, "saturation_mask.fits")
        mask = Mask(sky_subtractor.saturation_mask)
        mask.saveto(saturation_mask_path)

    if sky_subtractor.stars_mask is not None:
        stars_mask_path = fs.join(sky_path, "stars_mask.fits")
        mask = Mask(sky_subtractor.stars_mask)
        mask.saveto(stars_mask_path)

    if sky_subtractor.extra_mask is not None:
        extra_mask_path = fs.join(sky_path, "extra_mask.fits")
        mask = Mask(sky_subtractor.extra_mask)
        mask.saveto(extra_mask_path)

    # Write properties

    # Write the animation
    #if visualisation_path is not None:
        # Determine the path to the animation
        #path = fs.join(visualisation_path, time.unique_name(image.name + "_skysubtraction") + ".gif")
        # Save the animation
        #skysubtractor_animation.saveto(path)
        # ...
        #animation.add_image(image.primary)
        # Determine the path to the animation
        #path = fs.join(visualisation_path, time.unique_name(image.name + "_imagepreparation_subtractsky") + ".gif")
        # Save the animation
        #animation.saveto(path)

    # NEW: WRITE THE NOISE FRAME
    noise_frame_path = fs.join(sky_path, "noise.fits")
    sky_subtractor.noise_frame.saveto(noise_frame_path)

    # IMPORTANT: SET FLAG
    image.sky_subtracted = True

    # Return the noise frame
    #return sky_subtractor.noise_frame

    return sky_subtractor

# -----------------------------------------------------------------

def get_calibration_uncertainty_frame_from_magnitude(image, calibration_magn):

    """
    This fucntion ...
    :param image:
    :param calibration_magn:
    :return: 
    """

    # Convert the frame into AB magnitudes
    invalid = oldMask.is_zero_or_less(image.primary)
    ab_frame = unitconversion.jansky_to_ab(image.primary)
    # Set infinites to zero
    ab_frame[invalid] = 0.0

    # -----------------------------------------------------------------

    # The calibration uncertainty in AB magnitude
    #mag_error = calibration_magn.value
    mag_error = calibration_magn

    # a = image[mag] - mag_error
    a = ab_frame - mag_error

    # b = image[mag] + mag_error
    b = ab_frame + mag_error

    # Convert a and b to Jy
    a = unitconversion.ab_mag_zero_point.to("Jy").value * np.power(10.0, -2. / 5. * a)
    b = unitconversion.ab_mag_zero_point.to("Jy").value * np.power(10.0, -2. / 5. * b)

    # c = a[Jy] - image[Jy]
    # c = a - jansky_frame
    c = a - image.primary

    # d = image[Jy] - b[Jy]
    # d = jansky_frame - b
    d = image.primary - b

    # ----------------------------------------------------------------- BELOW: if frame was not already in Jy

    # Calibration errors = max(c, d)
    # calibration_errors = np.maximum(c, d)  # element-wise maxima
    # calibration_errors[invalid] = 0.0 # set zero where AB magnitude could not be calculated

    # relative_calibration_errors = calibration_errors / jansky_frame
    # relative_calibration_errors[invalid] = 0.0

    # Check that there are no infinities or nans in the result
    # assert not np.any(np.isinf(relative_calibration_errors)) and not np.any(np.isnan(relative_calibration_errors))

    # The actual calibration errors in the same unit as the data
    # calibration_frame = self.image.frames.primary * relative_calibration_errors

    # -----------------------------------------------------------------

    calibration_frame = np.maximum(c, d)  # element-wise maxima
    calibration_frame[invalid] = 0.0  # set zero where AB magnitude could not be calculated

    new_invalid = oldMask.is_nan(calibration_frame) + oldMask.is_inf(calibration_frame)
    calibration_frame[new_invalid] = 0.0

    # Check that there are no infinities or nans in the result
    assert not np.any(np.isinf(calibration_frame)) and not np.any(np.isnan(calibration_frame))

    # Make frame from numpy array
    calibration_frame = Frame(calibration_frame)

    # Return the frame
    return calibration_frame

#-----------------------------------------------------------------

def get_principal_shape(galaxy_region):

    """
    This function ...
    :param galaxy_region:
    :return:
    """

    if galaxy_region is None: return None

    largest_shape = None

    # Loop over all the shapes in the galaxy region
    for shape in galaxy_region:

        # Skip single coordinates
        if isinstance(shape, PixelCoordinate): continue

        if shape.label is not None and "principal" in shape.label: return shape
        if "text" in shape.meta and "principal" in shape.meta["text"]: return shape

        if not isinstance(shape, PixelEllipseRegion) and not isinstance(shape, PixelPointRegion): return shape

        semimajor_axis_length = shape.semimajor
        if largest_shape is None or semimajor_axis_length > largest_shape.semimajor: largest_shape = shape

    # Return the largest shape
    return largest_shape

# -----------------------------------------------------------------

def get_principal_shape_sky(galaxy_region, wcs):

    """
    This function ...
    :param galaxy_region: 
    :param wcs: 
    :return: 
    """

    shape_pixel = get_principal_shape(galaxy_region)
    principal_shape_sky = shape_pixel.to_sky(wcs)
    return principal_shape_sky

# -----------------------------------------------------------------

def get_principal_shape_sky_from_sources_path(sources_path, wcs):

    """
    This function ...
    :param sources_path: 
    :param wcs: 
    :return: 
    """

    path = fs.join(sources_path, "galaxies.reg")
    regions = PixelRegionList.from_file(path)
    return get_principal_shape_sky(regions, wcs)

# -----------------------------------------------------------------

def get_saturation_regions_sky_from_sources_path(sources_path, wcs):

    """
    This function ...
    :param sources_path: 
    :param wcs: 
    :return: 
    """

    # Determine path
    path = fs.join(sources_path, "saturation.reg")

    if not fs.is_file(path): return None

    regions = PixelRegionList.from_file(path)

    # Get the saturation region in sky coordinates
    saturation_regions_sky = regions.to_sky(wcs)

    # Return the saturation regions in sky coordinates
    return saturation_regions_sky

# -----------------------------------------------------------------

def get_star_regions_sky_from_sources_path(sources_path, wcs):

    """
    This funciton ...
    :param sources_path: 
    :param wcs: 
    :return: 
    """

    # Determine path
    path = fs.join(sources_path, "stars.reg")

    if not fs.is_file(path): return None

    regions = PixelRegionList.from_file(path)

    # Get in sky coordinates
    star_regions_sky = regions.to_sky(wcs)

    # Return
    return star_regions_sky

# -----------------------------------------------------------------

def get_noise_frame_from_sky_path(sky_path):

    """
    This function ...
    :param sky_path: 
    :return: 
    """

    path = fs.join(sky_path, "noise.fits")
    return Frame.from_file(path)

# -----------------------------------------------------------------
