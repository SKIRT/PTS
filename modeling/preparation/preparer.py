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
from ...core.tools import time
from ...magic.sources.extractor import SourceExtractor
from ...magic.sky.skysubtractor import SkySubtractor
from ...dustpedia.core.properties import DustPediaProperties
from ...magic.core.mask import Mask
from ...magic.basics.mask import Mask as oldMask
from ...magic.core.frame import Frame, sum_frames_quadratically
from . import unitconversion
from ...magic.region.point import PixelPointRegion
from ...magic.basics.coordinate import PixelCoordinate
from ...magic.region.ellipse import PixelEllipseRegion

# -----------------------------------------------------------------

initialized_name = "initialized.fits"
extracted_name = "extracted.fits"
corrected_name = "corrected_for_extinction.fits"
subtracted_name = "sky_subtracted.fits"
with_errors_name = "with_errormaps.fits"
result_name = "result.fits"

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

        # The path lists
        self.initialized_paths = NamedFileList()
        self.extracted_paths = NamedFileList()
        self.corrected_paths = NamedFileList()
        self.subtracted_paths = NamedFileList()
        self.with_errormaps_paths = NamedFileList()
        self.result_paths = NamedFileList()

        # The FWHM of the reference image
        self.reference_fwhm = None

        # The Aniano kernels service
        self.aniano = None

        # The prepared frames
        #self.frames = None

        # The prepared dataset
        self.prepared_dataset = DataSet()

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

        # Create the prepared dataset
        self.create_dataset()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataPreparer, self).setup(**kwargs)

        # Set options for the image preparer
        #self.set_preparer_options()

        # Setup the remote PTS launcher
        #if self.config.remote is not None: self.launcher.setup(self.config.remote)
        #else: self.preparer = ImagePreparer(self.preparer_config)

        # Get paths
        self.get_paths()

        # Create the galactic extinction calculator
        self.extinction = GalacticExtinction(self.galaxy_center)

        # Create the DustPedia properties
        self.properties = DustPediaProperties()

    # -----------------------------------------------------------------

    def get_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Collecting and sorting the image paths ...")

        # Loop over all images of the initial dataset
        for name in self.initial_dataset.names:

            # Get path of initial image
            image_path = self.initial_dataset.paths[name]

            # Determine preparation directory for this image
            path = fs.directory_of(image_path)

            # Check
            check_initialized(name, path)

            # Sort
            label, filepath = sort_image(name, path)
            if label == "initialized": self.initialized_paths.append(name, filepath)
            elif label == "extracted": self.extracted_paths.append(name, filepath)
            elif label == "corrected": self.corrected_paths.append(name, filepath)
            elif label == "subtracted": self.subtracted_paths.append(name, filepath)
            elif label == "with_errors": self.with_errormaps_paths.append(name, filepath)
            elif label == "result": self.result_paths.append(name, filepath)
            else: raise RuntimeError("Invalid answer for label: " + label)

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return: 
        """

        # Infomr the user
        log.info("Extracting the sources ...")

        # Loop over the images
        for name in self.initialized_paths.names:

            # Get the current path
            path = self.initialized_paths.pop(name)

            # Get directry path
            directory_path = fs.directory_of(path)

            # Get sources path
            sources_path = fs.join(directory_path, "sources")

            # Load the image
            image = Image.from_file(path)

            config = dict()
            config["only_foreground"] = True

            # Extract the sources
            extract_sources(image, config, sources_path)

            # Determine the new path
            new_path = fs.join(directory_path, extracted_name)

            # Save the image
            image.saveto(new_path)

            # Pop the path and add to extracted
            self.extracted_paths.append(name, new_path)

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
            sources_path = fs.join(directory_path, "sources")

            # load the image
            image = Image.from_file(path)

            config = dict()

            principal_shape = get_principal_shape_sky_from_sources_path(sources_path, image.wcs)
            saturation_regions = get_saturation_regions_sky_from_sources_path(sources_path, image.wcs)

            # Subtract
            subtract_sky(image, config, principal_shape, saturation_regions)

            # Determine the new path
            new_path = fs.join(directory_path, subtracted_name)

            # Save the image
            image.saveto(new_path)

            # Pop the path and add to subtracted
            self.subtracted_paths.append(name, new_path)

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
            sky_path = fs.join(directory_path, "sky")

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

            # Add the Poisson errors
            if "poisson_errors" in image.frames: error_maps.append(image.frames["poisson_errors"])

            # Load noise frame
            noise_frame = get_noise_frame_from_sky_path(sky_path)

            # Add the sky errors
            error_maps.append(noise_frame)

            # Add the calibration errors
            error_maps.append(image.frames["calibration_errors"])

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

            # Convert to Jansky, except for the Ha image
            if name != "Ha": image.convert_to("Jy")

            # Pop the path and add to result paths
            self.result_paths.append(name, path)

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

    initialized_path = fs.join(path, "initialized.fits")
    sources_path = fs.join(path, "sources")

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

def sort_image(name, path):

    """
    This function ...
    :param name: 
    :param path: 
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
    #converted_path = fs.join(path, "converted_units.fits")
    result_path = fs.join(path, result_name)

    # Sky directory path
    sky_path = fs.join(path, "sky")

    # -----------------------------------------------------------------

    # STEPS:

    # Source extraction
    # Correct for galactic extinction
    # Sky SUBTRACTION
    # CALCULATE ERROR MAPS
    # CONVERT UNITS

    # -----------------------------------------------------------------

    # ALREADY COMPLETELY PREPARED

    # Check if a prepared image is already present
    if fs.is_file(result_path):

        #self.result_paths.append(name, result_path)
        return "result", result_path

    # -----------------------------------------------------------------

    # ALREDAY WITH ERROR MAPS
    if fs.is_file(with_errors_path):

        return "with_errors", with_errors_path

    # ALREADY SKY-SUBTRACTED
    if fs.is_file(subtracted_path):

        # Check whether the sky directory is present
        if not fs.is_directory(sky_path): raise IOError("The sky subtraction output directory is not present for the '" + name + "' image")

        # Add the path of the sky-subtracted image
        #self.preparation_dataset.add_path(prep_name, subtracted_path)

        # Check whether keywords are set to True in image header ?

        return "subtracted", subtracted_path

    # -----------------------------------------------------------------

    # ALREADY EXTINCTION CORRECTED

    # Check if the extinction-corrected image is present
    elif fs.is_file(corrected_path):

        # Add the path of the extinction-corrected image
        #self.preparation_dataset.add_path(prep_name, corrected_path)

        return "corrected", corrected_path

    # ALREADY SOURCE-EXTRACTED

    # Check if the source-extracted image is present
    elif fs.is_file(extracted_path):

        # Add the path of the source-extracted image
        #self.preparation_dataset.add_path(prep_name, extracted_path)

        return "extracted", extracted_path

    # -----------------------------------------------------------------

    # NO STEPS PERFORMED YET, START FROM INITIALIZED IMAGE

    else:

        # Add the path to the initialized image to the dataset
        #self.preparation_dataset.add_path(prep_name, image_path)

        return "initialized", initialized_path

        # -----------------------------------------------------------------

    # If all images have already been prepared
    #if len(self.preparation_dataset) == 0: log.success("All images are already prepared")

# -----------------------------------------------------------------

def extract_sources(image, config, sources_path, visualisation_path=None):

    """
    This function ...
    :param image:
    :param config:
    :param sources_path:
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
    image.add_mask(extractor.mask, "sources")

    # Get the principal shape in sky coordinates
    #principal_shape_sky = extractor.principal_shape.to_sky(image.wcs)

    # Get the saturation region in sky coordinates
    #saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

    # IMPORTANT: SET FLAG
    image.source_extracted = True

    # Return
    #return principal_shape_sky, saturation_region_sky

# -----------------------------------------------------------------

def subtract_sky(image, sky_path, config, principal_sky_region, saturation_sky_region=None, visualisation_path=None):

    """
    This function ...
    :param image:
    :param config:
    :param principal_sky_region:
    :param saturation_sky_region:
    :return:
    """

    # Inform the user
    log.info("Subtracting the sky ...")

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
                       extra_mask=extra_mask, saturation_region=saturation_region, animation=skysubtractor_animation)

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
    sky_subtractor.apertures_frame.saveto(apertures_frame_path)

    # Write the apertures mean frame
    apertures_mean_path = fs.join(sky_path, "apertures_values.fits")
    sky_subtractor.apertures_values_frame.saveto(apertures_mean_path)

    # Write the apertures noise frame
    apertures_noise_path = fs.join(sky_path, "apertures_noise.fits")
    sky_subtractor.apertures_noise_frame.saveto(apertures_noise_path)

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

def get_noise_frame_from_sky_path(sky_path):

    """
    This function ...
    :param sky_path: 
    :return: 
    """

    path = fs.join(sky_path, "noise.fits")
    return Frame.from_file(path)

# -----------------------------------------------------------------
