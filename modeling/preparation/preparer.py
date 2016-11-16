#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.datapreparation Contains the DataPreparer class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...magic.region.list import PixelRegionList
from .component import PreparationComponent
from ...magic.prepare.preparer import ImagePreparer
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.region import tools as regions
from ...magic.misc.kernels import AnianoKernels, aniano_names
from ...magic.misc.calibration import CalibrationError
from ...magic.misc.extinction import GalacticExtinction
from ...core.basics.filter import Filter

# -----------------------------------------------------------------

class DataPreparer(PreparationComponent):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(DataPreparer, self).__init__(config)

        # -- Attributes --

        # The paths to the initialized images
        self.paths = []

        # Information about the images
        self.attenuations = dict()

        # The FWHM of the reference image
        self.reference_fwhm = None

        # The coordinate of the center of the galaxy
        self.center_coordinate = None

        # The Aniano kernels service
        self.aniano = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new DataPreparer instance
        preparer = cls(arguments.config)

        # Whether to write the results of intermediate steps
        preparer.config.preparation.write_steps = arguments.steps

        # Set the reference image
        if arguments.reference is not None: preparer.reference_image = arguments.reference

        # Set the modeling path
        preparer.config.path = arguments.path

        # A single image can be specified so the preparation is only run with that image
        preparer.config.single_image = arguments.image

        # Make visualisations
        preparer.config.visualise = arguments.visualise

        # Return the new instance
        return preparer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function runs the data preparation ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Check which images can be prepared
        self.check_images()

        # If all images have already been prepared, break
        if len(self.paths) == 0:
            log.success("All images are already prepared")
            return

        # 3. Get attenuations
        self.get_attenuations()

        # 4. Prepare the images
        self.prepare_images()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # -- Children --

        # Create the preparation object
        self.add_child("image_preparer", ImagePreparer, self.config.preparation)

        # -- Setup of the base class --

        # Call the setup function of the base class
        super(DataPreparer, self).setup()

        # -- Fixed properties for the image preparer (valid for all target images)

        # Set the path to the reference image for the rebinning
        reference_path = fs.join(self.prep_paths[self.config.reference_image], "initialized.fits")

        # Set the path of the rebinning reference path and the kernel image
        self.image_preparer.config.rebinning.rebin_to = reference_path

        # Get the FWHM of the reference image
        reference_frame = Frame.from_file(reference_path)
        self.reference_fwhm = reference_frame.fwhm

        # Get the center coordinate of the galaxy
        self.center_coordinate = reference_frame.coordinate_range[0]

        # Create the Aniano kernels service
        self.aniano = AnianoKernels()

    # -----------------------------------------------------------------

    def check_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the initialized images ...")

        # Loop over all subdirectories of the preparation directory
        for path in fs.directories_in_path(self.prep_path):

            # Debugging
            log.debug("Opening " + path + " ...")

            # Look if an initialized image file is present
            image_path = fs.join(path, "initialized.fits")
            if not fs.is_file(image_path):

                log.warning("Initialized image could not be found for " + path)
                continue

            # Look if the 'sources' directory is present
            sources_path = fs.join(path, "sources")
            if not fs.is_directory(sources_path):

                log.warning("Sources directory could not be found for " + path)
                continue

            # Check if a prepared image is already present
            result_path = fs.join(path, "result.fits")
            if fs.is_file(result_path): continue

            # Add the path to the initialized image to the list
            self.paths.append(image_path)

    # -----------------------------------------------------------------

    def get_attenuations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galactic extinction values for the different images ...")

        # Create the galactic extinction calculator
        extinction = GalacticExtinction(self.center_coordinate)

        # Loop over all image paths
        for image_path in self.paths:

            # Get the filter name
            filter_name = fs.name(fs.directory_of(image_path))

            # Create a filter instance
            fltr = Filter.from_string(filter_name)

            # Get the exintinction
            self.attenuations[filter_name] = extinction.extinction_for_filter(fltr)

    # -----------------------------------------------------------------

    def prepare_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the images ...")

        # Loop over the image paths
        for image_path in self.paths:

            # Get the directory containing this image = the output path for that image
            output_path = fs.directory_of(image_path)

            # Get the directory containing the output from the SourceFinder
            sources_path = fs.join(output_path, "sources")

            # Get the name
            name = fs.name(output_path)

            # Inform the user
            log.info("Starting preparation of " + name + " image ...")

            # -----------------------------------------------------------------

            # Load the initialized image
            image = Image.from_file(image_path)
            image.name = name

            # Load the output of the source finder
            galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments = load_sources(sources_path)

            # -----------------------------------------------------------------

            # Reset all flags to True
            self.enable_all_preparation_steps()

            # Set options for the ImagePreparation class
            self.set_preparation_options(image, output_path)

            # Check if the intermediate results have already been produced for this image and saved to the
            # corresponding preparation subdirectory
            extracted_path = fs.join(output_path, "extracted.fits")
            corrected_path = fs.join(output_path, "corrected_for_extinction.fits")
            converted_path = fs.join(output_path, "converted_unit.fits")
            convolved_path = fs.join(output_path, "convolved.fits")
            rebinned_path = fs.join(output_path, "rebinned.fits")
            subtracted_path = fs.join(output_path, "sky_subtracted.fits")

            ## CURRENT ORDER OF STEPS IN IMAGEPREPARER:

            # 1. Setup

            # 2. Extract stars and galaxies from the image
            # extract_sources

            # 3. If requested, correct for galactic extinction
            # correct_for_extinction

            # 4. If requested, convert the unit
            # convert_unit

            # 5. If requested, convolve
            # convolve

            # 6. If requested, rebin
            # rebin

            # 7. If requested, subtract the sky
            # subtract_sky

            # 8. Calculate the calibration uncertainties
            # calculate_calibration_uncertainties

            # 9. If requested, set the uncertainties
            # set_uncertainties

            ##

            # Check if the sky-subtracted image is present
            if fs.is_file(subtracted_path):

                # Disable all steps preceeding and including the sky subtraction
                self.image_preparer.config.extract_sources = False
                self.image_preparer.config.correct_for_extinction = False
                self.image_preparer.config.convert_unit = False
                self.image_preparer.config.convolve = False
                self.image_preparer.config.rebin = False
                self.image_preparer.config.subtract_sky = False

                # Set the principal ellipse and saturation region in sky coordinates
                self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(self.image.wcs)
                self.image_preparer.saturation_region_sky = saturation_region.to_sky(self.image.wcs) if saturation_region is not None else None

                # Load the sky-subtracted image
                image = Image.from_file(subtracted_path)
                image.name = name

            # Check if the rebinned image is present
            elif fs.is_file(rebinned_path):

                # Disable all steps preceeding and including the rebinning
                self.image_preparer.config.extract_sources = False
                self.image_preparer.config.correct_for_extinction = False
                self.image_preparer.config.convert_unit = False
                self.image_preparer.config.convolve = False
                self.image_preparer.config.rebin = False

                # Set the principal ellipse and saturation region in sky coordinates
                self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(image.wcs)
                self.image_preparer.saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

                # Load the rebinned image
                image = Image.from_file(rebinned_path)
                image.name = name

            # Check if the convolved image is present
            elif fs.is_file(convolved_path):

                # Disable all steps preceeding and including the convolution
                self.image_preparer.config.extract_sources = False
                self.image_preparer.config.correct_for_extinction = False
                self.image_preparer.config.convert_unit = False
                self.image_preparer.config.convolve = False

                # Set the principal ellipse and saturation region in sky coordinates
                self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(image.wcs)
                self.image_preparer.saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

                # Load the convolved image
                image = Image.from_file(convolved_path)
                image.name = name

            # Check if the converted image is present
            elif fs.is_file(converted_path):

                # Disable all steps preceeding and including the unit conversion
                self.image_preparer.config.extract_sources = False
                self.image_preparer.config.correct_for_extinction = False
                self.image_preparer.config.convert_unit = False

                # Set the principal ellipse and saturation region in sky coordinates
                self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(image.wcs)
                self.image_preparer.saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

                # Load the converted image
                image = Image.from_file(converted_path)
                image.name = name

            # Check if the extinction-corrected image is present
            elif fs.is_file(corrected_path):

                # Disable all steps preceeding and including the correction for extinction
                self.image_preparer.config.extract_sources = False
                self.image_preparer.config.correct_for_extinction = False

                # Set the principal ellipse and saturation region in sky coordinates
                self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(image.wcs)
                self.image_preparer.saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

                # Load the extinction-corrected image
                image = Image.from_file(corrected_path)
                image.name = name

            # Check if the source-extracted image is present
            elif fs.is_file(extracted_path):

                # Disable all steps preceeding and including the source extraction
                self.image_preparer.config.extract_sources = False

                # Set the principal ellipse and saturation region in sky coordinates
                self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(image.wcs)
                self.image_preparer.saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

                # Load the extracted image
                image = Image.from_file(extracted_path)
                image.name = name

            # -----------------------------------------------------------------

            # Write out sky annuli frames
            sky_path = fs.join(output_path, "sky")
            if not fs.is_directory(sky_path): fs.create_directory(sky_path)
            self.image_preparer.config.write_sky_apertures = True
            self.image_preparer.config.sky_apertures_path = sky_path

            # Set the visualisation path for the image preparer
            visualisation_path = self.visualisation_path if self.config.visualise else None

            # -----------------------------------------------------------------

            # Run the image preparation
            self.image_preparer.run(image, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments, visualisation_path)

            # -----------------------------------------------------------------

            # Inform the user
            log.success("Preparation of " + name + " image finished")

            # Clear the image preparer
            self.image_preparer.clear()

    # -----------------------------------------------------------------

    def enable_all_preparation_steps(self):

        """
        This function ...
        :return:
        """

        self.image_preparer.config.calculate_calibration_uncertainties = True
        self.image_preparer.config.extract_sources = True
        self.image_preparer.config.correct_for_extinction = True
        self.image_preparer.config.convert_unit = True
        self.image_preparer.config.convolve = True
        self.image_preparer.config.rebin = True
        self.image_preparer.config.subtract_sky = True
        self.image_preparer.config.set_uncertainties = True

    # -----------------------------------------------------------------

    def set_preparation_options(self, image, output_path):

        """
        This function ...
        :param image:
        :param output_path:
        :return:
        """

        # Set the attenuation value
        self.image_preparer.config.attenuation = self.attenuations[image.name]

        # If this image is not the reference image, set the appropriate options for rebinning and convolution
        # or this image does not need to be convolved (e.g. SPIRE images)
        if image.name == self.config.reference_image or aniano_names[image.name] is None:

            self.image_preparer.config.rebin = False
            self.image_preparer.config.convolve = False

        # Images that do need to be convolved
        else:

            # Debugging information
            log.debug("Setting the path to the convolution kernel ...")

            # Get the path to the local convolution kernel file
            this_aniano_name = aniano_names[image.name]
            reference_aniano_name = aniano_names[self.config.reference_image]
            kernel_file_path = self.aniano.get_kernel_path(this_aniano_name, reference_aniano_name)

            # Set the kernel path and FWHM
            self.image_preparer.config.convolution.kernel_path = kernel_file_path    # set kernel path
            self.image_preparer.config.convolution.kernel_fwhm = self.reference_fwhm # set kernel FWHM (is a quantity here)

            # Set flags to True
            self.image_preparer.config.rebin = True
            self.image_preparer.config.convolve = True

        # Convolve the SDSS images remotely
        if "SDSS" in image.name: self.image_preparer.config.convolution.remote = "nancy"
        else: self.image_preparer.config.convolution.remote = None

        # Check whether the image has to be sky subtracted
        if image.frames.primary.sky_subtracted:
            log.debug("The " + image.name + " image has already been sky subtracted")
            self.image_preparer.config.subtract_sky = False
        else: self.image_preparer.config.subtract_sky = True # Has yet to be sky subtracted

        # Set the calibration error
        self.image_preparer.config.uncertainties.calibration_error = CalibrationError.from_filter(image.filter)

        # Set the output directory
        self.image_preparer.config.output_path = output_path

        # -----------------------------------------------------------------

        # The units of the Halpha image don't have to be converted
        if "Halpha" in image.name: self.image_preparer.config.convert_unit = False
        else: self.image_preparer.config.convert_unit = True

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
    segments = Image.from_file(segments_path, no_filter=True)

    # Get the different segmentation frames
    galaxy_segments = segments.frames.galaxies
    star_segments = segments.frames.stars if "stars" in segments.frames else None
    other_segments = segments.frames.other_sources if "other_sources" in segments.frames else None

    # Return the regions and segmentation maps
    return galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments

# -----------------------------------------------------------------
