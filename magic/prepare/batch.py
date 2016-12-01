#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.batch Contains the BatchImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from multiprocessing import Pool
import numpy as np

# Import astronomical modules
from astropy.units import Unit
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.basics.configurable import Configurable
from ..misc.calibration import CalibrationError
from ..misc.extinction import GalacticExtinction
from ..core.frame import Frame, sum_frames_quadratically
from ..core.image import Image
from ..core.dataset import DataSet
from ..region.list import SkyRegionList
from ...magic.misc.kernels import AnianoKernels
from ...magic.sources.extractor import SourceExtractor
from ...magic.sky.skysubtractor import SkySubtractor
from ...core.basics.animation import Animation
from ...magic.animation.sourceextraction import SourceExtractionAnimation
from ...magic.animation.imageblink import ImageBlinkAnimation
from ...core.tools import time
from ...magic.core.remote import RemoteImage
from ...magic.core.kernel import ConvolutionKernel
from ...modeling.preparation import unitconversion
from ..basics.mask import Mask
from ...core.basics.composite import SimplePropertyComposite

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

        self.convolution_filter = kwargs.pop("convolution_filter")
        self.rebinning_filter = kwargs.pop("rebinning_filter")
        self.not_rebinned = kwargs.pop("not_rebinned")
        self.not_convolved = kwargs.pop("not_convolved")

# -----------------------------------------------------------------

class BatchImagePreparer(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(BatchImagePreparer, self).__init__(config)

        # The images
        self.images = dict()

        # The attenuations
        self.attenuations = dict()

        # The output paths
        self.output_paths = dict()

        # The process pool
        self.pool = None

        # Rebinning wcs and convolution filter
        self.rebinning_filter = None
        self.rebinning_wcs = None
        self.convolution_filter = None
        self.convolution_fwhm = None

        # The Aniano kernel service
        self.aniano = AnianoKernels()

        # Don't convolve or rebin
        self.dont_convolve = []
        self.dont_rebin = []

        # Set the principal ellipse and saturation region in sky coordinates
        #self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(self.image.wcs)
        #self.image_preparer.saturation_region_sky = saturation_region.to_sky(self.image.wcs) if saturation_region is not None else None
        self.principal_sky_regions = dict()
        self.saturation_sky_regions = dict()

        # The noise maps
        self.noise_maps = dict()

        # The output dataset
        self.dataset = DataSet()

    # -----------------------------------------------------------------

    def add_image(self, name, image, output_path=None):

        """
        This function ...
        :param name:
        :param image:
        :param output_path:
        :return:
        """

        # Check if name not already used
        if name in self.images: raise ValueError("Already an image with the name " + name)

        # Set the image
        self.images[name] = image

        # Check the output path
        if self.config.write and output_path is None: raise ValueError("Output path must be specified if 'write' option is enabled in configuration")

        # Set the output path
        if output_path is not None: self.output_paths[name] = output_path

    # -----------------------------------------------------------------

    @property
    def min_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the images
        for name in self.images:
            if fwhm is None or self.images[name].fwhm < fwhm: fwhm = self.images[name].fwhm

        # Return the minimum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the images
        for name in self.images:
            if fwhm is None or self.images[name].fwhm > fwhm: fwhm = self.images[name].fwhm

        # Return the maximum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.images:

            wcs = self.images[name].wcs
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.images:

            wcs = self.images[name].wcs
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegionList()

        # Add the bounding boxes as sky rectangles
        for name in self.images: boxes_region.append(self.images[name].wcs.bounding_box)

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

    # -----------------------------------------------------------------

    @property
    def center_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.bounding_box.center

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the galactic attenuation
        self.get_extinction()

        # 3. Prepare the images
        self.prepare()

        # 4. Writing
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchImagePreparer, self).setup()

        # Initialize the process pool
        self.pool = Pool(processes=self.config.nprocesses)

        # Load the images (from config or input kwargs)
        if "images" in kwargs: self.images = kwargs.pop("images")
        elif "dataset" in kwargs:
            dataset = kwargs.pop("dataset")
            self.images = dataset.get_images()
            self.output_paths = dataset.get_directory_paths()
        else: self.load_images()

        # Set rebinning wcs
        self.set_rebinning_wcs()

        # Set convolution filter
        self.set_convolution_filter()

    # -----------------------------------------------------------------

    def set_rebinning_wcs(self):

        """
        This function ...
        :return:
        """

        max_pixelscale_wcs = None
        max_pixelscale_filter = None

        # Loop over the images
        for name in self.images:

            # Get the pixelscale
            pixelscale = self.images[name].average_pixelscale

            # Ignore this frame if its pixelscale is larger than the minimum pixelscale that is specified by the user
            if self.config.max_pixelscale is not None and pixelscale > self.config.max_pixelscale:
                self.dont_rebin.append(name)
                continue

            # Check if the pixelscale is greater than that of the previous frame (but still below the maximum)
            if max_pixelscale_wcs is None or pixelscale > max_pixelscale_wcs.average_pixelscale:
                max_pixelscale_wcs = self.images[name].wcs
                max_pixelscale_filter = self.images[name].filter

        # Now set the rebinning filter and wcs
        self.rebinning_filter = max_pixelscale_filter
        self.rebinning_wcs = max_pixelscale_wcs

    # -----------------------------------------------------------------

    def set_convolution_filter(self):

        """
        This function ...
        :return:
        """

        max_fwhm = None
        max_fwhm_filter = None

        # Loop over the images
        for name in self.images:

            # Get the FWHM
            fwhm = self.images[name].fwhm

            # Ignore this frame if its FWHM is higher than the maximum FWHM specified by the user
            if self.config.max_fwhm is not None and fwhm > self.config.max_fwhm:
                self.dont_convolve.append(name)
                continue

            # Check if the FWHM is greater than that of the previous frame (but still below the maximum)
            if max_fwhm is None or fwhm > max_fwhm:
                max_fwhm = fwhm
                max_fwhm_filter = self.images[name].filter

        # Set the convolution filter and FWHM
        self.convolution_filter = max_fwhm_filter
        self.convolution_fwhm = max_fwhm

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Create new dataset
        if self.config.dataset.endswith(".fits"):

            # Load the image
            image = Image.from_file(self.config.dataset)

            # Determine the name for this image
            name = str(image.filter)

            # Determine the directory path
            directory_path = fs.directory_of(self.config.dataset)

            # Add the frame
            self.add_image(name, image, output_path=directory_path)

        # Load dataset from file
        elif self.config.dataset.endswith(".dat"):

            # Get the dataset
            dataset = DataSet.from_file(self.config.dataset)

            # Get the images
            self.images = dataset.get_images()

            # Get the directory paths
            self.output_paths = dataset.get_directory_paths()

        # Invalid value for 'dataset'
        else: raise ValueError("Parameter 'dataset' must be filename of a dataset file (.dat) or a FITS file (.fits)")

    # -----------------------------------------------------------------

    def get_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galactic extinction ...")

        # Create the galactic extinction calculator
        extinction = GalacticExtinction(self.center_coordinate)

        # Loop over the images
        for label in self.images:

            # Get the filter
            fltr = self.images[label].filter

            # Get the extinction
            self.attenuations[label] = extinction.extinction_for_filter(fltr)

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the images ...")

        # 1. Extract stars and galaxies from the image
        if self.config.extract_sources: self.extract_sources()

        # 2. If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # 3. If requested, convert the unit
        if self.config.convert_units: self.convert_units()

        # 4. If requested, convolve
        if self.config.convolve: self.convolve()

        # 5. Rebin
        if self.config.rebin: self.rebin()

        # 6. Subtract the sky
        if self.config.subtract_sky: self.subtract_sky()

        # 7. Calculate the calibration uncertainties
        if self.config.calculate_calibration_uncertainties: self.calculate_calibration_uncertainties()

        # 8. Calculate the error maps
        if self.config.create_errormaps: self.create_errormaps()

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting sources ...")

        # Loop over the images that have to be source-extracted
        for label, image in self.to_be_extracted():

            # Check if sources directory exists
            sources_path = fs.join(self.output_paths[label], "sources")
            if not fs.is_directory(sources_path): raise IOError("Sources directory not present for '" + label + "' image")

            # Execute
            # result = self.pool.apply_async(_extract_sources, args=(,))  # All simple types (strings) ?

            # Get the configuration for the source extractor
            config = self.config.source_extraction

            # Extract the sources
            principal_sky_region, saturation_sky_region = _extract_sources(image, config, sources_path, visualisation_path=self.config.visualisation_path)

            # Add to dict
            self.principal_sky_regions[label] = principal_sky_region
            self.saturation_sky_regions[label] = saturation_sky_region

            # Save intermediate result if requested
            if self.config.write_steps: image.save(fs.join(self.output_paths[label], "extracted.fits"))

        # Inform the user
        log.success("Sources extracted")

        # Multiprocessing:

            # Add the result
            # results.append(result)

        #self.dust_masses = [result.get() for result in results]

        # Close and join the process pool
        #self.pool.close()
        #self.pool.join()

    # -----------------------------------------------------------------

    def to_be_extracted(self):

        """
        This function ...
        :return:
        """

        # Loop over the image
        for label in self.images:

            # Check if source extraction is required
            if self._needs_extinction_correction(label): yield label, self.images[label]

    # -----------------------------------------------------------------

    def _needs_source_extraction(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return not self.images[label].source_extracted

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting for galactic extinction ...")

        # Loop over the images that have to be corrected for extinction
        for label, image in self.to_be_corrected():

            # Do the correction
            # Correct all data frames for galactic extinction (primary, poisson error frame, ...)
            image *= 10 ** (0.4 * self.attenuations[label])

            # IMPORTANT: SET FLAG
            image.corrected_for_extinction = True

            # Save intermediate result if requested
            if self.config.write_steps: image.save(fs.join(self.output_paths[label], "corrected_for_extinction.fits"))

        # Inform the user
        log.success("Corrected for extinction")

    # -----------------------------------------------------------------

    def to_be_corrected(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for label in self.images:

            # Check whether unit conversion is necessary
            if self._needs_unit_conversion(label): yield label, self.images[label]

    # -----------------------------------------------------------------

    def _needs_extinction_correction(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return not self.images[label].extinction_corrected

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting the units ...")

        # Loop over the images that have to be converted in units
        for label, image in self.to_be_converted():

            # Debugging
            log.debug("Converting the units of the '" + label + "' image ...")

            # Check the units of the image
            assert image.unit == Unit("Jy/pix")

            # Get pixelscale
            pixelscale = image.average_pixelscale

            # Conversion from Jy / pix to MJy / pix
            conversion_factor = 1e-6

            # Conversion from MJy / pix to MJy / sr
            conversion_factor *= (1.0 / pixelscale ** 2).to("pix2/sr").value

            # Multiply the image with the conversion factor
            image *= conversion_factor

            # We can only do unit conversion to MJy/sr at the moment
            assert self.config.unit_conversion.to_unit == Unit("MJy/sr")

            # Set the new unit
            image.unit = self.config.unit_conversion.to_unit

            # Save intermediate result if requested
            if self.config.write_steps: image.save(fs.join(self.output_paths[label], "converted_unit.fits"))

        # Inform the user
        log.success("Units converted")

    # -----------------------------------------------------------------

    def to_be_converted(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for label in self.images:

            # Check whether unit conversion is necessary
            if self._needs_unit_conversion(label): yield label, self.images[label]

    # -----------------------------------------------------------------

    def _needs_unit_conversion(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.images[label].unit != self.config.unit_conversion.to_unit

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Convolving the images to a common resolution ...")

        # Loop over the images that need to be convolved
        for label, image in self.to_be_convolved():

            # Debugging
            log.debug("Convolving the '" + label + "' image ...")

            # Filters
            from_filter = image.filter
            to_filter = self.convolution_filter

            # FWHMs
            from_fwhm = image.fwhm
            to_fwhm = self.convolution_fwhm

            # Get the kernel path
            kernel_file_path = self.aniano.get_kernel_path(from_filter, to_filter, from_fwhm=from_fwhm, to_fwhm=to_fwhm)

            # Set the kernel path and FWHM
            #self.image_preparer.config.convolution.kernel_path = kernel_file_path  # set kernel path
            #self.image_preparer.config.convolution.kernel_fwhm = self.reference_fwhm  # set kernel FWHM (is a quantity here)

            # Do convolution
            kernel = _convolve(image, kernel_file_path, self.convolution_fwhm, visualisation_path=self.config.visualisation_path)

            # Save intermediate result if requested
            if self.config.write_steps:

                # Save the kernel frame for manual inspection
                kernel_path = fs.join(self.output_paths[label], "kernel.fit")
                kernel.save(kernel_path)

                # Save the image
                image.save(fs.join(self.output_paths[label], "convolved.fits"))

        # Inform the user
        log.success("Convolution finished")

    # -----------------------------------------------------------------

    def to_be_convolved(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for label in self.images:

            # Check if convolution is necessary
            if self._needs_convolution(label): yield label, self.images[label]

    # -----------------------------------------------------------------

    def _needs_convolution(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # If this frame does not need to be convolved
        if label in self.dont_convolve: return False

        # Check the resolution
        if self.images[label].fwhm == self.convolution_fwhm: return False

        # Else, return True
        return True

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning to a common pixel grid ...")

        # Loop over the images that have to rebinned
        for label, image in self.to_be_rebinned():

            # Debugging
            log.debug("Rebinning the '" + label + "' image ...")

            # Do rebinning for this image
            _rebin(image, reference_wcs=self.rebinning_wcs, exact=self.config.rebinning.exact)

            # Save intermediate result if requested
            if self.config.write_steps: image.save(fs.join(self.output_paths[label], "rebinned.fits"))

        # Inform the user
        log.success("Rebinning finished")

    # -----------------------------------------------------------------

    def to_be_rebinned(self):

        """
        This function ...
        :return:
        """

        # Loop over the images
        for label in self.images:

            # Check if rebinning is necessary
            if not self._needs_rebinning(label): yield label, self.images[label]

    # -----------------------------------------------------------------

    def _needs_rebinning(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # If this frame does not need to be rebinned
        if label in self.dont_rebin: return False

        # Check the coordinate system
        if self.images[label].wcs == self.rebinning_wcs: return False

        # Else, return True
        return True

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the sky subtraction ...")

        # Loop over the images that have to be sky-subtracted
        for label, image in self.to_be_subtracted():

            # Debugging
            log.debug("Subtracting the sky from the '" + label + "' image ...")

            # Get regions
            principal_sky_region = self.principal_sky_regions[label]
            saturation_sky_region = self.saturation_sky_regions[label]

            # Get the configuration
            config = self.config.sky_subtraction

            # Do the sky subtraction
            noise_frame = _subtract_sky(image, config, principal_sky_region, saturation_sky_region, visualisation_path=self.config.visualisation_path)

            # Add the noise frame
            self.noise_maps[label] = noise_frame

            # Save intermediate result if requested
            if self.config.write_steps: image.save(fs.join(self.output_paths[label], "sky_subtracted.fits"))

        # Inform the user
        log.success("Sky subtraction finished")

    # -----------------------------------------------------------------

    def to_be_subtracted(self):

        """
        This function ...
        :return:
        """

        # Loop over the image
        for label in self.images:

            # Check whether the image has to be sky subtracted
            if self._needs_sky_subtraction(label): yield label, self.images[label]

    # -----------------------------------------------------------------

    def _needs_sky_subtraction(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        #if image.frames.primary.sky_subtracted:
        #    log.debug("The " + image.name + " image has already been sky subtracted")
        #    self.image_preparer.config.subtract_sky = False
        #else: self.image_preparer.config.subtract_sky = True  # Has yet to be sky subtracted

        return not self.images[label].sky_subtracted

    # -----------------------------------------------------------------

    def calculate_calibration_uncertainties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating calibration uncertainty maps ...")

        # Loop over the images
        for label in self.images:

            # Debugging
            log.debug("Calculating calibration uncertainties for '" + label + "' image ...")

            # Get the image
            image = self.images[label]

            # Get the calibration error
            calibration_error = CalibrationError.from_filter(image.filter)

            # Add the calibration uncertainty defined in (AB) magnitude
            if calibration_error.magnitude:

                # -----------------------------------------------------------------

                # Convert the frame into AB magnitudes
                invalid = Mask.is_zero_or_less(image.frames.primary)
                ab_frame = unitconversion.jansky_to_ab(image.frames.primary)
                # Set infinites to zero
                ab_frame[invalid] = 0.0

                # -----------------------------------------------------------------

                # The calibration uncertainty in AB magnitude
                mag_error = calibration_error.value

                # a = image[mag] - mag_error
                a = ab_frame - mag_error

                # b = image[mag] + mag_error
                b = ab_frame + mag_error

                # Convert a and b to Jy
                a = unitconversion.ab_mag_zero_point.to("Jy").value * np.power(10.0, -2. / 5. * a)
                b = unitconversion.ab_mag_zero_point.to("Jy").value * np.power(10.0, -2. / 5. * b)

                # c = a[Jy] - image[Jy]
                # c = a - jansky_frame
                c = a - image.frames.primary

                # d = image[Jy] - b[Jy]
                # d = jansky_frame - b
                d = image.frames.primary - b

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

                new_invalid = Mask.is_nan(calibration_frame) + Mask.is_inf(calibration_frame)
                calibration_frame[new_invalid] = 0.0

                # Check that there are no infinities or nans in the result
                assert not np.any(np.isinf(calibration_frame)) and not np.any(np.isnan(calibration_frame))

                # Make frame from numpy array
                calibration_frame = Frame(calibration_frame)

                # -----------------------------------------------------------------

            # The calibration uncertainty is expressed in a percentage (from the flux values)
            elif calibration_error.percentage:

                # Calculate calibration errors with percentage
                fraction = calibration_error.value * 0.01
                calibration_frame = image.frames.primary * fraction

            # Unrecognized calibration error (not a magnitude, not a percentage)
            else: raise ValueError("Unrecognized calibration error")

            # Add the calibration frame
            image.add_frame(calibration_frame, "calibration_errors")

        # Inform the user
        log.success("Calibration uncertainties calculated")

    # -----------------------------------------------------------------

    def create_errormaps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the error maps ...")

        # Loop over the frames
        for label, image in self.to_have_errormaps_created():

            # Inform the user
            log.debug("Calculating the error map for the '" + label + "' image ...")

            # Create a list to contain (the squares of) all the individual error contributions so that we can sum these arrays element-wise later
            error_maps = []

            # Add the Poisson errors
            if "poisson_errors" in image.frames: error_maps.append(image.frames.poisson_errors)

            # Get the noise frame
            noise_frame = self.noise_maps[label]

            # Add the sky errors
            error_maps.append(noise_frame)

            # Add the calibration errors
            error_maps.append(image.frames.calibration_errors)

            # Add additional error frames indicated by the user
            if self.config.error_frame_names is not None:
                for error_frame_name in self.config.error_frame_names: error_maps.append(image.frames[error_frame_name])

            # Calculate the final error map
            errors = sum_frames_quadratically(*error_maps)

            # Add the combined errors frame
            image.add_frame(errors, "errors")

        # Inform the user
        log.success("Error maps created")

    # -----------------------------------------------------------------

    def to_have_errormaps_created(self):

        """
        This function ...
        :return:
        """

        return self.images.items()

    # -----------------------------------------------------------------

    @lazyproperty
    def statistics(self):

        """
        This function ...
        :return:
        """

        # The statistics
        statistics = PreparationStatistics(convolution_filter=self.convolution_filter, rebinning_filter=self.rebinning_filter, not_convolved=self.dont_convolve, not_rebinned=self.dont_rebin)
        return statistics

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the images
        self.write_images()

        # Write the new dataset
        if self.config.writing.dataset_path is not None: self.write_dataset()

        # Write statistics
        if self.config.writing.statistics_path is not None: self.write_statistics()

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Loop over the images
        for name in self.images:

            # Get the output path
            if name not in self.output_paths: raise ValueError("Output directory not specified for '" + name + "' image")

            # Determine result path
            path = fs.join(self.output_paths[name], "result.fits")

            # Debugging
            log.debug("Writing the '" + name + "' image to '" + path + "' ...")

            # Add an entry to the output dataset
            self.dataset.add_path(name, path)

            # Save the image
            self.images[name].save(path)

    # -----------------------------------------------------------------

    def write_dataset(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dataset file to '" + self.config.writing.dataset_path + "' ...")

        # Write the dataset
        self.dataset.saveto(self.config.writing.dataset_path)

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the statistics to '" + self.config.writing.statistics_path + "' ...")

        # Write the statistics
        self.statistics.save(self.config.writing.statistics_path)

# -----------------------------------------------------------------

def _extract_sources(image, config, sources_path, visualisation_path=None):

    """
    This function ...
    :return:
    """

    # Inform the user
    log.info("Extracting the sources ...")

    # Create an animation to show the result of this step
    if visualisation_path is not None:

        # Create an animation
        animation = ImageBlinkAnimation()
        animation.add_image(image.frames.primary)

        # Create an animation to show the progress of the SourceExtractor
        source_extractor_animation = SourceExtractionAnimation(image.frames.primary)

    else:

        animation = None
        source_extractor_animation = None

    # Create the source extractor
    extractor = SourceExtractor(config)

    # Run the extractor OLD WAY
    #extractor.run(self.image.frames.primary, self.galaxy_region, self.star_region, self.saturation_region,
    #               self.other_region, self.galaxy_segments, self.star_segments, self.other_segments,
    #               source_extractor_animation)

    # Set configuration settings
    extractor.config.input = sources_path

    # Run the extraction
    special_region = None
    extractor.run(frame=image.frames.primary, animation=source_extractor_animation, special_region=special_region)

    # Get the saturation region
    saturation_region = extractor.saturation_region

    # Write the animation
    if visualisation_path is not None:

        # Determine the path to the animation
        path = fs.join(visualisation_path, time.unique_name(image.name + "_sourceextraction") + ".gif")

        # Save the animation
        source_extractor_animation.save(path)

        # ...
        animation.add_image(image.frames.primary)

        # Determine the path to the animation
        path = fs.join(visualisation_path, time.unique_name(image.name + "_imagepreparation_extractsources") + ".gif")

        # Save the animation
        animation.save(path)

    # Add the sources mask to the image
    image.add_mask(extractor.mask, "sources")

    # Get the principal shape in sky coordinates
    principal_shape_sky = extractor.principal_shape.to_sky(image.wcs)

    # Get the saturation region in sky coordinates
    saturation_region_sky = saturation_region.to_sky(image.wcs) if saturation_region is not None else None

    # IMPORTANT: SET FLAG
    image.source_extracted = True

    # Return
    return principal_shape_sky, saturation_region_sky

# -----------------------------------------------------------------

def _convolve(image, kernel_path, kernel_fwhm, visualisation_path=None, host_id=None):

    """
    This function ...
    :return:
    """

    # Inform the user
    log.info("Convolving the image with kernel " + kernel_path + " ...")

    # Create an animation to show the result of this step
    if visualisation_path is not None:

        # Create an animation
        animation = ImageBlinkAnimation()
        animation.add_image(image.frames.primary)

    else: animation = None

    # Open the convolution kernel
    kernel = ConvolutionKernel.from_file(kernel_path, fwhm=kernel_fwhm)

    # Prepare the kernel
    kernel.prepare_for(image)

    # Check whether the convolution has to be performed remotely
    if host_id is not None:

        # Inform the user
        log.info("Convolution will be performed remotely on host '" + host_id + "' ...")

        # Create remote image, convolve and make local again
        remote_image = RemoteImage.from_local(image, host_id)
        remote_image.convolve(kernel, allow_huge=True)
        new_image = remote_image.to_local()

        # Load the properties and content of the new image into the 'working' image
        image.load_image(new_image)

    # The convolution is performed locally. # Convolve the image (the primary and errors frame)
    else: image.convolve(kernel, allow_huge=True)

    # Write the animation
    if visualisation_path is not None:

        animation.add_image(image.frames.primary)

        # Determine the path to the animation
        path = fs.join(visualisation_path, time.unique_name(image.name + "_imagepreparation_convolve") + ".gif")

        # Save the animation
        animation.save(path)

    # Return the kernel
    return kernel

# -----------------------------------------------------------------

def _rebin(image, reference_wcs, exact, host_id=None):

    """
    This function ...
    :param image:
    :param reference_wcs:
    :return:
    """

    # Inform the user
    log.info("Rebinning the image to a pixelsize of " + str(reference_wcs.average_pixelscale) + "...")

    # Check whether the rebinning has to be performed remotely
    if host_id is not None:

        # Inform the user
        log.info("Rebinning will be performed remotely on host '" + host_id + "' ...")

        # Create remote image, rebin and make local again
        remote_image = RemoteImage.from_local(image, host_id)
        remote_image.rebin(reference_wcs, exact=exact)
        new_image = remote_image.to_local()

        # Load the properties and content of the new image into the 'working' image
        image.load_image(new_image)

    # Rebin the image locally (the primary and errors frame)
    else: image.rebin(reference_wcs)

# -----------------------------------------------------------------

def _subtract_sky(image, config, principal_sky_region, saturation_sky_region=None, visualisation_path=None, sky_apertures_path=None):

    """
    This function ...
    :param image:
    :param principal_sky_region:
    :param saturation_sky_region:
    :return:
    """

    # Inform the user
    log.info("Subtracting the sky ...")

    # Create an animation to show the result of this step
    if visualisation_path is not None:

        # Create an animation to show the result of the sky subtraction step
        animation = ImageBlinkAnimation()
        animation.add_image(image.frames.primary)

        # Create an animation to show the progress of the SkySubtractor
        skysubtractor_animation = Animation()

    else:

        animation = None
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
    sky_subtractor.run(image.frames.primary, principal_shape, image.masks.sources, extra_mask, saturation_region, skysubtractor_animation)

    # Add the sky frame to the image
    image.add_frame(sky_subtractor.sky_frame, "sky")

    # Add the mask that is used for the sky estimation
    image.add_mask(sky_subtractor.mask, "sky")

    # Add the sky noise frame
    image.add_frame(sky_subtractor.noise_frame, "sky_errors")

    # Write sky annuli maps if requested
    if sky_apertures_path is not None:

        # Write the sky region
        region_path = fs.join(sky_apertures_path, "annulus.reg")
        sky_subtractor.region.save(region_path)

        # Write the apertures frame
        apertures_frame_path = fs.join(sky_apertures_path, "apertures.fits")
        sky_subtractor.apertures_frame.save(apertures_frame_path)

        # Write the apertures mean frame
        apertures_mean_path = fs.join(sky_apertures_path, "apertures_mean.fits")
        sky_subtractor.apertures_mean_frame.save(apertures_mean_path)

        # Write the apertures noise frame
        apertures_noise_path = fs.join(sky_apertures_path, "apertures_noise.fits")
        sky_subtractor.apertures_noise_frame.save(apertures_noise_path)

    # Write the animation
    if visualisation_path is not None:

        # Determine the path to the animation
        path = fs.join(visualisation_path, time.unique_name(image.name + "_skysubtraction") + ".gif")

        # Save the animation
        skysubtractor_animation.save(path)

        # ...
        animation.add_image(image.frames.primary)

        # Determine the path to the animation
        path = fs.join(visualisation_path, time.unique_name(image.name + "_imagepreparation_subtractsky") + ".gif")

        # Save the animation
        animation.save(path)

    # IMPORTANT: SET FLAG
    image.sky_subtracted = True

    # Return the noise frame
    return sky_subtractor.noise_frame

# -----------------------------------------------------------------
