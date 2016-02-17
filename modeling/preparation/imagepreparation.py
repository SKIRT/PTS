#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.imagepreparation Contains the ImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant AstroMagic classes and modules
from ...magic.core import Frame, Source
from ...magic.basics import CoordinateSystem, Mask, SkyRegion
from ...magic.extract import Extractor
from ...magic.subtract import SkySubtractor

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from .unitconversion import UnitConverter
from ...core.tools.logging import log

# -----------------------------------------------------------------

class ImagePreparer(Configurable):

    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ImagePreparer, self).__init__(config, "modeling")

        # -- Attributes --

        # The image
        self.image = None

        # Temporary ...
        self.principal_ellipse_skycoord = None

    # -----------------------------------------------------------------

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup(image)

        # 2. Calculate the calibration uncertainties
        self.calculate_calibration_uncertainties()

        # 3. Extract stars and galaxies from the image
        if self.config.extract_sources: self.extract_sources()

        # 4. If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # 5. If requested, convert the unit
        if self.config.convert_unit: self.convert_unit()

        # 6. If requested, convolve
        if self.config.convolve: self.convolve()

        # 7. If requested, rebin
        if self.config.rebin: self.rebin()

        # 8. If requested, subtract the sky
        if self.config.subtract_sky: self.subtract_sky()

        # 9. If requested, set the uncertainties
        if self.config.set_uncertainties: self.set_uncertainties()

        # 10. If requested, crop to a smaller coordinate grid
        if self.config.crop: self.crop()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # -- Children --

        # Add extractor and sky subtractor
        self.add_child("extractor", Extractor, self.config.extraction)
        self.add_child("sky_subtractor", SkySubtractor, self.config.sky_subtraction)
        self.add_child("unit_converter", UnitConverter, self.config.unit_conversion)

        if self.config.write_steps:

            # -- Source extraction --

            self.extractor.config.write_catalogs = True
            self.extractor.config.write_statistics = True
            self.extractor.config.write_regions = True
            self.extractor.config.write_masked_frames = True

            # -- Sky subtraction --

        # Call the setup function of the base class
        super(ImagePreparer, self).setup()

        # Make a local reference to the passed image (with mask)
        self.image = image

    # -----------------------------------------------------------------

    def calculate_calibration_uncertainties(self):

        """
        This function ...
        :return:
        """

        from . import unitconversion

        jansky_frame = None
        invalid = None
        ab_frame = None

        # GALEX images
        if "GALEX" in self.image.filter.name:

            #FUV: mAB = -2.5 x log10(CPS) + 18.82
            #NUV: mAB = -2.5 x log10(CPS) + 20.08
            invalid = Mask.is_zero_or_less(self.image.frames.primary)
            magnitude_term = {"GALEX FUV": 18.82, "GALEX NUV": 20.08}
            ab_frame = -2.5 * np.log10(self.image.frames.primary) + magnitude_term[self.image.name]

            # Set infinites to zero
            ab_frame[invalid] = 0.0

            # Calculate data in Jansky
            jansky_frame = unitconversion.ab_mag_zero_point.to("Jy").value * np.power(10.0, -2./5. * ab_frame)

            # Add the frame with AB magnitudes and the mask with zeros
            #self.image.add_mask(zeros, "zeros")
            #self.image.add_frame(ab_frame, "abmag")

        # 2MASS images
        elif "2MASS" in self.image.filter.name:

            m_0 = self.image.frames.primary.zero_point
            f_0 = unitconversion.f_0_2mass[self.image.filter.name]
            to_jy_conversion_factor = f_0 * np.power(10.0, -m_0/2.5)

            jansky_frame = self.image.frames.primary * to_jy_conversion_factor

            # Pixels that are zero or less cannot be converted into magnitude scale
            invalid = Mask.is_zero_or_less(jansky_frame)
            ab_frame = -5./2. * np.log10(jansky_frame / unitconversion.ab_mag_zero_point.to("Jy").value)

            # Set infinites to zero
            ab_frame[invalid] = 0.0

            # Add the frame with AB magnitudes and the mask with zeros
            #self.image.add_mask(zeros, "zeros")
            #self.image.add_frame(ab_frame, "abmag")

        # Add the calibration uncertainty
        if "mag" in self.config.uncertainties.calibration_error:

            # Check that the frame has been calculated in Janskys
            assert jansky_frame is not None

            # The calibration uncertainty in AB magnitude
            mag_error = float(self.config.uncertainties.calibration_error.split("mag")[0])

            # a = image[mag] - mag_error
            a = ab_frame - mag_error

            # b = image[mag] + mag_error
            b = ab_frame + mag_error

            # Convert a and b to Jy
            a = unitconversion.ab_mag_zero_point.to("Jy").value * np.power(10.0, -2./5.*a)
            b = unitconversion.ab_mag_zero_point.to("Jy").value * np.power(10.0, -2./5.*b)

            # c = a[Jy] - image[Jy]
            c = a - jansky_frame

            # d = image[Jy] - b[Jy]
            d = jansky_frame - b

            # Calibration errors = max(c, d)
            calibration_errors = np.maximum(c, d)  # element-wise maxima
            calibration_errors[invalid] = 0.0 # set zero where AB magnitude could not be calculated

            ## CONVERT THE CALIBRATION ERRORS TO MJY/SR !! --> NO, TO THE ORIGINAL UNIT OF THE DATA

            relative_calibration_errors = calibration_errors / jansky_frame
            relative_calibration_errors[invalid] = 0.0

            # Check that there are no infinities or nans in the result
            assert not np.any(np.isinf(relative_calibration_errors)) and not np.any(np.isnan(relative_calibration_errors))

            # The actual calibration errors in the same unit as the data
            calibration_frame = self.image.frames.primary * relative_calibration_errors

        elif "%" in self.config.uncertainties.calibration_error:

            # Calculate calibration errors with percentage
            fraction = float(self.config.uncertainties.calibration_error.split("%")[0]) * 0.01
            calibration_frame = self.image.frames.primary * fraction

        else: raise ValueError("Unrecognized calibration error")

        # Add the calibration frame
        self.image.add_frame(calibration_frame, "calibration_errors")

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return:
        """

        # Open the exceptions file if specified
        if self.config.extraction.exceptions_path is not None:

            not_stars = []
            remove_stars = []
            not_saturation = []

            with open(self.config.extraction.exceptions_path, 'r') as exceptions_file:
                for line in exceptions_file:

                    if "not_stars" in line: not_stars = [int(item) for item in line.split(": ")[1].split(" ")]
                    if "remove_stars" in line: remove_stars = [int(item) for item in line.split(": ")[1].split(" ")]
                    if "not_saturation" in line: not_saturation = [int(item) for item in line.split(": ")[1].split(" ")]

            # Set manual indices
            self.extractor.config.stars.manual_indices.not_stars = not_stars
            self.extractor.config.stars.manual_indices.remove_stars = remove_stars
            self.extractor.config.stars.manual_indices.not_saturation = not_saturation

        # Run the extractor
        self.extractor.run(self.image)

        # Print the FWHM if it has been fitted by the extractor
        #print("image fwhm=", self.image.fwhm)
        #print("use frame fwhm=", self.extractor.star_extractor.config.use_frame_fwhm)
        fwhm_not_fitted = (self.extractor.star_extractor.config.use_frame_fwhm and self.image.fwhm is not None)
        if not fwhm_not_fitted:

            # Get the FWHM from the star extractor (in pixels -> in arcsec)
            fwhm = self.extractor.star_extractor.fwhm
            fwhm = fwhm * self.image.frames.primary.xy_average_pixelscale.to("arcsec/pix").value

            # Debug info
            log.debug("The fitted FWHM is " + str(fwhm) + " arcseconds")

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("extracted.fits")

        # Get ellipse of principal galaxy
        self.principal_ellipse_skycoord = self.extractor.galaxy_extractor.principal_sky_ellipse

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting image for galactic extinction ...")

        # Correct the primary frame for galactic extinction
        self.image.frames[self.config.primary] *= 10**(0.4 * self.config.attenuation)

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("corrected_for_extinction.fits")

    # -----------------------------------------------------------------

    def convert_unit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting image to surface brightness units ...")

        # Run the unit conversion
        self.unit_converter.run(self.image)

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("converted_unit.fits")

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Convolving the image with kernel " + self.config.convolution.kernel_path + " ...")

        # Open the kernel frame
        kernel = Frame.from_file(self.config.convolution.kernel_path)

        # Set the kernel FWHM
        kernel.fwhm = self.config.convolution.kernel_fwhm

        # Convolve the image (the primary and errors frame)
        self.image.convolve(kernel)

        # Save convolved frame
        if self.config.write_steps: self.write_intermediate_result("convolved.fits")

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning the image to the pixel grid of " + self.config.rebinning.rebin_to + " ...")

        # Get the coordinate system of the reference frame
        reference_system = CoordinateSystem.from_file(self.config.rebinning.rebin_to)

        # Add a mask to the image that covers the complete current pixel grid
        self.image.add_mask(Mask.full_like(self.image.frames.primary), "padded")

        # Rebin the image (the primary and errors frame)
        self.image.rebin(reference_system)

        # Invert the 'padded' mask -> this mask now covers pixels added to the frame after rebinning
        self.image.masks.padded = self.image.masks.padded.inverse().disk_dilation(radius=10)

        # Save rebinned frame
        if self.config.write_steps: self.write_intermediate_result("rebinned.fits")

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # Write out the principal_sky_ellipse and saturation_region in sky coordinates
        principal_sky_region = SkyRegion()
        principal_sky_region.append(self.principal_ellipse_skycoord)
        principal_sky_region_path = self.full_output_path("principal_skycoord.reg")
        principal_sky_region.save(principal_sky_region_path)

        saturation_sky_region_path = self.full_output_path("saturation_skycoord.reg")
        self.extractor.star_extractor.saturation_region.save(saturation_sky_region_path)

        # Run the sky extraction
        self.sky_subtractor.run(self.image, self.principal_ellipse_skycoord, self.extractor.star_extractor.saturation_region)

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("sky_subtracted.fits")

    # -----------------------------------------------------------------

    def set_uncertainties(self):

        """
        This function ...
        :return:
        """

        sky_mask = self.image.masks.sky

        photutils_estimated_sky = np.ma.masked_array(np.asarray(self.image.frames.phot_sky), mask=sky_mask)
        photutils_rms = np.ma.masked_array(np.asarray(self.image.frames.phot_rms), mask=sky_mask)

        #plotting.plot_box(photutils_rms, title="masked rms")

        # LARGE SCALE VARIATIONS = STANDARD DEVIATION OF ESTIMATED SKY PIXELS
        large_scale_variations_error = photutils_estimated_sky.std()

        # PIXEL TO PIXEL NOISE = MEAN OF RMS
        pixel_to_pixel_noise = np.ma.mean(photutils_rms) # Mean pixel-by-pixel variations

        # Add all the errors quadratically: existing error map + large scale variations + pixel to pixel noise + calibration error
        self.image.frames.errors = np.sqrt(self.image.frames.errors**2 + large_scale_variations_error**2 + \
                                           + pixel_to_pixel_noise**2 + self.image.frames.calibration_errors**2)

    # -----------------------------------------------------------------

    def crop(self):

        """
        This function ...
        :return:
        """

        # Get the cropping limits
        #x_min = self.config.cropping.limits[0]
        #x_max = self.config.cropping.limits[1]
        #y_min = self.config.cropping.limits[2]
        #y_max = self.config.cropping.limits[3]

        # Crop the image (the primary and errors frame)
        #self.image.crop(x_min, x_max, y_min, y_max)

        pass

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # If requested, write out the result
        if self.config.write_result: self.write_result()

    # -----------------------------------------------------------------

    def write_result(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(self.config.writing.result_path)

        # Inform the user
        log.info("Writing resulting image to " + path + " ...")

        # Write out the resulting image
        self.image.save(path, origin=self.name)

    # -----------------------------------------------------------------

    def write_intermediate_result(self, path):

        """
        This function ...
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(path)

        # Inform the user
        log.info("Writing intermediate result to " + path + " ...")

        # Write out the image
        self.image.save(path, origin=self.name)

# -----------------------------------------------------------------
