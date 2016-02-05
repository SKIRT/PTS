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
from ...magic.basics import Region, CoordinateSystem
from ...magic.extract import Extractor
from ...magic.subtract import SkySubtractor
from ...magic.tools import regions, cropping

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

    # -----------------------------------------------------------------

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup(image)

        # 3. Extract stars and galaxies from the image
        if self.config.extract_sources: self.extract_sources()

        # 4. If requested, subtract the sky
        #if self.config.subtract_sky: self.subtract_sky()

        # 5. If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # 6. If requested, convert the unit
        if self.config.convert_unit: self.convert_unit()

        # 7. If requested, convolve
        if self.config.convolve: self.convolve()

        # 8. If requested, rebin
        if self.config.rebin: self.rebin()

        # 4. If requested, subtract the sky
        if self.config.subtract_sky: self.subtract_sky()

        # 9. If requested, set the uncertainties
        #if self.config.set_uncertainties: self.set_uncertainties()

        # 10. If requested, crop
        #if self.config.crop: self.crop()

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

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("extracted.fits")

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

        # Rebin the image (the primary and errors frame)
        self.image.rebin(reference_system)

        # Save rebinned frame
        if self.config.write_steps: self.write_intermediate_result("rebinned.fits")

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # Run the sky extraction
        self.sky_subtractor.run(self.image, self.extractor.galaxy_extractor.principal_sky_ellipse, self.extractor.star_extractor.saturation_region)

        # Print the statistics of the sky frame
        log.info("Mean sky level = " + str(self.sky_subtractor.mean))
        log.info("Median sky level = " + str(self.sky_subtractor.median))
        log.info("Standard deviation of sky = " + str(self.sky_subtractor.stddev))

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("sky_subtracted.fits")

    # -----------------------------------------------------------------

    def set_uncertainties(self):

        """
        This function ...
        :return:
        """

        # Create noise region
        region = Region.from_file(self.config.uncertainties.noise_path, self.image.wcs)

        means = []
        stddevs = []

        for shape in region:

            # Create ellipse
            ellipse = regions.ellipse(shape)

            # Create source for ellipse
            source = Source.from_ellipse(self.image.frames.primary, ellipse, 1.5)

            # Calculate the mean and standard deviation in the box
            mean = np.ma.mean(np.ma.masked_array(source.cutout, mask=source.background_mask))
            stddev = np.median(np.ma.masked_array(source.cutout, mask=source.background_mask).compressed())

            # Add the mean and standard deviation to the appropriate list
            means.append(mean)
            stddevs.append(stddev)

        # Create numpy arrays
        means = np.array(means)
        stddevs = np.array(stddevs)

        # Calculate the global uncertainty
        uncertainty = np.sqrt(np.std(means)**2 + np.median(stddevs)**2)



        # Add the calibration uncertainty
        if "mag" in self.config.uncertainties.calibration_error:

            mag_error = float(self.config.uncertainties.calibration_error.split("mag")[0])

            # a = image[mag] - mag_error
            a = image_mag - mag_error

            # b = image[mag] + mag_error
            b = image_mag + mag_error

            # Convert a and b to Jy
            a = a
            b = b

            # c = a[Jy] - image[Jy]
            c = None

            # d = image[Jy] - b[Jy]
            d = None

            # Calibration errors = max(c, d)
            calibration_errors = Frame(np.zeros(self.image.shape))
            for x in self.image.xsize:
                for y in self.image.ysize:
                    calibration_errors = max(c[y,x], d[y,x])

        elif "%" in self.config.uncertainties.calibration_error:

            fraction = float(self.config.uncertainties.calibration_error.split("%")[0]) * 0.01

            calibration_errors = self.image.frames.primary * fraction

        else: raise ValueError("Unrecognized calibration error")

        # Add the uncertainty AND THE CALIBRATION ERRORS to the errors quadratically
        # IS THIS THE RIGHT WAY ??
        self.image.frames.errors = np.sqrt(np.power(self.image.frames.errors, 2) + uncertainty**2 + np.power(calibration_errors, 2))

    # -----------------------------------------------------------------

    def set_uncertainties_old(self):

        """
        This function ...
        :return:
        """

        # Create noise region
        region = Region.from_file(self.config.uncertainties.noise_path, self.image.frames[self.config.primary].wcs)

        # Initialize lists for the mean value and standard deviation in the different shapes
        means = []
        stddevs = []

        # Loop over all shapes
        for shape in region:


            # TODO: use angle

            x_center, y_center, x_radius, y_radius, angle = regions.ellipse_parameters(shape)

            box, x_min, x_max, y_min, y_max = cropping.crop(self.image.frames[self.config.primary], x_center, y_center, x_radius, y_radius)

            # Calculate the mean and standard deviation in the box
            mean = np.mean(box)
            stddev = np.std(box)

            means.append(mean)
            stddevs.append(stddev)

        means = np.array(means)
        stddevs = np.array(stddevs)

        # Calculate the global uncertainty
        uncertainty = np.sqrt(np.std(means)**2 + np.median(stddevs)**2)

        # If there is no errors frame
        if not self.config.errors in self.image.frames:

            wcs = self.image.frames[self.config.primary].wcs
            #pixelscale = self.image.frames[self.config.primary].pixelscale
            unit = self.image.frames[self.config.primary].unit
            filter = self.image.frames[self.config.primary].filter
            selected = True

            # Create the errors frame (select it)
            frame = Frame(np.full(self.image.frames[self.config.primary].shape, uncertainty), wcs, None, selected, unit, self.config.errors, filter)

            # Add the errors frame to the image
            self.image.add_frame(frame, self.config.errors)

        # If there is an errors frame, add the uncertainty to the errors quadratically
        else: self.image.frames[self.config.errors] = np.sqrt(np.power(self.image.frames[self.config.errors], 2) + uncertainty**2)

        ### CALIBRATION ERRORS

        if self.config.uncertainties.add_calibration_error:

            # Get calibration error
            if self.config.uncertainties.calibration_error_type == "abs":

                self.image.frames[self.config.errors] += self.config.uncertainties.calibration_error

            elif self.config.uncertainties.calibration_error_type == "rel":

                # Create frame for the calibration error map
                self.image.frames[self.config.errors] += self.config.uncertainties.calibration_error * self.image.frames[self.config.primary]

    # -----------------------------------------------------------------

    def crop(self):

        """
        This function ...
        :param frame:
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
