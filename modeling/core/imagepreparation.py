#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.imagepreparation Contains the ImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant AstroMagic classes and modules
from ...magic.core import Frame
from ...magic.basics import Region
from ...magic.extract import Extractor
from ...magic.subtract import SkySubtractor
from ...magic.tools import regions, cropping

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from .unitconversion import UnitConverter

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

        # 2. Create the mask for bad pixels
        self.create_mask()

        # 3. Extract stars and galaxies from the image
        if self.config.extract_sources: self.extract_sources()

        # 4. If requested, subtract the sky
        if self.config.subtract_sky: self.subtract_sky()

        # 5. If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # 6. If requested, convert the unit
        if self.config.convert_unit: self.convert_unit()

        # 7. If requested, convolve
        #if self.config.convolve: self.convolve()

        # 8. If requested, rebin
        #if self.config.rebin: self.rebin()

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

            self.extractor.config.write_result = True
            self.extractor.config.writing.result_path = "extracted.fits"

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

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # Run the sky extraction
        self.sky_subtractor.run(self.image, self.extractor.galaxy_extractor.principal_ellipse, self.extractor.star_extractor.saturation_contours)

        # Print the statistics of the sky frame
        self.log.info("Mean sky level = " + str(self.sky_subtractor.mean))
        self.log.info("Median sky level = " + str(self.sky_subtractor.median))
        self.log.info("Standard deviation of sky = " + str(self.sky_subtractor.stddev))

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("sky_subtracted.fits")

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Correcting image for galactic extinction ...")

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
        self.log.info("Converting image to surface brightness units ...")

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
        self.log.info("Convolving the image with kernel " + self.config.convolution.kernel_path + " ...")

        # Open the kernel frame
        kernel = Frame.from_file(self.config.convolution.kernel_path)

        # Convolve the image (the primary and errors frame)
        self.image.convolve(kernel)

        # Save convolved frame
        path = self.full_output_path("convolved.fits")
        if self.config.write_steps: self.image.frames["primary"].save(path)

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Rebinning the image to the pixel grid of " + self.config.rebinning.rebin_to + " ...")

        # Open the reference frame
        reference = Frame.from_file(self.config.rebinning.rebin_to)

        # Rebin the image (the primary and errors frame)
        self.image.rebin(reference)

        # Save rebinned frame
        path = self.full_output_path("rebinned.fits")
        if self.config.write_steps: self.image.frames["primary"].save(path)

    # -----------------------------------------------------------------

    def set_uncertainties(self):

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
            pixelscale = self.image.frames[self.config.primary].pixelscale
            unit = self.image.frames[self.config.primary].unit
            filter = self.image.frames[self.config.primary].filter
            selected = True

            # Create the errors frame (select it)
            frame = Frame(np.full(self.image.frames[self.config.primary].shape, uncertainty), wcs, pixelscale, None, selected, unit, self.config.errors, filter)

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
        x_min = self.config.cropping.limits[0]
        x_max = self.config.cropping.limits[1]
        y_min = self.config.cropping.limits[2]
        y_max = self.config.cropping.limits[3]

        # Crop the image (the primary and errors frame)
        self.image.crop(x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing ...")

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
        self.log.info("Writing resulting image to " + path + " ...")

        # Write out the resulting image
        self.image.save(path)

    # -----------------------------------------------------------------

    def write_intermediate_result(self, path):

        """
        This function ...
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(path)

        # Inform the user
        self.log.info("Writing intermediate result to " + path + " ...")

        # Write out the image
        self.image.save(path)

# -----------------------------------------------------------------
