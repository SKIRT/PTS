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
import os
import numpy as np

# Import astronomical modules
import astropy.units as u

# Import the relevant AstroMagic classes and modules
from ...magic.core import Frame
from ...magic.basics import Mask, Region
from ...magic.extract import Extractor
from ...magic.subtract import SkySubtractor
from ...magic.tools import regions, cropping

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

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

        # The mask for bad pixels
        self.bad_mask = None

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
        self.extract_sources()

        # 4. If requested, subtract the sky
        if self.config.subtract_sky: self.subtract_sky()

        # 5. If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # 6. If requested, convert the unit
        if self.config.convert_unit: self.convert_unit()

        # 7. If requested, convolve
        if self.config.convolve: self.convolve()

        # 8. If requested, rebin
        if self.config.rebin: self.rebin()

        # 9. If requested, set the uncertainties
        if self.config.set_uncertainties: self.set_uncertainties()

        # 10. If requested, crop
        if self.config.crop: self.crop()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, image):

        """
        This function ...
        :return:
        """

        # -- Children --

        # Add extractor and sky subtractor
        self.add_child("extractor", Extractor, self.config.galaxy_extraction)
        self.add_child("sky_subtractor", SkySubtractor, self.config.sky_extraction)

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

        # Make a local reference to the passed image
        self.image = image

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This fuction ...
        :return:
        """

        # Create a mask for the nans in the primary
        self.bad_mask = Mask.is_nan(self.image.frames["primary"])

        # Make a mask from the extra region
        if self.config.extra_region_path is not None:

            # Replace pixels in the 'extra' region with zeros
            extra_path = self.full_input_path(self.config.extra_region_path)

            extra_region = Region.from_file(extra_path)
            extra_mask = Mask.from_region(extra_region, self.image.frames["primary"].shape)

            # Add the extra mask to the bad pixel mask
            self.bad_mask += extra_mask

        # Set value of bad pixels to zero
        self.image.frames["primary"][self.bad_mask] = 0.0

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return:
        """

        # Run the extractor
        self.extractor.run(self.image.frames[self.config.primary], self.bad_mask)

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # Run the sky extraction
        self.sky_subtractor.run(self.image.frames[self.config.primary], self.extractor.mask + self.bad_mask)

        # Print the statistics of the sky frame
        self.log.info("Mean sky level = " + str(self.skysub.mean))
        self.log.info("Median sky level = " + str(self.skysub.median))
        self.log.info("Standard deviation of sky = " + str(self.skysub.stddev))

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # ...
        if self.image.wavelength < 1.0 * u.Unit("micron"):

            # Correct the primary frame for galactic extinction
            self.image.frames[self.config.primary] *= 10**(0.4 * self.config.attenuation)

        else: self.log.info("No extinction correction for this image")

    # -----------------------------------------------------------------

    def convert_unit(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # TODO: make this function automatic

        # Create a unit object
        unit = u.Unit(self.config.unit_conversion.to_unit)

        # Inform the user
        self.log.info("Converting image to " + str(unit) + ", but not yet automatic")

        # Convert the image to different units (primary and errors frame)
        #self.image.convert_to(unit)

        ####

        # THE GALEX FUV IMAGE
        if self.image.name == "GALEX FUV":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale.value

            # Get the wavelength of the image
            wavelength = self.image.frames.primary.filter.centerwavelength()
            wavelength = (wavelength * u.Unit("micron")).to("AA")

            # Speed of light in Angstrom per seconds
            c = (299792458.0 * u.Unit("m") / u.Unit("s")).to("AA/s")
            spectralfactor = wavelength.value**2 / c.value
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 1.4e-15 * spectralfactor * 1e17 * pixelfactor

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE 2MASS H IMAGE
        elif self.image.name == "2MASS H":

            # Conversion factor to magnitudes
            m0 = 20.6473999

            # Conversion factor back to fluxes
            F0 = 1024.0

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale.value

            # Calculate the conversion factor
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 10e-6 * pixelfactor
            factor *= F0 * np.power(10.0, -m0/2.5)

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE Halpha IMAGE
        elif self.image.name == "Ha":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale.value

            pixelfactor = (206264.806247 / pixelscale)**2

            # The wavelength (in meter)
            wavelength = 0.657894736 * 1e-6

            # The speed of light (in m / s)
            c = 299792458

            # The frequency (in Hz)
            frequency = c / wavelength

            # Calculate the conversion factor
            factor = 1e23 * 1e-6 * pixelfactor / frequency

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE IRAC I1 IMAGE IS ALREADY IN MJY/SR
        elif self.image.name == "IRAC I1": pass

        # THE MIPS 24 IMAGE IS ALREADY IN MJY/SR
        elif self.image.name == "MIPS 24": pass

        # THE PACS 70 IMAGE
        elif self.image.name == "PACS 70":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale.value

            # Calculate the conversion factor
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 1e-6 * pixelfactor

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE PACS 160 IMAGE
        elif self.image.name == "PACS 160":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale.value

            # Calculate the conversion factor
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 1e-6 * pixelfactor

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # UNKOWN IMAGE NAME
        else: raise ValueError("Unkown image: " + self.image.name)

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :param frame:
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
        :param frame:
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
