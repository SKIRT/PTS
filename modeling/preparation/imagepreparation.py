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
import tempfile
import numpy as np

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.basics.mask import Mask
from ...magic.sources.extractor import SourceExtractor
from ...magic.sky.skysubtractor import SkySubtractor
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from . import unitconversion
from ...core.basics.remote import Remote
from ...core.tools import time, filesystem

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

        # The regions
        self.galaxy_region = None
        self.star_region = None
        self.saturation_region = None
        self.other_region = None

        # The segmentation maps
        self.galaxy_segments = None
        self.star_segments = None
        self.other_segments = None

        # The principal ellipse and saturation region in sky coordinates
        self.principal_ellipse_sky = None
        self.saturation_region_sky = None

    # -----------------------------------------------------------------

    def run(self, image, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments):

        """
        This function ...
        :param image:
        :param galaxy_region:
        :param star_region:
        :param saturation_region:
        :param other_region:
        :param galaxy_segments:
        :param star_segments:
        :param other_segments:
        :return:
        """

        # 1. Call the setup function
        self.setup(image, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments)

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

        # 11. Save the result
        self.save_result()

    # -----------------------------------------------------------------

    def setup(self, image, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments):

        """
        This function ...
        :param image:
        :param galaxy_region:
        :param star_region:
        :param saturation_region:
        :param other_region:
        :param galaxy_segments:
        :param star_segments:
        :param other_segments:
        :return:
        """

        # -- Children --

        # Add extractor and sky subtractor
        self.add_child("extractor", SourceExtractor, self.config.extraction)
        self.add_child("sky_subtractor", SkySubtractor, self.config.sky_subtraction)

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

        # Make local references to the regions
        self.galaxy_region = galaxy_region
        self.star_region = star_region
        self.saturation_region = saturation_region
        self.other_region = other_region

        # Make local references to the segmentation maps
        self.galaxy_segments = galaxy_segments
        self.star_segments = star_segments
        self.other_segments = other_segments

    # -----------------------------------------------------------------

    def calculate_calibration_uncertainties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the calibration uncertainties ...")

        # Add the calibration uncertainty defined in (AB) magnitude
        if "mag" in self.config.uncertainties.calibration_error:

            # -----------------------------------------------------------------

            # Convert the frame into AB magnitudes
            invalid = Mask.is_zero_or_less(self.image.frames.primary)
            ab_frame = unitconversion.jansky_to_ab(self.image.frames.primary)
            # Set infinites to zero
            ab_frame[invalid] = 0.0

            # -----------------------------------------------------------------

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
            #c = a - jansky_frame
            c = a - self.image.frames.primary

            # d = image[Jy] - b[Jy]
            #d = jansky_frame - b
            d = self.image.frames.primary - b

            # ----------------------------------------------------------------- BELOW: if frame was not already in Jy

            # Calibration errors = max(c, d)
            #calibration_errors = np.maximum(c, d)  # element-wise maxima
            #calibration_errors[invalid] = 0.0 # set zero where AB magnitude could not be calculated

            #relative_calibration_errors = calibration_errors / jansky_frame
            #relative_calibration_errors[invalid] = 0.0

            # Check that there are no infinities or nans in the result
            #assert not np.any(np.isinf(relative_calibration_errors)) and not np.any(np.isnan(relative_calibration_errors))

            # The actual calibration errors in the same unit as the data
            #calibration_frame = self.image.frames.primary * relative_calibration_errors

            # -----------------------------------------------------------------

            calibration_frame = np.maximum(c, d) # element-wise maxima
            calibration_frame[invalid] = 0.0 # set zero where AB magnitude could not be calculated

            new_invalid = Mask.is_nan(calibration_frame) + Mask.is_inf(calibration_frame)
            calibration_frame[new_invalid] = 0.0

            # Check that there are no infinities or nans in the result
            assert not np.any(np.isinf(calibration_frame)) and not np.any(np.isnan(calibration_frame))

            # -----------------------------------------------------------------

        # The calibration uncertainty is expressed in a percentage (from the flux values)
        elif "%" in self.config.uncertainties.calibration_error:

            # Calculate calibration errors with percentage
            fraction = float(self.config.uncertainties.calibration_error.split("%")[0]) * 0.01
            calibration_frame = self.image.frames.primary * fraction

        # Unrecognized calibration error (not a magnitude, not a percentage)
        else: raise ValueError("Unrecognized calibration error")

        # Add the calibration frame
        self.image.add_frame(calibration_frame, "calibration_errors")

        # Inform the user
        log.success("Calibration uncertainties calculated")

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the sources ...")

        # Run the extractor
        self.extractor.run(self.image.frames.primary, self.galaxy_region, self.star_region, self.saturation_region,
                           self.other_region, self.galaxy_segments, self.star_segments, self.other_segments)

        # Add the sources mask to the image
        self.image.add_mask(self.extractor.mask, "sources")

        # Get the principal ellipse in sky coordinates
        self.principal_ellipse_sky = self.extractor.principal_ellipse.to_sky(self.image.wcs)

        # Get the saturation region in sky coordinates
        self.saturation_region_sky = self.saturation_region.to_sky(self.image.wcs) if self.saturation_region is not None else None

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("extracted.fits")

        # Inform the user
        log.success("Sources extracted")

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting for galactic extinction ...")

        # Correct the primary frame for galactic extinction
        self.image.frames.primary *= 10**(0.4 * self.config.attenuation)

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("corrected_for_extinction.fits")

        # Inform the user
        log.success("Galactic extinction correction finished")

    # -----------------------------------------------------------------

    def convert_unit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting to surface brightness units ...")

        # Run the unit conversion
        #self.unit_converter.run(self.image)

        assert self.image.unit == Unit("Jy/pix")

        # Get pixelscale
        pixelscale = self.image.xy_average_pixelscale

        # Conversion from Jy / pix to MJy / pix
        conversion_factor = 1e-6

        # Conversion from MJy / pix to MJy / sr
        conversion_factor *= (1.0/pixelscale**2).to("pix2/sr").value

        # Multiply the image with the conversion factor
        self.image *= conversion_factor

        # Set the new unit
        self.image.unit = Unit("MJy/sr")

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("converted_unit.fits")

        # Inform the user
        log.success("Units converted")

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Convolving the image with kernel " + self.config.convolution.kernel_path + " ...")

        # Check whether the convolution has to be performed remotely
        if self.config.convolution.remote:

            # Inform the user
            log.info("Convolution will be performed remotely on host nancy")

            # Debugging
            log.debug("Logging in to remote host ...")

            # Create a remote instance for Nancy
            remote = Remote()
            remote.setup("nancy")

            # If intermediate results are not saved, explicitly save the image as it is now
            if not self.config.write_steps: self.write_intermediate_result("converted_unit.fits")

            # Debugging
            log.debug("Creating temporary directory remotely ...")

            # Create a temporary directory to do the convolution
            remote_home_directory = remote.home_directory
            remote_temp_path = filesystem.join(remote_home_directory, time.unique_name("convolution"))
            remote.create_directory(remote_temp_path)

            # Debugging
            log.debug("Uploading the kernel to the remote directory ...")

            # Upload the kernel FITS file to the remote directory
            remote_kernel_path = filesystem.join(remote_temp_path, "kernel.fits")
            remote.upload(self.config.convolution.kernel_path, remote_kernel_path)

            # Debugging
            log.debug("Uploading the image to the remote directory ...")

            # Upload the image to the remote directory
            local_image_path = self.full_output_path("converted_unit.fits")
            remote_image_path = filesystem.join(remote_image_path, "image.fits")
            remote.upload(local_image_path, remote_image_path)

            # Debugging
            log.debug("Creating a python script to perform the convolution remotely ...")

            # Create a python script that does the convolution
            script_file = tempfile.NamedTemporaryFile()
            local_script_path = script_file.name
            with open(local_script_path, 'w') as script_file:

                script_file.write("#!/usr/bin/env python\n")
                script_file.write("# -*- coding: utf8 -*-\n")
                script_file.write("\n")
                script_file.write("# Import astronomical modules\n")
                script_file.write("from astropy.units import Unit\n")
                script_file.write("\n")
                script_file.write("# Import the relevant PTS classes and modules\n")
                script_file.write("from pts.magic.core.frame import Frame\n")
                script_file.write("from pts.magic.core.image import Image\n")
                script_file.write("\n")
                script_file.write("# Open the kernel frame\n")
                script_file.write("kernel = Frame.from_file(" + remote_kernel_path + ")\n")
                script_file.write("\n")
                script_file.write("# Set the FWHM of the kernel\n")
                script_file.write("fwmh = " + str(self.config.convolution.kernel_fwhm.to("arcsec")) + " * Unit('arcsec')\n")
                script_file.write("kernel.fwhm = fwhm\n")
                script_file.write("\n")
                script_file.write("# Open the image\n")
                script_file.write("image = Image.from_file(" + remote_image_path + ")\n")
                script_file.write("\n")
                script_file.write("# Do the convolution")
                script_file.write("image.convolve(kernel)\n")
                script_file.write("\n")
                script_file.write("# Save the image\n")
                script_file.write("image.save(" + remote_image_path + ")\n")

            # Debugging
            log.debug("Uploading the python script ...")

            # Upload the script file
            remote_script_path = filesystem.join(remote_temp_path, "convolve.py")
            remote.upload(local_script_path, remote_script_path)

            # Debugging
            log.debug("Executing the script remotely ...")

            # Execute the script file remotely
            remote.execute("python " + remote_script_path, output=False)

            # Debugging
            log.debug("Downloading the result ...")

            # Download the resulting FITS file (the convolved image)
            local_result_path = self.full_output_path("convolved.fits")
            remote.download(remote_image_path, local_result_path)

            # Load the result
            self.image = Image.from_file(local_result_path)

            # If the results of intermediate steps should not be kept, remove the (downloaded) convolved image file
            if not self.config.write_steps: filesystem.remove_file(local_result_path)

        else:

            # Open the kernel frame
            kernel = Frame.from_file(self.config.convolution.kernel_path)

            # Set the kernel FWHM
            kernel.fwhm = self.config.convolution.kernel_fwhm

            # Convolve the image (the primary and errors frame)
            self.image.convolve(kernel)

            # Save convolved frame
            if self.config.write_steps: self.write_intermediate_result("convolved.fits")

        # Inform the user
        log.success("Convolution finished")

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

        # Inform the user
        log.success("Rebinning finished")

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the sky ...")

        # Write out the principal_sky_ellipse and saturation_region in sky coordinates
        #principal_sky_region = SkyRegion()
        #principal_sky_region.append(self.principal_ellipse_skycoord)
        #principal_sky_region_path = self.full_output_path("principal_skycoord.reg")
        #principal_sky_region.save(principal_sky_region_path)

        #saturation_sky_region_path = self.full_output_path("saturation_skycoord.reg")
        #self.extractor.star_extractor.saturation_region.save(saturation_sky_region_path)

        # Convert the principal ellipse in sky coordinates into pixel coordinates
        principal_ellipse = self.principal_ellipse_sky.to_pixel(self.image.wcs)

        # Conver the saturation region in sky coordinate into pixel coordinates
        if self.saturation_region_sky is not None:
            saturation_region = self.saturation_region_sky.to_pixel(self.image.wcs)
        else: saturation_region = None

        # Create the 'extra' mask (bad and padded pixels)
        extra_mask = None
        if "bad" in self.image.masks:
            # Combine with padded mask or just use bad mask
            if "padded" in self.image.masks: extra_mask = self.image.masks.bad + self.image.masks.padded
            else: extra_mask = self.image.masks.bad
        elif "padded" in self.image.masks: extra_mask = self.image.masks.padded # just use padded mask

        # Run the sky subtractor
        self.sky_subtractor.run(self.image.frames.primary, principal_ellipse, self.image.masks.sources, extra_mask, saturation_region)

        # Add the sky map to the image
        self.image.add_frame(self.sky_subtractor.sky, "sky")
        self.image.add_frame(self.sky_subtractor.phot_sky, "phot_sky")
        self.image.add_frame(self.sky_subtractor.phot_rms, "phot_rms")

        # Write intermediate result
        if self.config.write_steps: self.write_intermediate_result("sky_subtracted.fits")

        # Inform the user
        log.success("Sky subtracted")

    # -----------------------------------------------------------------

    def set_uncertainties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the uncertainties ...")

        sky_mask = self.image.masks.sky

        photutils_estimated_sky = np.ma.masked_array(np.asarray(self.image.frames.phot_sky), mask=sky_mask)
        photutils_rms = np.ma.masked_array(np.asarray(self.image.frames.phot_rms), mask=sky_mask)

        #plotting.plot_box(photutils_rms, title="masked rms")

        # LARGE SCALE VARIATIONS = STANDARD DEVIATION OF ESTIMATED SKY PIXELS
        large_scale_variations_error = photutils_estimated_sky.std()

        # PIXEL TO PIXEL NOISE = MEAN OF RMS
        pixel_to_pixel_noise = np.ma.mean(photutils_rms) # Mean pixel-by-pixel variations

        # Add all the errors quadratically: existing error map + large scale variations + pixel to pixel noise + calibration error
        if "errors" in self.image.frames:

            self.image.frames.errors = np.sqrt(self.image.frames.errors**2 + large_scale_variations_error**2 + \
                                               + pixel_to_pixel_noise**2 + self.image.frames.calibration_errors**2)

        else:

            errors = np.sqrt(large_scale_variations_error**2 + pixel_to_pixel_noise**2 + self.image.frames.calibration_errors**2)
            self.image.add_frame(errors, "errors")

        # Inform the user
        log.success("Uncertainties calculated")

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

    def save_result(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path("result.fits")

        # Inform the user
        log.info("Writing resulting image to " + path + " ...")

        # Write out the resulting image
        self.image.save(path, origin=self.name)

    # -----------------------------------------------------------------

    def write_intermediate_result(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(path)

        # Inform the user
        log.info("Writing intermediate result to " + path + " ...")

        # Write out the image
        self.image.save(path, origin=self.name)

# -----------------------------------------------------------------
