#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.truncation Contains the Truncator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...magic.region.list import SkyRegionList
from ...magic.basics.mask import Mask as oldMask
from ...magic.core.mask import Mask as newMask
from .component import TruncationComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.dist_ellipse import distance_ellipse

# -----------------------------------------------------------------

class Truncator(TruncationComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(Truncator, self).__init__(config, interactive)

        # --- Attributes ---

        #self.images = dict()
        self.bad_masks = dict()
        self.padded_masks = dict()

        # The frames and error maps
        self.frames = dict()
        self.error_maps = dict()

        # The truncated images (keys are the different scale factors)
        self.truncated_images = dict()

        # Truncation ellipse
        self.ellipse = None

        # The truncation masks
        self.masks = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the images
        self.load_images()

        # 3. Find the best radius for the truncation
        self.find_radius()

        # 3. Truncate the images
        self.truncate()

        # Create the truncation ellipse
        #self.create_ellipse()

        # 4. Create the masks
        #self.create_masks()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Truncator, self).setup()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info(" ... ")

        # Loop over the test images
        for name in self.dataset.names:

            # Get mask names
            mask_names = self.dataset.masks_in_image(name)

            # Check if 'bad' or 'padded' mask is present
            if "bad" in mask_names or "padded" in mask_names:

                # Load the image
                image = self.dataset.get_image(name)

                # Add the frame
                frame = image.primary
                self.images[name] = frame

                # Add bad mask
                if "bad" in image.masks: self.bad_masks[name] = image.masks.bad

                # Add padded mask
                if "padded" in image.masks: self.padded_masks[name] = image.masks.padded

            # No 'bad' or 'padded' mask
            else: self.images[name] = self.dataset.get_frame(name, masked=False)

    # -----------------------------------------------------------------

    def find_radius(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding the best truncation radius ...")

        # Get the angle
        center = self.disk_ellipse.center  # in sky coordinates
        semimajor = self.disk_ellipse.semimajor
        semiminor = self.disk_ellipse.semiminor
        angle = self.disk_ellipse.angle

        # Detemrine the ratio of semimajor and semiminor
        ratio = semiminor / semimajor

        # Loop over all images
        for name in self.config.image_names:

            # Convert center to pixel coordinates
            center_pix = center.to_pixel(self.images[name].wcs)

            # Create distance-ellipse
            distance_frame = distance_ellipse(self.images[name].shape, center_pix, ratio, angle)

            radius_list = []
            signal_to_noise_list = []
            nbad_list = []

            # Loop over the radii
            min_distance = np.min(distance_frame)
            max_distance = np.max(distance_frame)
            for radius in np.linspace(min_distance,max_distance,num=int(max_distance-min_distance+1),dtype=int,endpoint=True):

                # Make a mask of the pixels corresponding to the current radius
                mask = distance_frame == radius

                # Calculate the mean signal to noise in the pixels
                signal_to_noises = self.frames[name][mask] / self.error_maps[name][mask]

                # Calcalute the mean signal to noise
                signal_to_noise = np.mean(signal_to_noises)

                # Add point
                radius_list.append(radius)
                signal_to_noise_list.append(radius)
                nbad_list.append(radius)

    # -----------------------------------------------------------------

    def truncate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making truncated images ...")

        # Loop over the different scale factors
        #for factor in (self.config.factor_range.linear(self.config.factor_nvalues, as_list=True) + [self.config.best_factor]):
        for factor in self.config.factor_range.linear(self.config.factor_nvalues, as_list=True):

            # Truncate the images with this factor
            truncated = self.make_truncated_images(factor)

            # Add the truncated frame to the dictionary
            self.truncated_images[factor] = truncated

    # -----------------------------------------------------------------

    def make_truncated_images(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        # Inform the user
        log.info("Making truncated image for a scale factor of " + str(factor) + " ...")

        # Truncated images
        truncated_images = dict()

        # Loop over the test images
        for name in self.config.image_names:

            # Create ellipse in image coordinates from ellipse in sky coordinates for this image
            ellipse = self.disk_ellipse.to_pixel(self.images[name].wcs) * factor

            # Create mask from ellipse
            mask = oldMask.from_shape(ellipse, self.images[name].xsize, self.images[name].ysize, invert=True)

            # Add masks of padded and bad pixels
            # for example for WISE, this mask even covers pixels within the elliptical region
            if name in self.bad_masks: mask += self.bad_masks[name]
            if name in self.padded_masks: mask += self.padded_masks[name]

            # Truncate the image
            truncated = self.images[name].copy()
            truncated[mask] = 0.0

            # Add the image
            truncated_images[name] = truncated

        # Return the dictionary of truncated images
        return truncated_images

    # -----------------------------------------------------------------

    def create_ellipse(self):

        """
        This function ...
        :return:
        """

        self.ellipse = self.disk_ellipse * self.config.best_factor

    # -----------------------------------------------------------------

    def create_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info(" ... ")

        # Loop over all images
        for name in self.images:

            # Create ellipse in image coordinates from ellipse in sky coordinates for this image
            ellipse = self.disk_ellipse.to_pixel(self.images[name].wcs) * self.config.best_factor

            # Create mask from ellipse
            mask = newMask(oldMask.from_shape(ellipse, self.images[name].xsize, self.images[name].ysize, invert=True))

            # Add masks of padded and bad pixels
            # for example for WISE, this mask even covers pixels within the elliptical region
            if name in self.bad_masks: mask += newMask(self.bad_masks[name])
            if name in self.padded_masks: mask += newMask(self.padded_masks[name])

            # Add the mask
            self.masks[name] = mask

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write curves of the signal-to-noise and the number of bad pixels
        self.write_curves()

        # Write the truncated images
        self.write_images()

        # Write the truncation ellipse
        self.write_ellipse()

        # Write the truncation masks
        self.write_masks()

        # Write the reference truncation mask
        self.write_reference_mask()

    # -----------------------------------------------------------------

    def write_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing curves ...")

        # Write signal to noise curves
        self.write_signal_to_noise_curves()

        # Write bad pixel curves
        self.write_bad_pixel_curves()

    # -----------------------------------------------------------------

    def write_signal_to_noise_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing curves of the signal-to-noise ...")

    # -----------------------------------------------------------------

    def write_bad_pixel_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing curves of the number of bad pixels ...")

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the truncated images ...")

        for factor in self.truncated_images:

            path = fs.create_directory_in(self.truncation_images_path, str(factor))

            # Save the images
            for name in self.truncated_images[factor]:

                image_path = fs.join(path, name + ".fits")
                self.truncated_images[factor][name].saveto(image_path)

    # -----------------------------------------------------------------

    def write_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the truncation ellipse region ...")

        # Determine the path to the region file
        path = fs.join(self.truncation_path, "ellipse.reg")

        # Write the ellipse region
        region = SkyRegionList()
        region.append(self.ellipse)
        region.saveto(path)

    # -----------------------------------------------------------------

    def write_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks ...")

        # Loop over the masks
        for name in self.masks:

            # Save the mask
            path = fs.join(self.truncation_masks_path, name + ".fits")
            self.masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_reference_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the reference truncation mask ...")

        # Save the mask created for the reference image as a seperate file ("reference.fits")
        self.masks[self.reference_image].saveto(self.reference_mask_path)

# -----------------------------------------------------------------
