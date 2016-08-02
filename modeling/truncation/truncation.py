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
from ...magic.basics.skyregion import SkyRegion
from ...magic.basics.mask import Mask as oldMask
from ...magic.core.mask import Mask as newMask
from .component import TruncationComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

# -----------------------------------------------------------------

# TODO: also crop the FITS files to the bounding box of the disk ellipse?

# -----------------------------------------------------------------

class Truncator(TruncationComponent):
    
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
        super(Truncator, self).__init__(config)

        # --- Attributes ---

        self.images = dict()
        self.bad_masks = dict()
        self.padded_masks = dict()

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

        # 3. Truncate the images
        self.truncate()

        # Create the truncation ellipse
        self.create_ellipse()

        # 4. Create the masks
        self.create_masks()

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
                frame = image.frames.primary
                self.images[name] = frame

                # Add bad mask
                if "bad" in image.masks: self.bad_masks[name] = image.masks.bad

                # Add padded mask
                if "padded" in image.masks: self.padded_masks[name] = image.masks.padded

            # No 'bad' or 'padded' mask
            else: self.images[name] = self.dataset.get_frame(name, masked=False)

    # -----------------------------------------------------------------

    def truncate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info(" ...")

        # Loop over the different scale factors
        for factor in (self.config.factor_range.linear(self.config.factor_nvalues, as_list=True) + [self.config.best_factor]):

            # Calculate the corrected 24 micron image
            truncated = self.make_truncated_images(factor)

            # Add the attenuation map to the dictionary
            self.truncated_images[factor] = truncated

    # -----------------------------------------------------------------

    def make_truncated_images(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        # Inform the user
        log.info("Making ...")

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

    def truncate_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the images ...")

        # Loop over all the images
        for image in self.images:

            # Inform the user
            log.debug("Truncating the " + image.name + " image ...")

            # Inform the user
            log.debug("Creating mask of truncated pixels ...")

            # Create ellipse in image coordinates from ellipse in sky coordinates for this image
            ellipse = self.disk_ellipse.to_pixel(image.wcs)

            # Create mask from ellipse
            inverted_mask = Mask.from_shape(ellipse, image.xsize, image.ysize, invert=True)

            # Get mask of padded (and bad) pixels, for example for WISE, this mask even covers pixels within the elliptical region
            mask = inverted_mask + image.masks.bad
            if "padded" in image.masks: mask += image.masks.padded

            # Truncate the primary frame
            image.frames.primary[mask] = 0.0

            # Truncate the errors frame
            image.frames.errors[mask] = 0.0
            image.frames.errors[np.isnan(image.frames.errors)] = 0.0

            # Remove all other frames
            image.remove_frames_except(["primary", "errors"])

            # Remove all masks
            image.remove_all_masks()

        # Truncate the bulge and disk images
        mask = Mask.from_shape(self.disk_ellipse.to_pixel(self.disk.wcs), self.disk.xsize, self.disk.ysize, invert=True)
        self.disk[mask] = 0.0
        self.bulge[mask] = 0.0
        self.model[mask] = 0.0

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the truncated images
        self.write_images()

        # Write the truncation ellipse
        self.write_ellipse()

        # Write the truncation masks
        self.write_masks()

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
                self.truncated_images[factor][name].save(image_path)

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
        region = SkyRegion()
        region.append(self.ellipse)
        region.save(path)

    # -----------------------------------------------------------------

    def write_masks(self):

        """
        This function ...
        :return:
        """

        # Loop over the masks
        for name in self.masks:

            # Save the mask
            path = fs.join(self.truncation_masks_path, name + ".fits")
            self.masks[name].save(path)

# -----------------------------------------------------------------
