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
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...magic.basics.skyregion import SkyRegion
from ...magic.basics.mask import Mask
from .component import TruncationComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

# -----------------------------------------------------------------

# TODO: also crop the FITS files to the bounding box of the disk ellipse

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

        # The list of images
        self.images = []

        # The disk ellipse
        self.disk_ellipse = None

        # The disk and bulge frames
        self.disk = None
        self.bulge = None
        self.model = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new PhotoMeter instance
        truncator = cls(arguments.config)

        # Set the input and output path
        truncator.config.path = arguments.path

        # Return the new instance
        return truncator

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the prepared images
        self.load_images()

        # 3. Get the disk ellipse
        self.get_disk_ellipse()

        # 4. Truncate the images
        self.truncate()

        # 7. Writing
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
        log.info("Loading the images ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Look for a file called 'result.fits'
            image_path = fs.join(directory_path, "result.fits")
            if not fs.is_file(image_path):
                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # If the truncated image is already present, skip it
            truncated_path = fs.join(self.truncation_path, directory_name + ".fits")
            if fs.is_file(truncated_path): continue

            # Open the prepared image
            image = Image.from_file(image_path)

            # Set the image name
            image.name = directory_name

            # Add the image to the list
            self.images.append(image)

        # Load the disk image
        disk_path = fs.join(self.components_path, "disk.fits")
        self.disk = Frame.from_file(disk_path)

        # Load the bulge image
        bulge_path = fs.join(self.components_path, "bulge.fits")
        self.bulge = Frame.from_file(bulge_path)

        # Load the model image
        model_path = fs.join(self.components_path, "model.fits")
        self.model = Frame.from_file(model_path)

    # -----------------------------------------------------------------

    def get_disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Get the path to the disk region
        path = fs.join(self.components_path, "disk.reg")

        # Open the region
        region = SkyRegion.from_file(path)

        # Get ellipse in sky coordinates
        scale_factor = 0.82
        self.disk_ellipse = region[0] * scale_factor

    # -----------------------------------------------------------------

    def truncate(self):

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

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Loop over all images
        for image in self.images:

            # Determine the path to the truncated image
            truncated_path = fs.join(self.truncation_path, image.name + ".fits")

            # Save the image
            image.save(truncated_path)

        # Write the bulge, disk and total model images
        disk_path = fs.join(self.truncation_path, "disk.fits")
        self.disk.save(disk_path)
        bulge_path = fs.join(self.truncation_path, "bulge.fits")
        self.bulge.save(bulge_path)
        model_path = fs.join(self.truncation_path, "model.fits")
        self.model.save(model_path)

# -----------------------------------------------------------------
