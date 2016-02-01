#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.imageimporter Contains the ImageImporter class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant AstroMagic classes and modules
from ..core import Image, Box
from ..basics import Mask, Region

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class ImageImporter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ImageImporter, self).__init__(config, "magic")

        # The image
        self.image = None

        # The path and name of the image
        self.directory_path = None
        self.image_path = None
        self.image_name = None

        # Region
        self.bad_region = None

        # Unit and FWHM
        self.unit = None
        self.fwhm = None

    # -----------------------------------------------------------------

    def run(self, path, bad_region_path=None, unit=None, fwhm=None):

        """
        This function ...
        :param path:
        :param bad_region_path:
        :param unit:
        :param fwhm:
        :return:
        """

        # 1. Call the setup function
        self.setup(path, bad_region_path, unit, fwhm)

        # 2. Load the image
        self.load_image()

        # 3. Set the mask
        self.set_mask()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        self.image = None
        self.directory_path = None
        self.image_path = None
        self.image_name = None
        self.bad_region = None
        self.unit = None
        self.fwhm = None

    # -----------------------------------------------------------------

    def setup(self, path, bad_region_path=None, unit=None, fwhm=None):

        """
        This function ...
        :param path:
        :param bad_region_path:
        :param unit:
        :param fwhm:
        :return:
        """

        # Call the setup function of the base class
        super(ImageImporter, self).setup()

        # Set the image path and name
        self.directory_path = os.path.dirname(path)
        self.image_path = path
        self.image_name = os.path.splitext(os.path.basename(self.image_path))[0]

        # Set the bad region
        if bad_region_path is not None: self.bad_region = Region.from_file(bad_region_path)

        # Set the unit and FWHM
        self.unit = unit
        self.fwhm = fwhm

    # -----------------------------------------------------------------

    def load_image(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Importing image from " + self.image_path + " ...")

        # Open the image
        self.image = Image(self.image_path)

        # Load error frame
        self.load_error_frame()

        # Select the appropriate frames
        self.select_frames()

        # Set the image unit and FWHM
        self.image.set_unit(self.unit)
        self.image.set_fwhm(self.fwhm)

    # -----------------------------------------------------------------

    def load_error_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Looking for an error frame ...")

        # If no error map was found in the FITS file, try to find a seperate FITS file containing error data
        if self.image.frames.errors is None:

            # self.directory_path is the same directory as where the image is located
            error_path = os.path.join(self.directory_path, self.image_name + " error.fits")
            if os.path.isfile(error_path): self.image.load_frames(error_path, 0, "errors", "the error map")

        # Still no errors frame
        if self.image.frames.errors is None: self.log.warning("No error data found for " + self.image_name + ".fits")

    # -----------------------------------------------------------------

    def select_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Selecting the appropriate frames ...")

        # Select the primary and errors frame (if present)
        self.image.deselect_all()
        self.image.frames.primary.select()
        if "errors" in self.image.frames: self.image.frames.primary.select()

    # -----------------------------------------------------------------

    def set_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Creating a mask to cover bad pixels ...")

        # Create a mask for the nans in the primary
        nan_mask = Mask.is_nan(self.image.frames.primary)

        # Sometimes, saturated stars have a few pixels that are nan. In this case, we certainly don't want to ignore these pixels
        # because we want to remove the star and its diffraction spikes. So, we want to remove these little blobs of nan from the nan_mask,
        # and interpolate the image there ... So here we seperate the nan_mask into a new mask without the little blobs, and a mask consisting
        # of only these blobs.

        from ..tools import plotting
        #plotting.plot_box(nan_mask[1200:1300, 700:790], "nan mask")

        #eroded_nan_mask = nan_mask.eroded()
        #plotting.plot_box(eroded_nan_mask[1200:1300, 700:790], "eroded nan mask")

        #blob_mask = nan_mask - eroded_nan_mask
        #plotting.plot_box(blob_mask[1200:1300, 700:790], "blob mask")
        #plotting.plot_box(blob_mask, "blob mask")

        from photutils import detect_sources
        segments = detect_sources(nan_mask.astype(float), 0.1, 1).data

        # Check where the nan_mask hits the boundary
        hits_boundary, where = nan_mask.hits_boundary(where=True)

        blob_mask = nan_mask.copy()
        if hits_boundary:

            #print(where[0])

            for pixel in where:

                index = segments[pixel.y, pixel.x]

                #plotting.plot_box(segments == index, "segments == index")

                blob_mask[segments == index] = False
                segments[segments == index] = False

        #plotting.plot_box(segments, "segmentation map")
        #plotting.plot_box(blob_mask[1200:1300, 700:790], "blob mask")

        # Interpolate the frame over the blobs

        # Get a list of contours for the blobs (the oversaturated pixels)
        from ..analysis import sources
        contours = sources.find_contours(segments, segments, 10.0)

        #print(contours)

        # Create a file
        #f = open(os.path.join(self.directory_path, "oversaturated.reg"), 'w')
        # Initialize the region string
        #print("# Region file format: DS9 version 4.1", file=f)
        #for contour in contours:
        #    print("image;ellipse({},{},{},{})".format(contour.center.x+1, contour.center.y+1, contour.radius.x, contour.radius.y), file=f)
        #f.close()

        for contour in contours:

            # Create source object
            #source = Source.from_ellipse(self.image.frames.primary, contour, 1.5)
            #source.plot()

            cutout = Box.from_ellipse(self.image.frames.primary, contour)

            #cutout.plot()

            cutout_segment = segments[cutout.y_slice, cutout.x_slice]

            #plotting.plot_box(cutout_segment)

            import numpy as np
            cutout[np.isnan(cutout)] = 0.0

            where_is_cutout_segment = cutout_segment.astype(bool)

            interpolated_box = cutout.interpolated(where_is_cutout_segment, "local_mean")

            #plotting.plot_box(interpolated_box)

            # Replace frame pixels
            self.image.frames.primary[cutout.y_slice, cutout.x_slice][where_is_cutout_segment] = interpolated_box[where_is_cutout_segment]
            nan_mask[cutout.y_slice, cutout.x_slice][where_is_cutout_segment] = False

        # Add the mask
        self.image.add_mask(nan_mask, "bad")

        # Add the bad mask
        if self.bad_region is not None:

            bad_mask = Mask.from_region(self.bad_region, self.image.frames.primary.shape)
            self.image.masks.bad += bad_mask

# -----------------------------------------------------------------
