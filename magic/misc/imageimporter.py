#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.imageimporter Contains the ImageImporter class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..core.image import Image
from ..core.frame import Frame
from ..core.cutout import Cutout
from ..basics.mask import Mask
from ..region.list import PixelRegionList
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ImageImporter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ImageImporter, self).__init__(*args, **kwargs)

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

    def _run(self, **kwargs):

        """
        This function ...
        :param path:
        :param bad_region_path:
        :param unit:
        :param fwhm:
        :param find_error_frame:
        :return:
        """

        # 2. Load the image
        self.load_image(**kwargs)

        # 3. Set the mask of bad pixels
        if "bad" not in self.image.masks: self.set_mask()

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

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get arguments
        path = kwargs.pop("path")
        bad_region_path = kwargs.pop("bad_region_path", None)
        unit = kwargs.pop("unit", None)
        fwhm = kwargs.pop("fwhm", None)

        # Call the setup function of the base class
        super(ImageImporter, self).setup(**kwargs)

        # Set the image path and name
        self.directory_path = fs.directory_of(path)
        self.image_path = path
        self.image_name = fs.strip_extension(fs.name(self.image_path))

        # Set the bad region
        if bad_region_path is not None: self.bad_region = Region.from_file(bad_region_path)

        # Set the unit and FWHM
        self.unit = unit
        self.fwhm = fwhm

    # -----------------------------------------------------------------

    def load_image(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Importing image from " + self.image_path + " ...")

        # Get flag
        find_error_frame = kwargs.pop("find_error_frame", "True")

        # Open the image
        self.image = Image.from_file(self.image_path)

        # Load error frame
        if find_error_frame: self.load_error_frame()

        # Set the image unit and FWHM
        if self.unit is not None: self.image.unit = self.unit
        if self.fwhm is not None: self.image.fwhm = self.fwhm

    # -----------------------------------------------------------------

    def load_error_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for an error frame ...")

        # If no error map was found in the FITS file, try to find a seperate FITS file containing error data
        if "errors" not in self.image.frames:

            # self.directory_path is the same directory as where the image is located
            error_path = fs.join(self.directory_path, self.image_name + "_Error.fits")
            #if os.path.isfile(error_path): self.image.load_frames(error_path, 0, "errors", "the error map") # is now possible, i added the "rebin_to_wcs flag" ...

            # Check if the errors frame exists
            if fs.is_file(error_path):

                # Open the errors frame
                error_frame = Frame.from_file(error_path, name="errors", description="the error map")

                # Check if the shape of the error frame matches the shape of the image
                if self.image.shape != error_frame.shape:

                    # Inform the user
                    log.warning("The error frame does not have the same shape as the image, errors frame will be rebinned")

                    # Check if the unit is a surface brightness unit
                    if error_frame.unit != Unit("MJy/sr"): raise ValueError("Cannot rebin since unit " + str(error_frame.unit) + " is not recognized as a surface brightness unit")

                    # Do the rebinning
                    error_frame = error_frame.rebinned(self.image.primary.wcs)

                # Add the error frame
                self.image.add_frame(error_frame, "errors")

    # -----------------------------------------------------------------

    def set_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating a mask to cover bad pixels ...")

        # Create a mask for the nans in the primary
        nan_mask = Mask.is_nan(self.image.primary)

        # Sometimes, saturated stars have a few pixels that are nan. In this case, we certainly don't want to ignore these pixels
        # because we want to remove the star and its diffraction spikes. So, we want to remove these little blobs of nan from the nan_mask,
        # and interpolate the image there ... So here we seperate the nan_mask into a new mask without the little blobs, and a mask consisting
        # of only these blobs.

        #from ..tools import plotting
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

            indices = set()

            for pixel in where:
                index = segments[pixel.y, pixel.x]
                indices.add(index)

            #print("indices=", indices)

            for index in indices:

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

            cutout = Cutout.from_ellipse(self.image.primary, contour)

            #cutout.plot()

            cutout_segment = segments[cutout.y_slice, cutout.x_slice]

            #plotting.plot_box(cutout_segment)

            import numpy as np
            cutout[np.isnan(cutout)] = 0.0

            where_is_cutout_segment = cutout_segment.astype(bool)

            interpolated_box = cutout.interpolated(where_is_cutout_segment, "local_mean")

            #plotting.plot_box(interpolated_box)

            # Replace frame pixels
            self.image.primary[cutout.y_slice, cutout.x_slice][where_is_cutout_segment] = interpolated_box[where_is_cutout_segment]
            nan_mask[cutout.y_slice, cutout.x_slice][where_is_cutout_segment] = False

        # Add the mask
        self.image.add_mask(nan_mask, "bad")

        # Add the bad mask
        if self.bad_region is not None:

            bad_mask = Mask.from_region(self.bad_region, self.image.xsize, self.image.ysize)
            self.image.masks.bad += bad_mask

# -----------------------------------------------------------------
