#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.skysubtractor Contains the SkySubtractor class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import astronomical modules
from photutils.background import Background

# Import the relevant AstroMagic classes and modules
from ..core.frame import Frame
from ..basics.mask import Mask
from ..tools import interpolation, plotting, statistics, fitting

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log

# -----------------------------------------------------------------

class SkySubtractor(Configurable):

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
        super(SkySubtractor, self).__init__(config, "magic")

        # -- Attributes --

        # The image
        self.image = None

        # The galaxy and saturation region
        self.principal_ellipse = None
        self.saturation_region = None

        # The output mask (combined input + bad mask + galaxy annulus mask + expanded saturation mask + sigma-clipping mask)
        self.mask = None

        # The estimated sky frame
        self.sky = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new SkySubtractor instance
        if arguments.config is not None: subtractor = cls(arguments.config)
        elif arguments.settings is not None: subtractor = cls(arguments.settings)
        else: subtractor = cls()

        # Debug mode
        if arguments.debug:

            subtractor.config.logging.level = "DEBUG"
            subtractor.config.logging.cascade = True

        # Report file
        if arguments.report: subtractor.config.logging.path = "log.txt"

        # Set the input and output path
        if arguments.input_path is not None: subtractor.config.input_path = arguments.input_path
        if arguments.output_path is not None: subtractor.config.output_path = arguments.output_path

        # Return the new instance
        return subtractor

    # -----------------------------------------------------------------

    def run(self, image, principal_ellipse, saturation_region=None):

        """
        This function ...
        :param image:
        :param principal_ellipse:
        :param saturation_region:
        :return:
        """

        # 1. Call the setup function
        self.setup(image, principal_ellipse, saturation_region)

        # 2. Create mask
        self.create_mask()

        # 2. Perform sigma-clipping
        if self.config.sigma_clip_mask: self.sigma_clip()

        # 3. Estimate the sky
        self.estimate()

        # 4. If requested, subtract the sky
        if self.config.subtract: self.subtract()

        # Update the image mask
        self.update_mask()

        # Add the sky frame
        self.add_sky_frame()

        # Set zero outside
        #self.set_zero_outside()

        # 5. Write out the results
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the sky extractor ...")

        # Set all attributes to None
        self.image = None
        self.principal_ellipse = None
        self.saturation_region = None
        self.mask = None
        self.sky = None

    # -----------------------------------------------------------------

    def setup(self, image, principal_ellipse, saturation_region=None):

        """
        This function ...
        :param image:
        :param principal_ellipse:
        :param saturation_region:
        :return:
        """

        # Call the setup function of the base class
        super(SkySubtractor, self).setup()

        # Set the paths to the resulting frame and the total mask
        self.config.write_result = True
        self.config.write_masked_frame = True
        self.config.writing.result_path = "subtracted.fits"
        self.config.writing.masked_frame_path = "masked_sky_frame.fits"

        # Make a local reference to the image
        self.image = image

        # Convert the principal ellipse and saturation region to image coordinates
        self.principal_ellipse = principal_ellipse.to_pixel(self.image.wcs)
        self.saturation_region = saturation_region.to_pixel(self.image.wcs) if saturation_region is not None else None

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Create a mask from the ellipse
        annulus_outer_factor = 3.0
        annulus_inner_factor = 1.0
        annulus_mask = Mask.from_shape(self.principal_ellipse * annulus_outer_factor, self.image.xsize, self.image.ysize).inverse() + \
                       Mask.from_shape(self.principal_ellipse * annulus_inner_factor, self.image.xsize, self.image.ysize)

        # Set the mask, make a copy of the input mask initially
        self.mask = self.image.masks.sources + annulus_mask

        # Add the bad pixels mask
        if "bad" in self.image.masks: self.mask += self.image.masks.bad

        # Check whether saturation contours are defined
        if self.saturation_region is not None:

            # Expand all contours
            expanded_region = self.saturation_region * 1.5

            # Create the saturation mask
            saturation_mask = expanded_region.to_mask(self.image.xsize, self.image.ysize)
            self.mask += saturation_mask

        # Add the mask of padded pixels (during rebinning)
        if "padded" in self.image.masks: self.mask += self.image.masks.padded

    # -----------------------------------------------------------------

    def sigma_clip(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing sigma-clipping on the pixel values ...")

        ### TEMPORARY: WRITE OUT MASK BEFORE CLIPPING

        # Create a frame where the objects are masked
        #frame = copy.deepcopy(self.image.frames.primary)
        #frame[self.mask] = float(self.config.writing.mask_value)

        # Save the masked frame
        #frame.save("masked_sky_frame_notclipped.fits")

        ###

        # Create the sigma-clipped mask
        self.mask = statistics.sigma_clip_mask(self.image.frames.primary, self.config.sigma_clipping.sigma_level, self.mask)

    # -----------------------------------------------------------------

    def estimate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky by using " + self.config.estimation.method + " ...")

        # If the mean sky level should be used
        if self.config.estimation.method == "mean":

            # Create a frame filled with the mean value
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            self.sky = Frame(np.full(self.image.shape, self.mean),
                             wcs=self.image.wcs,
                             name="sky",
                             description="estimated sky",
                             unit=self.image.unit,
                             zero_point=self.image.frames.primary.zero_point,
                             filter=self.image.filter,
                             sky_subtracted=False,
                             fwhm=self.image.fwhm)

        # If the median sky level should be used
        elif self.config.estimation.method == "median":

            # Create a frame filled with the median value
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            self.sky = Frame(np.full(self.image.shape, self.median),
                             wcs=self.image.wcs,
                             name="sky",
                             description="estimated sky",
                             unit=self.image.unit,
                             zero_point=self.image.frames.primary.zero_point,
                             filter=self.image.filter,
                             sky_subtracted=False,
                             fwhm=self.image.fwhm)

        # If the sky should be estimated by using low-resolution interpolation
        elif self.config.estimation.method == "low-res-interpolation":

            # Estimate the sky
            data = interpolation.low_res_interpolation(self.image.frames.primary, self.config.estimation.downsample_factor, self.mask)

            # Create sky map
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            self.sky = Frame(data,
                             wcs=self.image.wcs,
                             name="sky",
                             description="estimated sky",
                             unit=self.image.unit,
                             zero_point=self.image.frames.primary.zero_point,
                             filter=self.image.filter,
                             sky_subtracted=False,
                             fwhm=self.image.fwhm)

        # if the sky should be approximated by a polynomial function
        elif self.config.estimation.method == "polynomial":

            polynomial = fitting.fit_polynomial(self.image.frames.primary, 3, mask=self.mask)

            # Evaluate the polynomial
            data = fitting.evaluate_model(polynomial, 0, self.image.frames.primary.xsize, 0, self.image.frames.primary.ysize)

            plotting.plot_box(data, title="background")

            # Create sky map
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            self.sky = Frame(data,
                             wcs=self.image.wcs,
                             name="sky",
                             description="estimated sky",
                             unit=self.image.unit,
                             zero_point=self.image.frames.primary.zero_point,
                             filter=self.image.filter,
                             sky_subtracted=False,
                             fwhm=self.image.fwhm)

        elif self.config.estimation.method == "photutils":

            #bkg = Background(self.image.frames.primary, (20, 20), filter_shape=(3, 3), method='median', mask=self.mask)

            #method = "sextractor" # default
            method = "mean"

            bkg = Background(self.image.frames.primary, (20, 20), filter_shape=(3, 3), filter_threshold=None, mask=self.mask,
                      method=method, backfunc=None, interp_order=3, sigclip_sigma=3.0, sigclip_iters=10)

            plotting.plot_box(bkg.background, title="background")
            plotting.plot_box(bkg.background_low_res, title="low-res background")
            plotting.plot_box(bkg.background_rms, title="background rms")
            plotting.plot_box(bkg.background_rms_low_res, title="low-res rms")

            print("background median = ", bkg.background_median)
            print("background rms median = ", bkg.background_rms_median)

            # Create sky map
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            self.sky = Frame(bkg.background,
                             wcs=self.image.wcs,
                             name="sky",
                             description="estimated sky",
                             unit=self.image.unit,
                             zero_point=self.image.frames.primary.zero_point,
                             filter=self.image.filter,
                             sky_subtracted=False,
                             fwhm=self.image.fwhm)

        elif self.config.estimation.method == "hybrid":

            bkg = Background(self.image.frames.primary, (50, 50), filter_shape=(3, 3), filter_threshold=None, mask=self.mask,
                      method="sextractor", backfunc=None, interp_order=3, sigclip_sigma=3.0, sigclip_iters=10)

            # Masked background
            masked_background = np.ma.masked_array(bkg.background, mask=self.mask)
            #plotting.plot_box(masked_background, title="masked background")

            mean_sky = np.ma.mean(masked_background)
            median_sky = np.median(masked_background.compressed())

            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            phot_sky_frame = Frame(bkg.background,
                                   wcs=self.image.wcs,
                                   name="phot_sky",
                                   description="photutils background",
                                   unit=self.image.unit,
                                   zero_point=self.image.frames.primary.zero_point,
                                   filter=self.image.filter,
                                   sky_subtracted=False,
                                   fwhm=self.image.fwhm)
            self.image.add_frame(phot_sky_frame, "phot_sky")

            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            phot_rms_frame = Frame(bkg.background_rms,
                                   wcs=self.image.wcs,
                                   name="phot_rms",
                                   description="photutils rms",
                                   unit=self.image.unit,
                                   zero_point=self.image.frames.primary.zero_point,
                                   filter=self.image.filter,
                                   sky_subtracted=False,
                                   fwhm=self.image.fwhm)
            self.image.add_frame(phot_rms_frame, "phot_rms")

            # Create sky map of median sky level
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            self.sky = Frame(np.full(self.image.shape, median_sky),
                             wcs=self.image.wcs,
                             name="sky",
                             description="estimated sky",
                             unit=self.image.unit,
                             zero_point=self.image.frames.primary.zero_point,
                             filter=self.image.filter,
                             sky_subtracted=False,
                             fwhm=self.image.fwhm)

        # Unkown estimation method
        else: raise ValueError("Unkown sky estimation method")

    # -----------------------------------------------------------------

    def subtract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the sky from the frame ...")

        # Subtract the estimated sky from the image frame
        self.image.frames.primary -= self.sky

    # -----------------------------------------------------------------

    def update_mask(self):

        """
        This function ...
        :return:
        """

        self.image.add_mask(self.mask, "sky")

    # -----------------------------------------------------------------

    def add_sky_frame(self):

        """
        This function ...
        :return:
        """

        self.image.add_frame(self.sky, "sky")

    # -----------------------------------------------------------------

    def set_zero_outside(self):

        """
        This function ...
        :return:
        """

        # Create a mask from the principal galaxy region
        annulus_outer_factor = 1.2
        mask = Mask.from_shape(self.principal_ellipse * annulus_outer_factor, self.image.xsize, self.image.ysize).inverse()

        # Set the primary frame zero outside the principal ellipse
        self.image.frames.primary[mask] = 0.0

        # Add mask
        self.image.add_mask(mask, "outside")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write out the result
        #if self.config.write_result: self.write_result()

        # Write out a histogram of the sky pixels
        #if self.config.write_histogram: self.write_histogram()

        # If requested, write out the masked frame
        if self.config.write_masked_frame: self.write_masked_frame()

        # If requested, write out the sky statistics
        if self.config.write_statistics: self.write_statistics()

    # -----------------------------------------------------------------

    def write_result(self, header=None):

        """
        This function ...
        :param header:
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(self.config.writing.result_path)

        # Inform the user
        log.info("Writing resulting frame to " + path + " ...")

        # Write out the resulting frame
        self.image.frames.primary.save(path, header)

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the masked frame file
        path = self.full_output_path(self.config.writing.masked_frame_path)

        # Inform the user
        log.info("Writing masked frame to " + path + " ...")

        # Create a frame where the objects are masked
        frame = self.image.frames.primary.copy()
        frame[self.mask] = float(self.config.writing.mask_value)

        # Write out the masked frame
        frame.save(path)

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the statistics file
        path = self.full_output_path(self.config.writing.statistics_path)

        # Inform the user
        log.info("Writing statistics to " + path + " ...")

    # -----------------------------------------------------------------

    def write_histogram(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing sky histogram to " + self.config.writing.histogram_path +  " ...")

        # Create a masked array
        masked = np.ma.masked_array(self.image.frames.primary, mask=self.mask)
        masked_clipped = np.ma.masked_array(self.image.frames.primary, mask=self.clipped_mask)

        # Create a figure
        fig = plt.figure()

        min = self.mean - 4.0 * self.stddev
        max = self.mean + 4.0 * self.stddev

        # Plot the histograms
        #b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
        plt.subplot(211)
        plt.hist(masked.compressed(), 200, range=(min,max), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='not clipped')
        if self.config.histogram.log_scale: plt.semilogy()

        plt.subplot(212)
        plt.hist(masked_clipped.compressed(), 200, range=(min,max), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='clipped')
        if self.config.histogram.log_scale: plt.semilogy()

        # Save the figure
        plt.savefig(self.config.writing.histogram_path, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    # -----------------------------------------------------------------

    @property
    def mean(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped mean
        return np.ma.mean(np.ma.masked_array(self.image.frames.primary, mask=self.mask))

    # -----------------------------------------------------------------

    @property
    def median(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped median
        return np.median(np.ma.masked_array(self.image.frames.primary, mask=self.mask).compressed())

    # -----------------------------------------------------------------

    @property
    def stddev(self):

        """
        This function ...
        :return:
        """

        # Return the standard deviation of the sigma-clipped frame
        return np.ma.masked_array(self.image.frames.primary, mask=self.mask).std()

# -----------------------------------------------------------------
