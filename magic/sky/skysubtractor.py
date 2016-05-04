#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
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

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ..basics.mask import Mask
from ..tools import interpolation, plotting, statistics, fitting
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

        # The image frame
        self.frame = None

        # The mask of sources
        self.sources_mask = None

        # The extra mask
        self.extra_mask = None

        # The principal ellipse
        self.principal_ellipse = None

        # The region of saturated stars
        self.saturation_region = None

        # The output mask (combined input + bad mask + galaxy annulus mask + expanded saturation mask + sigma-clipping mask)
        self.mask = None

        # The estimated sky frame
        self.sky = None

        # Photutils results
        self.phot_sky = None
        self.phot_rms = None

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

        # Return the new instance
        return subtractor

    # -----------------------------------------------------------------

    def run(self, frame, principal_ellipse, sources_mask, extra_mask=None, saturation_region=None):

        """
        This function ...
        :param frame:
        :param principal_ellipse:
        :param sources_mask:
        :param extra_mask:
        :param saturation_region:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, principal_ellipse, sources_mask, extra_mask, saturation_region)

        # 2. Create mask
        self.create_mask()

        # 2. Perform sigma-clipping
        if self.config.sigma_clip_mask: self.sigma_clip()

        # 3. Estimate the sky
        self.estimate()

        # 4. If requested, subtract the sky
        self.subtract()

        # Set zero outside
        #self.set_zero_outside()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the sky subtractor ...")

        # Set all attributes to None
        self.frame = None
        self.principal_ellipse = None
        self.saturation_region = None
        self.mask = None
        self.sky = None

    # -----------------------------------------------------------------

    def setup(self, frame, principal_ellipse, sources_mask, extra_mask=None, saturation_region=None):

        """
        This function ...
        :param frame:
        :param principal_ellipse:
        :param sources_mask:
        :param extra_mask:
        :param saturation_region:
        :return:
        """

        # Call the setup function of the base class
        super(SkySubtractor, self).setup()

        # Make a local reference to the image frame
        self.frame = frame

        # Make a reference to the principal ellipse
        self.principal_ellipse = principal_ellipse

        # Set the masks
        self.sources_mask = sources_mask
        self.extra_mask = extra_mask

        # Set the saturation_region
        self.saturation_region = saturation_region

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Create a mask from the ellipse
        annulus_outer_factor = self.config.annulus_outer_factor
        annulus_inner_factor = self.config.annulus_inner_factor
        annulus_mask = Mask.from_shape(self.principal_ellipse * annulus_outer_factor, self.frame.xsize, self.frame.ysize).inverse() + \
                       Mask.from_shape(self.principal_ellipse * annulus_inner_factor, self.frame.xsize, self.frame.ysize)

        # Set the mask, make a copy of the input mask initially
        self.mask = self.sources_mask + annulus_mask

        # Add the extra mask (if specified)
        if self.extra_mask is not None: self.mask += self.extra_mask

        # Check whether saturation contours are defined
        if self.saturation_region is not None:

            # Expand all contours
            expanded_region = self.saturation_region * 1.5

            # Create the saturation mask
            saturation_mask = expanded_region.to_mask(self.frame.xsize, self.frame.ysize)
            self.mask += saturation_mask

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
        #frame = copy.deepcopy(self.frame)
        #frame[self.mask] = float(self.config.writing.mask_value)

        # Save the masked frame
        #frame.save("masked_sky_frame_notclipped.fits")

        ###

        # Create the sigma-clipped mask
        self.mask = statistics.sigma_clip_mask(self.frame, self.config.sigma_clipping.sigma_level, self.mask)

    # -----------------------------------------------------------------

    def estimate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky ...")

        # Estimate the sky by taking the mean value of all pixels that are not masked
        if self.config.estimation.method == "mean": self.estimate_sky_mean()

        # Estimate the sky by taking the median value of all pixels that are not masked
        elif self.config.estimation.method == "median": self.estimate_sky_median()

        # The sky should be estimated by fitting a polynomial function to the pixels
        elif self.config.estimation.method == "polynomial": self.estimate_sky_polynomial()

        # Use photutils to estimate the sky and sky noise
        elif self.config.estimation.method == "photutils": self.estimate_sky_photutils()

        # Use our own method to estimate the sky and sky noise
        elif self.config.estimation.method == "pts": self.estimate_sky_pts()

        # Unkown sky estimation method
        else: raise ValueError("Unknown sky estimation method")

    # -----------------------------------------------------------------

    def estimate_sky_mean(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky by calculating the mean value of all non-masked pixels ...")

        # Create a frame filled with the mean value
        self.sky = Frame(np.full(self.frame.shape, self.mean),
                         wcs=self.frame.wcs,
                         name="sky",
                         description="estimated sky",
                         unit=self.frame.unit,
                         zero_point=self.frame.zero_point,
                         filter=self.frame.filter,
                         sky_subtracted=False,
                         fwhm=self.frame.fwhm)

    # -----------------------------------------------------------------

    def estimate_sky_median(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky by calculating the median value of all non-masked pixels ...")

        # Create a frame filled with the median value
        self.sky = Frame(np.full(self.frame.shape, self.median),
                         wcs=self.frame.wcs,
                         name="sky",
                         description="estimated sky",
                         unit=self.frame.unit,
                         zero_point=self.frame.zero_point,
                         filter=self.frame.filter,
                         sky_subtracted=False,
                         fwhm=self.frame.fwhm)

    # -----------------------------------------------------------------

    def estimate_sky_polynomial(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky by fitting a polynomial function to all non-masked pixels ...")

        polynomial = fitting.fit_polynomial(self.frame, 3, mask=self.mask)

        # Evaluate the polynomial
        data = fitting.evaluate_model(polynomial, 0, self.frame.xsize, 0, self.frame.ysize)

        plotting.plot_box(data, title="background")

        # Create sky map
        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        self.sky = Frame(data,
                         wcs=self.frame.wcs,
                         name="sky",
                         description="estimated sky",
                         unit=self.frame.unit,
                         zero_point=self.frame.zero_point,
                         filter=self.frame.filter,
                         sky_subtracted=False,
                         fwhm=self.frame.fwhm)

    # -----------------------------------------------------------------

    def estimate_sky_photutils(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky and sky noise by using photutils ...")

        bkg = Background(self.frame, (50, 50), filter_shape=(3, 3), filter_threshold=None, mask=self.mask,
                  method="sextractor", backfunc=None, interp_order=3, sigclip_sigma=3.0, sigclip_iters=10)

        # Masked background
        masked_background = np.ma.masked_array(bkg.background, mask=self.mask)
        #plotting.plot_box(masked_background, title="masked background")

        mean_sky = np.ma.mean(masked_background)
        median_sky = np.median(masked_background.compressed())

        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        self.phot_sky = Frame(bkg.background,
                               wcs=self.frame.wcs,
                               name="phot_sky",
                               description="photutils background",
                               unit=self.frame.unit,
                               zero_point=self.frame.zero_point,
                               filter=self.frame.filter,
                               sky_subtracted=False,
                               fwhm=self.frame.fwhm)

        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        self.phot_rms = Frame(bkg.background_rms,
                               wcs=self.frame.wcs,
                               name="phot_rms",
                               description="photutils rms",
                               unit=self.frame.unit,
                               zero_point=self.frame.zero_point,
                               filter=self.frame.filter,
                               sky_subtracted=False,
                               fwhm=self.frame.fwhm)

        # Create sky map of median sky level
        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        self.sky = Frame(np.full(self.frame.shape, median_sky),
                         wcs=self.frame.wcs,
                         name="sky",
                         description="estimated sky",
                         unit=self.frame.unit,
                         zero_point=self.frame.zero_point,
                         filter=self.frame.filter,
                         sky_subtracted=False,
                         fwhm=self.frame.fwhm)

    # -----------------------------------------------------------------

    def estimate_sky_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky and sky noise by using or own procedures ...")

        # Get an array of the indices of all the pixels that are not masked
        indices = np.where(not self.mask)

        print(indices)

    # -----------------------------------------------------------------

    def subtract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the sky from the frame ...")

        # Subtract the estimated sky from the image frame
        self.frame[:] = self.frame - self.sky

    # -----------------------------------------------------------------

    def set_zero_outside(self):

        """
        This function ...
        :return:
        """

        # Create a mask from the principal galaxy region
        annulus_outer_factor = 1.2
        mask = Mask.from_shape(self.principal_ellipse * annulus_outer_factor, self.frame.xsize, self.frame.ysize).inverse()

        # Set the primary frame zero outside the principal ellipse
        self.frame[mask] = 0.0

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
        return np.ma.mean(np.ma.masked_array(self.frame, mask=self.mask))

    # -----------------------------------------------------------------

    @property
    def median(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped median
        return np.median(np.ma.masked_array(self.frame, mask=self.mask).compressed())

    # -----------------------------------------------------------------

    @property
    def stddev(self):

        """
        This function ...
        :return:
        """

        # Return the standard deviation of the sigma-clipped frame
        return np.ma.masked_array(self.frame, mask=self.mask).std()

# -----------------------------------------------------------------
