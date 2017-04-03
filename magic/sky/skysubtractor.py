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
import io
import math
import imageio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.interpolate import CloughTocher2DInterpolator as intp
from scipy.interpolate import SmoothBivariateSpline
from scipy.ndimage import zoom

# Test
# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.linear_model import LinearRegression
# from sklearn.pipeline import Pipeline

# Import astronomical modules
from photutils.background import Background2D
from photutils import SigmaClip
from photutils import SExtractorBackground
from astropy.modeling import models
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter, SLSQPLSQFitter
from astropy import stats
from astropy.utils import lazyproperty
from photutils.background import BkgZoomInterpolator, BkgIDWInterpolator
from photutils.background.core import MeanBackground, MedianBackground, ModeEstimatorBackground, MMMBackground, SExtractorBackground, BiweightLocationBackground
from photutils.background.core import StdBackgroundRMS, MADStdBackgroundRMS, BiweightMidvarianceBackgroundRMS

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ..basics.mask import Mask
from ..core.detection import Detection
from ..basics.coordinate import PixelCoordinate, SkyCoordinate
from ..region.circle import PixelCircleRegion
from ..region.composite import PixelCompositeRegion
from ..region.list import PixelRegionList, SkyRegionList
from ..tools import plotting, statistics, fitting, plotting
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
#from ...core.basics.distribution import Distribution
from ..misc import chrisfuncs
from ..core.mask import Mask as newMask
from ..core.cutout import CutoutMask

# -----------------------------------------------------------------

finishing_steps = ["mean", "median", "polynomial", "interpolation"]
interpolation_methods = ["griddata", "zoom", "IDW", "spline"]

# -----------------------------------------------------------------

estimators = ["mean", "median", "mode", "MMM", "sextractor", "biweight_location"]
noise_estimators = ["stddev", "MAD", "biweight_midvariance"]

# -----------------------------------------------------------------

class SkySubtractor(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SkySubtractor, self).__init__(config, interactive)

        # -- Attributes --

        # The image frame
        self.frame = None

        # The mask of sources
        self.sources_mask = None

        # The extra mask
        self.extra_mask = None

        # The principal shape
        self.principal_shape = None

        # The region of saturated stars
        self.saturation_region = None

        # The animation
        self.animation = None

        # The sky region
        self.region = None

        # The output mask (combined input + bad mask + galaxy annulus mask + expanded saturation mask + sigma-clipping mask)
        self.mask = None

        # The estimated sky (a single floating point value or a Frame, depending on the estimation method)
        self.sky = None

        # The estimated sky noise
        self.noise = None

        # Relevant for when estimation method is 'photutils'
        self.phot_sky = None
        self.phot_rms = None

        # Relevant for when estimation method is 'pts'
        self.apertures_frame = None
        self.apertures_values_frame = None
        self.apertures_noise_frame = None
        self.apertures_mask = None

        # Aperture properties
        self.napertures = None
        self.aperture_radius = None
        self.aperture_centers = None
        self.aperture_values = None
        self.aperture_noise_values = None
        self.aperture_masks = None

        # The interpolator
        self.interpolator = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the sky region
        self.create_region()

        # 3. Create mask
        self.create_mask()

        # 4. Do an extra sigma-clipping step on the data
        if self.config.sigma_clip_mask: self.sigma_clip()

        # 5. Estimate the sky (and sky noise)
        self.estimate()

        # 6. Subtract the sky
        #self.subtract()

        # 7. Set the frame to zero outside of the principal galaxy
        if self.config.set_zero_outside: self.set_zero_outside()

        # 8. Eliminate negative values from the frame, set them to zero
        if self.config.eliminate_negatives: self.eliminate_negatives()

        # 9. Write
        if self.config.write: self.write()

        # 10. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the sky subtractor ...")

        # Set default values for all attributes
        self.frame = None
        self.sources_mask = None
        self.extra_mask = None
        self.principal_shape = None
        self.saturation_region = None
        self.animation = None
        self.mask = None
        self.sky = None
        self.noise = None
        self.phot_sky = None
        self.phot_rms = None
        self.apertures_frame = None
        self.apertures_values_frame = None
        self.apertures_noise_frame = None
        self.apertures_mask = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SkySubtractor, self).setup(**kwargs)

        # Get the frame
        if "frame" in kwargs: self.frame = kwargs.pop("frame")
        elif self.config.image is not None: self.frame = Frame.from_file(self.config.image)
        else: raise ValueError("Frame must be given as input or image path should be set in configuration")

        # Get other required
        self.sources_mask = kwargs.pop("sources_mask")

        # NOT REQUIRED ANYMORE
        if "principal_shape" in kwargs: self.principal_shape = kwargs.pop("principal_shape")

        # Get optional input
        self.extra_mask = kwargs.pop("extra_mask", None)
        self.saturation_region = kwargs.pop("saturation_region", None)
        self.animation = kwargs.pop("animation", None)

    # -----------------------------------------------------------------

    def create_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sky region ...")

        # If the sky region has to be loaded from file
        if self.config.sky_region is not None:

            sky_region = SkyRegionList.from_file(self.config.sky_region)
            self.region = sky_region.to_pixel(self.frame.wcs)

        # If no region file is given by the user, create an annulus from the principal ellipse
        elif self.principal_shape is not None:

            # Create the sky annulus
            annulus_outer_factor = self.config.mask.annulus_outer_factor
            annulus_inner_factor = self.config.mask.annulus_inner_factor
            inner_shape = self.principal_shape * annulus_inner_factor
            outer_shape = self.principal_shape * annulus_outer_factor

            # Create the annulus
            annulus = PixelCompositeRegion(outer_shape, inner_shape)

            # Create the sky region consisting of only the annulus
            self.region = PixelRegionList()
            self.region.append(annulus)

        #else: log.warning("No central region or sky regions have been defined")

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sky mask ...")

        masks = []

        if self.region is not None:

            # Create a mask from the pixels outside of the sky region
            outside_mask = self.region.to_mask(self.frame.xsize, self.frame.ysize).inverse()
            masks.append(outside_mask)

        if self.principal_shape is not None:

            # Create a mask from the principal shape
            principal_mask = self.principal_shape.to_mask(self.frame.xsize, self.frame.ysize)
            masks.append(principal_mask)

        # Sources mask
        masks.append(self.sources_mask)

        #plotting.plot_mask(outside_mask, title="outside mask")
        #plotting.plot_mask(principal_mask, title="principal mask")
        #plotting.plot_mask(self.sources_mask, title="sources mask")

        # NEW
        self.mask = newMask.union(*masks)

        # Set the mask, make a copy of the input mask initially
        #self.mask = self.sources_mask + outside_mask + principal_mask

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
        #frame.saveto("masked_sky_frame_notclipped.fits")

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

        # Use the CAAPR method
        elif self.config.estimation.method == "caapr": self.estimate_caapr()

        # Unkown sky estimation method
        else: raise ValueError("Unknown sky estimation method: '" + self.config.estimation.method + "'")

    # -----------------------------------------------------------------

    def estimate_sky_mean(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky by calculating the mean value of all non-masked pixels ...")

        # Create a frame filled with the mean value
        self.sky = self.mean

        # Determine the global noise
        self.noise = self.stddev_subtracted

        # Debugging
        log.debug("The estimated sky value is " + str(self.sky))
        log.debug("The mean sky value after subtraction is " + str(self.mean_subtracted))
        log.debug("The median sky value after subtraction is " + str(self.median_subtracted))
        log.debug("The global noise is " + str(self.noise))

    # -----------------------------------------------------------------

    def estimate_sky_median(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky by calculating the median value of all non-masked pixels ...")

        # Create a frame filled with the median value
        self.sky = self.median

        # Determine the global noise
        self.noise = self.stddev_subtracted

        # Debugging
        log.debug("The estimated sky value is " + str(self.sky))
        log.debug("The mean sky value after subtraction is " + str(self.mean_subtracted))
        log.debug("The median sky value after subtraction is " + str(self.median_subtracted))
        log.debug("The global noise is " + str(self.noise))

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

        #plotting.plot_box(data, title="estimated sky")

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

        # Set the noise
        self.noise = self.stddev_subtracted

        # Debugging
        log.debug("The mean sky value after subtraction is " + str(self.mean_subtracted))
        log.debug("The median sky value after subtraction is " + str(self.median_subtracted))
        log.debug("The global noise is " + str(self.noise))

    # -----------------------------------------------------------------

    def estimate_sky_photutils(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky and sky noise by using photutils ...")

        # Determine box size
        if self.config.estimation.photutils_fixed_width is not None: box_shape = (self.config.estimation.photutils_fixed_width, self.config.estimation.photutils_fixed_width)
        else:
            # Determine radius
            self.determine_aperture_radius()
            integer_radius = int(math.ceil(self.aperture_radius))
            box_shape = (2 * integer_radius, 2 * integer_radius)

        # Determine filter size
        filter_size = (self.config.estimation.photutils_filter_size, self.config.estimation.photutils_filter_size)

        # NEW
        sigma_clip = SigmaClip(sigma=3., iters=10)
        # bkg_estimator = MedianBackground()
        bkg_estimator = SExtractorBackground()
        bkg = Background2D(self.frame, box_shape, filter_size=filter_size, sigma_clip=sigma_clip,
                           bkg_estimator=bkg_estimator, mask=self.mask, filter_threshold=None)

        # Masked background
        masked_background = np.ma.masked_array(bkg.background, mask=self.mask)
        #plotting.plot_box(masked_background, title="masked background")

        # Masked background rms
        masked_background_rms = np.ma.masked_array(bkg.background_rms, mask=self.mask)

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

        # Use global values
        if self.config.estimation.photutils_global is not None:

            mean_sky = np.ma.mean(masked_background)
            median_sky = np.median(masked_background.compressed())

            # Median
            if self.config.estimation.photutils_global == "median":

                # Set sky level
                self.sky = median_sky

                # Set noise
                self.noise = np.ma.mean(masked_background_rms)

            # Mean
            elif self.config.estimation.photutils_global == "mean":

                # Set sky level
                self.sky = mean_sky

                # Set noise
                self.noise = np.ma.mean(masked_background_rms)

            # Invalid
            else: raise ValueError("Invalid option for 'photutils_global'")

            # Debugging
            log.debug("The estimated sky value is " + str(self.sky))
            log.debug("The mean value after subtraction is " + str(self.mean_subtracted))
            log.debug("The median value after subtraction is " + str(self.median_subtracted))
            log.debug("The standard deviation after subtraction is " + str(self.stddev_subtracted))
            log.debug("The estimated noise is " + str(self.noise))

        # Use frames
        else:

            # Set original photutils frames
            self.sky = self.phot_sky
            self.noise = self.phot_rms

            # Debugging
            log.debug("The mean value after subtraction is " + str(self.mean_subtracted))
            log.debug("The median value after subtraction is " + str(self.median_subtracted))
            log.debug("The standard deviation after subtraction is " + str(self.stddev_subtracted))
            log.debug("The estimated noise is " + str(self.noise))

    # -----------------------------------------------------------------

    def estimate_sky_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky and sky noise by randomly placing apertures ...")

        # 1. Determine the aperture radius
        self.determine_aperture_radius()

        # 2. Determine the number of apertures to use
        self.determine_number_of_apertures()

        # 3. Create the apertures
        self.create_apertures()

        # 4. Create aperture frames
        self.create_aperture_frames()

        # 5. Finish
        self.finish_sky()

    # -----------------------------------------------------------------

    def determine_aperture_radius(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the aperture radius ...")

        # Whether or not user has chosen this value
        if self.config.estimation.aperture_radius is not None: radius = self.config.estimation.aperture_radius
        else:

            # Check whether the FWHM is defined for the frame
            if self.frame.fwhm is None: raise RuntimeError("The FWHM of the frame is not defined: sky apertures cannot be generated")

            # Determine the radius for the sky apertures
            fwhm_pix = self.frame.fwhm_pix
            radius = self.config.estimation.aperture_fwhm_factor * fwhm_pix

        # Debugging
        log.debug("Using sky apertures with a radius of " + str(radius) + " pixels")

        # Set the aperture radius
        self.aperture_radius = radius

    # -----------------------------------------------------------------

    @lazyproperty
    def aperture_area(self):

        """
        This function ...
        :return:
        """

        circle_area = np.pi * self.aperture_radius ** 2
        return circle_area

    # -----------------------------------------------------------------

    def determine_number_of_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the number of apertures ...")

        npixels = np.sum(self.mask.inverse())

        # Assuming optimal hexagonal packing, get an estimate of the maximum number of circles of given radius
        # can fit in the area covered by the pixels that are not masked. This is obviously a significant overestimation
        # especially in the case where the radius becomes of the same order of magnitude as the radius of the
        # galaxy annulus (the hexagonal packing assumes a rectangular area or at least rectangular-like edges)
        # With perfect hexagonal packing, the area of the rectangle that will be covered by the circles is π/(2√3),
        # which is approximately equal to 0.907
        # See: https://www.quora.com/How-many-3-75-inch-circles-will-fit-inside-a-17-inch-square
        coverable_area = 0.907 * npixels
        optimal_number_of_apertures = coverable_area / self.aperture_area

        # Debugging
        log.debug("The upper limit to the number of apertures that fit in the part of the frame that is not masked "
                  "(assuming hexagonal packing) is " + str(optimal_number_of_apertures))

        # Determine the number of apertures that are going to be used
        napertures = int(optimal_number_of_apertures * self.config.estimation.relative_napertures_max)

        # Don't take less than 'min_napertures'
        napertures = max(napertures, self.config.min_apertures)

        # Debugging
        log.debug("A total of " + str(napertures) + " apertures are going to be used to estimate the sky ...")

        # Set the number of apertures
        self.napertures = napertures

    # -----------------------------------------------------------------

    def create_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the apertures ...")

        # Generate the apertures
        aperture_centers, aperture_values, aperture_noise_values = self.generate_apertures()

        # Remove outliers
        self.aperture_centers, self.aperture_values, self.aperture_noise_values = self.remove_aperture_outliers(aperture_centers, aperture_values, aperture_noise_values)

    # -----------------------------------------------------------------

    def generate_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the apertures ...")

        # Get arrays of the coordinates of all pixels that are not masked
        pixels_y, pixels_x = np.where(self.mask.inverse())

        # Get the number of pixels that are not masked (also the area of the frame not masked)
        npixels = pixels_x.size

        # Create a mask that tags all pixels that have been covered by one of the apertures
        apertures_mask = Mask.empty_like(self.frame)

        # Counter to keep track of the number of 'succesful' apertures that have been used
        current_napertures = 0

        # Initialize lists to contain the mean sky levels and noise levels in each of the apertures
        aperture_centers = []
        aperture_means = []
        aperture_stddevs = []
        aperture_masks = []

        # Draw random apertures
        while True:

            # Draw a random pixel index
            index = np.random.randint(npixels)

            # Get the x and y coordinate of the pixel
            x = pixels_x[index]
            y = pixels_y[index]

            # Create a coordinate for the center of the aperture
            center = PixelCoordinate(x, y)

            # Create a circular aperture
            circle = PixelCircleRegion(center, self.aperture_radius)

            # Create a Source from the frame
            source = Detection.from_shape(self.frame, circle, 1.3)

            # Get a mask of the pixels that overlap with the sky mask
            sky_mask_cutout = self.mask[source.y_slice, source.x_slice]
            overlapping = sky_mask_cutout * source.mask

            # Calculate the overlap fraction with the sky mask
            number_of_overlapping_pixels = np.sum(overlapping)
            overlap_fraction = number_of_overlapping_pixels / self.aperture_area

            # If the overlap fraction is larger than 50% for this aperture, skip it
            if overlap_fraction >= 0.5:
                log.debug("For this aperture, an overlap fraction of more than 50% was found with the sky mask, skipping ...")
                continue

            # Get a mask of the pixels that overlap with the apertures mask
            apertures_mask_cutout = apertures_mask[source.y_slice, source.x_slice]
            overlapping = apertures_mask_cutout * source.mask

            # Calculate the overlap fraction with the apertures mask
            number_of_overlapping_pixels = np.sum(overlapping)
            overlap_fraction = number_of_overlapping_pixels / self.aperture_area

            # If the overlap fraction is larger than 10% for this aperture, skip it
            if overlap_fraction >= 0.1:
                log.debug("For this aperture, an overlap fraction of more than 10% was found with other apertures, skipping ...")
                continue

            # Add the aperture area to the mask
            apertures_mask[source.y_slice, source.x_slice] += source.mask

            # Debugging
            log.debug("Placed aperture " + str(current_napertures+1) + " of " + str(self.napertures) + " ({0:.2f}%)".format((current_napertures+1)/self.napertures*100.))

            # Add a frame to the animation
            if self.animation is not None:
                plt.figure()
                plt.imshow(apertures_mask, origin="lower")
                plt.title("Aperture mask")
                buf = io.BytesIO()
                plt.savefig(buf, format='png')
                buf.seek(0)
                im = imageio.imread(buf)
                buf.close()
                self.animation.add_frame(im)

            aperture_mask = sky_mask_cutout + source.background_mask
            cutout_mask = CutoutMask(aperture_mask, source.x_min, source.x_max, source.y_min, source.y_max)

            # Calculate the mean sky value in this aperture
            masked_array_cutout = np.ma.MaskedArray(source.cutout, mask=aperture_mask)

            # plotting.plot_box(masked_array_cutout)

            aperture_mean = np.ma.mean(masked_array_cutout)
            #aperture_median = np.ma.median(masked_array_cutout)
            # aperture_median2 = np.median(masked_array_cutout.compressed()) # same result, but unnecessary compressed step
            aperture_stddev = np.std(masked_array_cutout)

            # print("aperture mean:", aperture_mean)
            # print("aperture median:", aperture_median, aperture_median2)
            # print("aperture stddev:", aperture_std)

            # Add properties of the aperture
            aperture_centers.append(center)
            aperture_means.append(aperture_mean)
            aperture_stddevs.append(aperture_stddev)
            aperture_masks.append(cutout_mask)

            # Another succesful aperture
            current_napertures += 1

            # Stop when we have reached the desired number of apertures
            if current_napertures == self.napertures: break

        # Create Numpy arrays from the aperture means and standard deviations
        aperture_means = np.array(aperture_means)
        aperture_stddevs = np.array(aperture_stddevs)

        # Return the aperture properties
        return aperture_centers, aperture_means, aperture_stddevs

    # -----------------------------------------------------------------

    def remove_aperture_outliers(self, aperture_centers, aperture_means, aperture_stddevs):

        """
        This function ...
        :param aperture_centers:
        :param aperture_means:
        :param aperture_stddevs:
        :return:
        """

        # Inform the user
        log.info("Removing aperture outliers ...")

        #means_distribution = Distribution.from_values(aperture_means, bins=50)
        #stddevs_distribution = Distribution.from_values(aperture_stddevs, bins=50)

        #means_distribution.plot("Aperture means before sigma-clipping")
        #stddevs_distribution.plot("Aperture stddevs before sigma-clipping")

        clip_mask = stats.sigma_clip(aperture_stddevs, sigma=3.0, iters=None, copy=False).mask

        clipped_aperture_centers = []
        for i in range(len(clip_mask)):
            if clip_mask[i]: continue
            else: clipped_aperture_centers.append(aperture_centers[i])

        aperture_centers = clipped_aperture_centers
        aperture_means = np.ma.MaskedArray(aperture_means, clip_mask).compressed()
        aperture_stddevs = np.ma.MaskedArray(aperture_stddevs, clip_mask).compressed()

        #means_distribution = Distribution.from_values(aperture_means, bins=50)
        #stddevs_distribution = Distribution.from_values(aperture_stddevs, bins=50)

        #means_distribution.plot("Aperture means after sigma-clipping")
        #stddevs_distribution.plot("Aperture stddevs after sigma-clipping")

        # Return the sigma-clipped aperture properties
        return aperture_centers, aperture_means, aperture_stddevs

    # -----------------------------------------------------------------

    def create_aperture_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating aperture frames ...")

        # Initialize frames and mask
        self.apertures_frame = Frame.nans_like(self.frame)
        self.apertures_mean_frame = Frame.nans_like(self.frame)
        self.apertures_noise_frame = Frame.nans_like(self.frame)
        self.apertures_mask = newMask.empty_like(self.frame)

        # Loop over the apertures
        for i in range(len(self.aperture_centers)):

            # Create mask
            center = self.aperture_centers[i]
            circle = PixelCircleRegion(center, self.aperture_radius)
            mask = Mask.from_shape(circle, self.frame.xsize, self.frame.ysize)

            # Set
            self.apertures_frame[mask] = self.frame[mask]
            self.apertures_mean_frame[mask] = self.aperture_values[i]
            self.apertures_noise_frame[mask] = self.aperture_noise_values[i]
            self.apertures_mask[mask][self.aperture_masks[i]] = True

    # -----------------------------------------------------------------

    def finish_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finishing the sky map ...")

        # Finishing step
        if self.config.estimation.finishing_step == "mean": self.set_sky_mean()
        elif self.config.estimation.finishing_step == "median": self.set_sky_median()
        elif self.config.estimation.finishing_step == "polynomial": self.fit_polynomial_to_apertures()
        elif self.config.estimation.finishing_step == "interpolation": self.interpolate_apertures()
        else: raise ValueError("Invalid finishing step")

    # -----------------------------------------------------------------

    def set_sky_mean(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the sky to the mean of the aperture fluxes ...")

        # Determine the mean sky level
        self.sky = np.mean(self.aperture_values)

        # Determine the global noise level
        self.noise = self.calculate_global_noise_level_constant_sky()

        # Debugging
        log.debug("The estimated sky level is " + str(self.sky))
        log.debug("The estimated sky noise level is " + str(self.noise))

    # -----------------------------------------------------------------

    def set_sky_median(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the sky to the median of the aperture fluxes ...")

        # Determine the median sky level
        self.sky = np.median(self.aperture_values)

        # Calculate the global noise level
        self.noise = self.calculate_global_noise_level_constant_sky()

        # Debugging
        log.debug("The estimated sky level is " + str(self.sky))
        log.debug("The estimated sky noise level is " + str(self.noise))

    # -----------------------------------------------------------------

    def calculate_global_noise_level_constant_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the global noise value for methods that calculate a single, constant sky value for the entire image ...")

        # Calculate large scale variation and pixel to pixel noise
        large_scale_variations_error, pixel_to_pixel_noise = self.calculate_large_scale_variation_and_pixel_to_pixel_noise()

        # Determine the noise by quadratically adding the large scale variation and the mean pixel-by-pixel noise
        value = np.sqrt(large_scale_variations_error ** 2 + pixel_to_pixel_noise ** 2)

        # Return the noise value
        return value

    # -----------------------------------------------------------------

    def calculate_large_scale_variation_and_pixel_to_pixel_noise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the large scale variation and the pixel to pixel noise of the apertures ...")

        # Calculate the large-scale variation level
        large_scale_variations_error = self.aperture_values.std()

        # Calculate the mean pixel-by-pixel noise over all apertures
        pixel_to_pixel_noise = np.mean(self.aperture_noise_values)

        # Return the result
        return large_scale_variations_error, pixel_to_pixel_noise

    # -----------------------------------------------------------------

    def calculate_global_noise_level_variate_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the global noise value for methods that calculate a variable sky map ...")

        # Calculate large scale variation and pixel to pixel noise
        large_scale_variations_error, pixel_to_pixel_noise = self.calculate_large_scale_variation_and_pixel_to_pixel_noise()

        # Return the average pixel-to-pixel noise
        return pixel_to_pixel_noise

    # -----------------------------------------------------------------

    def fit_polynomial_to_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting a polynomial to the values in each aperture ...")

        # Get list of x and y coordinates
        x_values = [center.x for center in self.aperture_centers]
        y_values = [center.y for center in self.aperture_centers]

        # -- Fit polynomial --

        # Fit polynomial to aperture means
        poly_init = models.Polynomial2D(degree=self.config.estimation.polynomial_degree)
        #fit_model = LevMarLSQFitter()
        #fit_model = NonLinearLSQFitter()
        fit_model = SLSQPLSQFitter()
        polynomial = fit_model(poly_init, x_values, y_values, self.aperture_values)

        # Create x and y meshgrid for evaluating
        y_grid, x_grid = np.mgrid[:self.frame.ysize, :self.frame.xsize]

        # Evaluate the model
        data = polynomial(x_grid, y_grid)

        # Set the sky
        self.sky = Frame(data)

        # plotting.plot_box(data)

        # -- Fit spline --

        #f = interpolate.interp2d(x_values, y_values, aperture_means, kind='cubic')


        #x_grid = np.array(range(self.frame.xsize))
        #y_grid = np.array(range(self.frame.ysize))

        #data = f(x_grid, y_grid)

        # Set new sky frame
        #self.sky = Frame(data)

        ## NEW: NOISE

        # Calculate the global noise level
        self.noise = self.calculate_global_noise_level_variate_sky()

    # -----------------------------------------------------------------

    def interpolate_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating between the values of each aperture to fill the sky frame ...")

        # Choose method
        if self.config.estimation.interpolation_method == "griddata": self.interpolate_griddata()
        elif self.config.estimation.interpolation_method == "zoom": self.interpolate_zoom()
        elif self.config.estimation.interpolation_method == "IDW": self.interpolate_idw()
        elif self.config.estimation.interpolation_method == "spline": self.interpolate_spline()
        else: raise ValueError("Invalid interpolation method: " + self.config.estimation.interpolation_method)

    # -----------------------------------------------------------------

    def interpolate_griddata(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the aperture fluxes using the 'griddata' method")

        # Get lists of x and y values
        x_values = np.array([center.x for center in self.aperture_centers])
        y_values = np.array([center.y for center in self.aperture_centers])

        x_ticks = np.arange(0, self.frame.xsize, 1)
        y_ticks = np.arange(0, self.frame.ysize, 1)
        z_grid = mlab.griddata(x_values, y_values, aperture_means, x_ticks, y_ticks)

        # Set the sky frame
        self.sky = Frame(z_grid)

        # Calculate the global noise level
        self.noise = self.calculate_global_noise_level_variate_sky()

    # -----------------------------------------------------------------

    def interpolate_zoom(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the aperture fluxes using the 'zoom' method ...")

        # UNFORTUNATELY, something like this is not possible

        # Create the interpolator
        #self.interpolator = BkgZoomInterpolator()

        # Interpolate
        #data = self.interpolator(self.background_mesh, self)
        #self.sky = Frame(data)

        # IMPLEMENTATION OF BKGZOOMINTERPOLATOR:

        mesh = np.asanyarray(mesh)
        if np.ptp(mesh) == 0:
            return np.zeros_like(bkg2d_obj.data) + np.min(mesh)

        #if bkg2d_obj.edge_method == 'pad':
        if False:

            # The mesh is first resized to the larger padded-data size
            # (i.e. zoom_factor should be an integer) and then cropped
            # back to the final data size.

            zoom_factor = (int(bkg2d_obj.nyboxes * bkg2d_obj.box_size[0] /
                               mesh.shape[0]),
                           int(bkg2d_obj.nxboxes * bkg2d_obj.box_size[1] /
                               mesh.shape[1]))
            result = zoom(mesh, zoom_factor, order=self.order, mode=self.mode,
                          cval=self.cval)

            return result[0:bkg2d_obj.data.shape[0], 0:bkg2d_obj.data.shape[1]]

        else:

            # The mesh is resized directly to the final data size.
            zoom_factor = (float(bkg2d_obj.data.shape[0] / mesh.shape[0]),
                           float(bkg2d_obj.data.shape[1] / mesh.shape[1]))

            return zoom(mesh, zoom_factor, order=self.order, mode=self.mode,
                        cval=self.cval)

    # -----------------------------------------------------------------

    def interpolate_idw(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the aperture fluxes using the 'IDW' method ...")

        # Create the interpolator
        self.interpolator = BkgIDWInterpolator()

        # Interpolate
        data = self.interpolator(self.background_mesh, self)
        self.sky = Frame(data)

        # IMPLEMENTATION OF BKGIDWINTERPOLATOR:

        mesh = np.asanyarray(mesh)
        if np.ptp(mesh) == 0:
            return np.zeros_like(bkg2d_obj.data) + np.min(mesh)

        mesh1d = mesh[bkg2d_obj.mesh_yidx, bkg2d_obj.mesh_xidx]
        f = ShepardIDWInterpolator(bkg2d_obj.yx, mesh1d,
                                   leafsize=self.leafsize)
        data = f(bkg2d_obj.data_coords, n_neighbors=self.n_neighbors,
                 power=self.power, reg=self.reg)

        data = data.reshape(bkg2d_obj.data.shape)



    # -----------------------------------------------------------------

    def interpolate_spline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the aperture fluxes using the 'spline' method ...")

        # Get lists of the x and y values
        x_values = np.array([center.x for center in self.aperture_centers])
        y_values = np.array([center.y for center in self.aperture_centers])

        #X, Y = np.meshgrid(x_values, y_values)

        X = x_values
        Y = y_values
        Z = self.aperture_means

        #print(X, Y, Z)
        #print(len(X), len(Y), len(Z))

        #C = intp((X, Y), Z)
        x_space = np.linspace(0, self.frame.xsize, 1)
        y_space = np.linspace(0, self.frame.ysize, 1)
        xi, yi = np.meshgrid(x_space, y_space)
        #zi = C(xi, yi)

        #self.sky = Frame(zi)

        from scipy.interpolate import LSQBivariateSpline

        spline = SmoothBivariateSpline(X, Y, Z, kx=1, ky=1)
        #spline = LSQBivariateSpline(X, Y, Z, X, Y)
        #zi = spline(xi, yi)
        #self.sky = Frame(zi)

        from scipy.interpolate import griddata


        #x_space = np.linspace(0.3*self.frame.xsize, 0.7*self.frame.xsize)
        #y_space = np.linspace(0.3*self.frame.ysize, 0.7*self.frame.ysize)

        x_space = np.array(range(int(0.3*self.frame.xsize), int(0.7*self.frame.xsize)))
        y_space = np.array(range(int(0.3*self.frame.ysize), int(0.7*self.frame.ysize)))

        znew = griddata((X, Y), Z, (x_space[None,:], y_space[:,None]), method='cubic')

        plt.figure()
        levels = np.linspace(min(Z), max(Z), 15)
        plt.ylabel('Y', size=15)
        plt.xlabel('X', size=15)
        cmap = plt.cm.jet_r
        cs = plt.contourf(x_space, y_space, znew, levels=levels, cmap=cmap)
        cbar = plt.colorbar(cs)
        cbar.set_label('Z', rotation=90, fontsize=15)  # gas fraction
        plt.show()

        self.sky = Frame.zeros_like(self.frame)
        self.sky[int(0.3*self.frame.ysize):int(0.3*self.frame.ysize)+len(y_space), int(0.3*self.frame.xsize):int(0.3*self.frame.xsize)+len(x_space)] = znew

        #self.sky = Frame(znew)

    # -----------------------------------------------------------------

    def plot_interpolated(self, aperture_centers, aperture_means):

        """
        This function ...
        :param aperture_centers:
        :param aperture_means:
        :return:
        """

        x_values = np.array([center.x for center in aperture_centers])
        y_values = np.array([center.y for center in aperture_centers])

        x_ticks = np.arange(0, self.frame.xsize, 1)
        y_ticks = np.arange(0, self.frame.ysize, 1)
        z_grid = mlab.griddata(x_values, y_values, aperture_means, x_ticks, y_ticks)

        self.sky = Frame(z_grid)

        from matplotlib.backends import backend_agg as agg
        from matplotlib import cm

        # plot
        #fig = Figure()  # create the figure
        fig = plt.figure()
        agg.FigureCanvasAgg(fig)  # attach the rasterizer
        ax = fig.add_subplot(1, 1, 1)  # make axes to plot on
        ax.set_title("Interpolated Contour Plot of Experimental Data")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        cmap = cm.get_cmap("hot")  # get the "hot" color map
        contourset = ax.contourf(x_ticks, y_ticks, z_grid, 10, cmap=cmap)

        cbar = fig.colorbar(contourset)
        cbar.set_ticks([0, 100])
        fig.axes[-1].set_ylabel("Z")  # last axes instance is the colorbar

        plt.show()

    # -----------------------------------------------------------------

    def estimate_caapr(self):

        """
        This function ...
        :return:
        """

        # Run pod through function that removes large-scale sky using a 2-dimensional polynomial filter
        #pod = CAAPR.CAAPR_Pipeline.PolySub(pod, 2.0 * pod['semimaj_initial_pix'], pod['opt_axial_ratio'],
        #                                   pod['opt_angle'],
        #                                   instant_quit=max([not kwargs_dict['polysub'], pod['band_exclude']]))

        if pod['verbose']: print('[' + pod['id'] + '] Determining if (and how) background is significantly variable.')

        # Define Keflavich function to downsample an array
        def Downsample(myarr, factor, estimator=np.nanmean):
            ys, xs = myarr.shape
            crarr = myarr[:ys - (ys % int(factor)), :xs - (xs % int(factor))]
            dsarr = estimator(np.concatenate([[crarr[i::factor, j::factor]
                                               for i in range(factor)]
                                              for j in range(factor)]), axis=0)
            return dsarr

        # If polynomial background subraction not wanted, immediately return everything unchanged
        if instant_quit:
            pod['sky_poly'] = False
            return pod

        # If image has pixels smaller than some limit, downsample image to improve processing time
        pix_size = pod['pix_arcsec']
        pix_size_limit = 2.0
        if pix_size < pix_size_limit:
            downsample_factor = int(np.ceil(pix_size_limit / pix_size))
        else:
            downsample_factor = 1
        image_ds = Downsample(pod['cutout'], downsample_factor)

        # Downsample related values accordingly
        mask_semimaj_pix = mask_semimaj_pix / downsample_factor
        centre_i = int(round(float((0.5 * pod['centre_i']) - 1.0)))
        centre_j = int(round(float((0.5 * pod['centre_j']) - 1.0)))

        # Find cutoff for excluding bright pixels by sigma-clipping map
        clip_value = chrisfuncs.SigmaClip(image_ds, tolerance=0.01, sigma_thresh=3.0, median=True)
        noise_value = clip_value[0]
        field_value = clip_value[1]
        cutoff = field_value + (cutoff_sigma * noise_value)

        # Mask all image pixels in masking region around source
        image_masked = image_ds.copy()
        ellipse_mask = chrisfuncs.EllipseMask(image_ds, mask_semimaj_pix, mask_axial_ratio, mask_angle, centre_i, centre_j)
        image_masked[np.where(ellipse_mask == 1)] = np.nan

        # Mask all image pixels identified as being high SNR
        image_masked[np.where(image_masked > cutoff)] = np.nan

        # Use astropy to fit 2-dimensional polynomial to the image
        image_masked[np.where(np.isnan(image_masked) == True)] = field_value
        poly_model = astropy.modeling.models.Polynomial2D(degree=poly_order)
        i_coords, j_coords = np.mgrid[:image_masked.shape[0], :image_masked.shape[1]]
        fitter = astropy.modeling.fitting.LevMarLSQFitter()
        i_coords = i_coords.flatten()
        j_coords = j_coords.flatten()
        image_flattened = image_masked.flatten()
        good = np.where(np.isnan(image_flattened) == False)
        i_coords = i_coords[good]
        j_coords = j_coords[good]
        image_flattened = image_flattened[good]
        fit = fitter(poly_model, i_coords, j_coords, image_flattened)

        # Create final polynomial filter (undoing downsampling using lorenzoriano GitHub script)
        i_coords, j_coords = np.mgrid[:image_ds.shape[0], :image_ds.shape[1]]
        poly_fit = fit(i_coords, j_coords)
        poly_full = scipy.ndimage.interpolation.zoom(poly_fit,
                                                     [float(pod['cutout'].shape[0]) / float(poly_fit.shape[0]),
                                                      float(pod['cutout'].shape[1]) / float(poly_fit.shape[1])],
                                                     mode='nearest')  # poly_full = congrid.congrid(poly_fit, (pod['cutout'].shape[0], pod['cutout'].shape[1]), minusone=True)

        # Establish background variation before application of filter
        sigma_thresh = 3.0
        clip_in = chrisfuncs.SigmaClip(pod['cutout'], tolerance=0.005, median=True, sigma_thresh=sigma_thresh)
        bg_in = pod['cutout'][np.where(pod['cutout'] < clip_in[1])]
        spread_in = np.mean(np.abs(bg_in - clip_in[1]))

        # How much reduction in background variation there was due to application of the filter
        image_sub = pod['cutout'] - poly_full
        clip_sub = chrisfuncs.SigmaClip(image_sub, tolerance=0.005, median=True, sigma_thresh=sigma_thresh)
        bg_sub = image_sub[np.where(image_sub < clip_sub[1])]
        spread_sub = np.mean(np.abs(bg_sub - clip_sub[1]))
        spread_diff = spread_in / spread_sub

        # If the filter made significant difference, apply to image and return it; otherwise, just return the unaltered map
        if spread_diff > 1.1:
            if pod['verbose']: print('[' + pod['id'] + '] Background is significantly variable; removing polynomial background fit.')
            pod['cutout_nopoly'] = pod['cutout'].copy()
            pod['cutout'] = image_sub
            pod['sky_poly'] = poly_model
        else:
            if pod['verbose']: print('[' + pod['id'] + '] Background is not significantly variable; leaving image unaltered.')
            pod['sky_poly'] = False

        #return pod

    # -----------------------------------------------------------------

    #def subtract(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Inform the user
        #log.info("Subtracting the sky from the frame ...")

        # Subtract the estimated sky from the image frame
        #self.frame -= self.sky

    # -----------------------------------------------------------------

    def set_zero_outside(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the frame zero outside of the principal galaxy ...")

        # Create a mask from the principal galaxy region
        factor = self.config.zero_outside.factor
        mask = Mask.from_shape(self.principal_shape * factor, self.frame.xsize, self.frame.ysize).inverse()

        # Set the primary frame zero outside the principal ellipse
        self.frame[mask] = 0.0

    # -----------------------------------------------------------------

    def eliminate_negatives(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting pixels with negative value to zero ...")

        # Set all negative pixels to zero
        self.frame[self.frame <= 0.] = 0.0

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write histrogram
        #self.write_histogram()

        # Write apertures frame
        if self.config.estimation.method == "pts": self.write_apertures()

        # Write aperture means
        if self.config.estimation.method == "pts": self.write_apertures_value()

        # Write aperture noise
        if self.config.estimation.method == "pts": self.write_apertures_noise()

        # Write region
        self.write_region()

        # Write mask
        self.write_mask()

        # Write sky
        self.write_sky()

        # Write noise
        self.write_noise()

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

    def write_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the aperture frame ...")

        # Determine the path
        path = self.output_path_file("apertures.fits")

        # Save
        self.apertures_frame.saveto(path)

    # -----------------------------------------------------------------

    def write_apertures_value(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the aperture values frame ...")

        # Determine the path
        path = self.output_path_file("aperture_values.fits")

        # Save
        self.apertures_mean_frame.saveto(path)

    # -----------------------------------------------------------------

    def write_apertures_noise(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Writing the aperture noise frame ...")

        # Determine the path
        path = self.output_path_file("aperture_noise.fits")

        # Save
        self.apertures_noise_frame.saveto(path)

    # -----------------------------------------------------------------

    def write_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the region ...")

        # Determine the path
        path = self.output_path_file("sky.reg")

        # Save
        self.region.saveto(path)

    # -----------------------------------------------------------------

    def write_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mask ...")

        # Determine the path
        path = self.output_path_file("mask.fits")

        # Save
        self.mask.saveto(path)

    # -----------------------------------------------------------------

    def write_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sky ...")

        # Determine the path
        path = self.output_path_file("sky.fits")

        # Save
        self.sky_frame.saveto(path)

    # -----------------------------------------------------------------

    def write_noise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the noise ...")

        # Determine the path
        path = self.output_path_file("noise.fits")

        # Save
        self.noise_frame.saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def sky_frame(self):

        """
        This function ...
        :return:
        """

        if isinstance(self.sky, Frame): return self.sky
        else: return Frame(np.full(self.frame.shape, self.sky))

    # -----------------------------------------------------------------

    @property
    def noise_frame(self):

        """
        This function ...
        :return:
        """

        if isinstance(self.noise, Frame): return self.noise
        else: return Frame(np.full(self.frame.shape, self.noise))

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

    @lazyproperty
    def subtracted(self):

        """
        This function ...
        :return:
        """

        return self.frame - self.sky

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_subtracted(self):

        """
        This function ...
        :return:
        """

        return np.ma.mean(np.ma.masked_array(self.subtracted, mask=self.mask))

    # -----------------------------------------------------------------

    @lazyproperty
    def median_subtracted(self):

        """
        This function ...
        :return:
        """

        return np.median(np.ma.masked_array(self.subtracted, mask=self.mask).compressed())

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_subtracted(self):

        """
        This function ...
        :return:
        """

        return np.ma.masked_array(self.subtracted, mask=self.mask).std()

# -----------------------------------------------------------------
