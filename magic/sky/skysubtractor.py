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
import imageio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from scipy.interpolate import CloughTocher2DInterpolator as intp
from scipy.interpolate import SmoothBivariateSpline

# Test
# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.linear_model import LinearRegression
# from sklearn.pipeline import Pipeline

# Import astronomical modules
from photutils.background import Background
from astropy.modeling import models
from astropy.modeling.fitting import LevMarLSQFitter
from astropy import stats

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ..basics.mask import Mask
from ..core.source import Source
from ..basics.geometry import Coordinate, Circle, Composite
from ..basics.region import Region
from ..basics.skyregion import SkyRegion
from ..tools import plotting, statistics, fitting, plotting
from ...core.basics.configurable import OldConfigurable
from ...core.tools.logging import log
from ...core.basics.distribution import Distribution

# -----------------------------------------------------------------

class SkySubtractor(OldConfigurable):

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
        self.apertures_mean_frame = None
        self.apertures_noise_frame = None

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

    def run(self, frame, principal_shape, sources_mask, extra_mask=None, saturation_region=None, animation=None):

        """
        This function ...
        :param frame:
        :param principal_shape:
        :param sources_mask:
        :param extra_mask:
        :param saturation_region:
        :param animation:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, principal_shape, sources_mask, extra_mask, saturation_region, animation)

        # 2. Create the sky region
        self.create_region()

        # 3. Create mask
        self.create_mask()

        # 4. Do an extra sigma-clipping step on the data
        if self.config.sigma_clip_mask: self.sigma_clip()

        # 5. Estimate the sky (and sky noise)
        self.estimate()

        # 6. Subtract the sky
        self.subtract()

        # 7. Set the frame to zero outside of the principal galaxy
        if self.config.set_zero_outside: self.set_zero_outside()

        # 8. Eliminate negative values from the frame, set them to zero
        if self.config.eliminate_negatives: self.eliminate_negatives()

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
        self.apertures_mean_frame = None
        self.apertures_noise_frame = None

    # -----------------------------------------------------------------

    def setup(self, frame, principal_shape, sources_mask, extra_mask=None, saturation_region=None, animation=None):

        """
        This function ...
        :param frame:
        :param principal_shape:
        :param sources_mask:
        :param extra_mask:
        :param saturation_region:
        :param animation:
        :return:
        """

        # Call the setup function of the base class
        super(SkySubtractor, self).setup()

        # Make a local reference to the image frame
        self.frame = frame

        # Make a reference to the principal shape
        self.principal_shape = principal_shape

        # Set the masks
        self.sources_mask = sources_mask
        self.extra_mask = extra_mask

        # Set the saturation_region
        self.saturation_region = saturation_region

        # Make a reference to the animation
        self.animation = animation

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
            sky_region = SkyRegion.from_file(self.config.sky_region)
            self.region = sky_region.to_pixel(self.frame.wcs)

        # If no region file is given by the user, create an annulus from the principal ellipse
        else:

            # Create the sky annulus
            annulus_outer_factor = self.config.mask.annulus_outer_factor
            annulus_inner_factor = self.config.mask.annulus_inner_factor
            inner_shape = self.principal_shape * annulus_inner_factor
            outer_shape = self.principal_shape * annulus_outer_factor

            # Create the annulus
            annulus = Composite(outer_shape, inner_shape)

            # Create the sky region consisting of only the annulus
            self.region = Region()
            self.region.append(annulus)

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sky mask ...")

        # Create a mask from the pixels outside of the sky region
        outside_mask = self.region.to_mask(self.frame.xsize, self.frame.ysize).inverse()

        # Create a mask from the principal shape
        principal_mask = self.principal_shape.to_mask(self.frame.xsize, self.frame.ysize)

        #plotting.plot_mask(outside_mask, title="outside mask")
        #plotting.plot_mask(principal_mask, title="principal mask")
        #plotting.plot_mask(self.sources_mask, title="sources mask")

        # Set the mask, make a copy of the input mask initially
        self.mask = self.sources_mask + outside_mask + principal_mask

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
        self.sky = self.mean

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

        # Set sky level
        self.sky = median_sky

    # -----------------------------------------------------------------

    def estimate_sky_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky and sky noise by using or own procedures ...")

        # Check whether the FWHM is defined for the frame
        if self.frame.fwhm is None: raise RuntimeError("The FWHM of the frame is not defined: sky apertures cannot be generated")

        # Determine the aperture radius
        aperture_radius = self.determine_aperture_radius()

        # Determine the number of apertures to use
        napertures = self.determine_number_of_apertures(aperture_radius)

        # Generate the apertures
        aperture_centers, aperture_means, aperture_stddevs = self.generate_apertures(aperture_radius, napertures)

        # Remove outliers
        aperture_centers, aperture_means, aperture_stddevs = self.remove_aperture_outliers(aperture_centers, aperture_means, aperture_stddevs)

        # Calculate the large-scale variation level
        large_scale_variations_error = aperture_means.std()

        # Calculate the mean pixel-by-pixel noise over all apertures
        pixel_to_pixel_noise = np.mean(aperture_stddevs)

        # Determine the median sky level
        self.sky = np.median(aperture_means)

        # Determine the noise by quadratically adding the large scale variation and the mean pixel-by-pixel noise
        self.noise = np.sqrt(large_scale_variations_error**2 + pixel_to_pixel_noise**2)

        # Debugging
        log.debug("The estimated sky level is " + str(self.sky))
        log.debug("The estimated sky noise level is " + str(self.noise))

        # Create aperture frames
        self.create_aperture_frames(aperture_centers, aperture_means, aperture_stddevs, aperture_radius)

        # Finishing step
        if self.config.estimation.finishing_step is None: pass
        elif self.config.estimation.finishing_step == "polynomial": self.fit_polynomial_to_apertures(aperture_centers, aperture_means)
        elif self.config.estimation.finishing_step == "interpolation": self.interpolate_apertures(aperture_centers, aperture_means)
        else: raise ValueError("Invalid finishing step")

        #self.plot_interpolated(aperture_centers, aperture_means)
        #self.try_to_interpolate_smart(aperture_centers, aperture_means)

    # -----------------------------------------------------------------

    def determine_aperture_radius(self):

        """
        This function ...
        :return:
        """

        # Determine the radius for the sky apertures
        fwhm_pix = self.frame.fwhm_pix
        radius = 4.0 * fwhm_pix

        # Debugging
        log.debug("Using sky apertures with a radius of " + str(radius) + " pixels")

        # Return the aperture radius
        return radius

    # -----------------------------------------------------------------

    def determine_number_of_apertures(self, radius):

        """
        This function ...
        :param radius:
        :return:
        """

        npixels = np.sum(self.mask.inverse())

        # Assuming optimal hexagonal packing, get an estimate of the maximum number of circles of given radius
        # can fit in the area covered by the pixels that are not masked. This is obviously a significant overestimation
        # especially in the case where the radius becomes of the same order of magnitude as the radius of the
        # galaxy annulus (the hexagonal packing assumes a rectangular area or at least rectangular-like edges)
        # With perfect hexagonal packing, the area of the rectangle that will be covered by the circles is π/(2√3),
        # which is approximately equal to 0.907
        # See: https://www.quora.com/How-many-3-75-inch-circles-will-fit-inside-a-17-inch-square
        coverable_area = 0.907 * npixels
        circle_area = np.pi * radius ** 2
        optimal_number_of_apertures = coverable_area / circle_area

        # Debugging
        log.debug("The upper limit to the number of apertures that fit in the part of the frame that is not masked "
                  "(assuming hexagonal packing) is " + str(optimal_number_of_apertures))

        # Determine the number of apertures that are going to be used, take a third of the upper limit
        napertures = int(optimal_number_of_apertures / 3.)

        # Debugging
        log.debug("A total of " + str(napertures) + " apertures are going to be used to estimate the sky ...")

        # Return the number of apertures
        return napertures

    # -----------------------------------------------------------------

    def generate_apertures(self, radius, napertures):

        """
        This function ...
        :param radius:
        :param napertures:
        :return:
        """

        circle_area = np.pi * radius ** 2

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

        # Draw 100 random coordinates
        while True:

            # Draw a random pixel index
            index = np.random.randint(npixels)

            # Get the x and y coordinate of the pixel
            x = pixels_x[index]
            y = pixels_y[index]

            # Create a coordinate for the center of the aperture
            center = Coordinate(x, y)

            # Create a circular aperture
            circle = Circle(center, radius)

            # Create a Source from the frame
            source = Source.from_shape(self.frame, circle, 1.3)

            # Get a mask of the pixels that overlap with the sky mask
            sky_mask_cutout = self.mask[source.y_slice, source.x_slice]
            overlapping = sky_mask_cutout * source.mask

            # Calculate the overlap fraction with the sky mask
            number_of_overlapping_pixels = np.sum(overlapping)
            overlap_fraction = number_of_overlapping_pixels / circle_area

            # If the overlap fraction is larger than 50% for this aperture, skip it
            if overlap_fraction >= 0.5:
                log.debug(
                    "For this aperture, an overlap fraction of more than 50% was found with the sky mask, skipping ...")
                continue

            # Get a mask of the pixels that overlap with the apertures mask
            apertures_mask_cutout = apertures_mask[source.y_slice, source.x_slice]
            overlapping = apertures_mask_cutout * source.mask

            # Calculate the overlap fraction with the apertures mask
            number_of_overlapping_pixels = np.sum(overlapping)
            overlap_fraction = number_of_overlapping_pixels / circle_area

            # If the overlap fraction is larger than 10% for this aperture, skip it
            if overlap_fraction >= 0.1:
                log.debug(
                    "For this aperture, an overlap fraction of more than 10% was found with other apertures, skipping ...")

            # Add the aperture area to the mask
            apertures_mask[source.y_slice, source.x_slice] += source.mask

            # Debugging
            log.debug("Placed aperture " + str(current_napertures+1) + " of " + str(napertures) + " ({0:.2f}%)".format((current_napertures+1)/napertures*100.))

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

            # Calculate the mean sky value in this aperture
            masked_array_cutout = np.ma.MaskedArray(source.cutout, mask=sky_mask_cutout + source.background_mask)

            # plotting.plot_box(masked_array_cutout)

            aperture_mean = np.ma.mean(masked_array_cutout)
            #aperture_median = np.ma.median(masked_array_cutout)
            # aperture_median2 = np.median(masked_array_cutout.compressed()) # same result, but unnecessary compressed step
            aperture_stddev = np.std(masked_array_cutout)

            # print("aperture mean:", aperture_mean)
            # print("aperture median:", aperture_median, aperture_median2)
            # print("aperture stddev:", aperture_std)

            aperture_centers.append(center)
            aperture_means.append(aperture_mean)
            aperture_stddevs.append(aperture_stddev)

            # Another succesful aperture
            current_napertures += 1

            # Stop when we have reached the desired number of apertures
            if current_napertures == napertures: break

        # Create Numpy arrays from the aperture means and standard deviations
        aperture_means = np.array(aperture_means)
        aperture_stddevs = np.array(aperture_stddevs)

        # Return the aperture properties
        return aperture_centers, aperture_means, aperture_stddevs

    # -----------------------------------------------------------------

    def remove_aperture_outliers(self, aperture_centers, aperture_means, aperture_stddevs):

        """
        This function ...
        :return:
        """

        means_distribution = Distribution.from_values(aperture_means, bins=50)
        stddevs_distribution = Distribution.from_values(aperture_stddevs, bins=50)

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

        means_distribution = Distribution.from_values(aperture_means, bins=50)
        stddevs_distribution = Distribution.from_values(aperture_stddevs, bins=50)

        #means_distribution.plot("Aperture means after sigma-clipping")
        #stddevs_distribution.plot("Aperture stddevs after sigma-clipping")

        # Return the sigma-clipped aperture properties
        return aperture_centers, aperture_means, aperture_stddevs

    # -----------------------------------------------------------------

    def create_aperture_frames(self, aperture_centers, aperture_means, aperture_stddevs, aperture_radius):

        """
        This function ...
        :param aperture_centers:
        :param aperture_means:
        :param aperture_stddevs:
        :param aperture_radius:
        :return:
        """

        self.apertures_frame = Frame.nans_like(self.frame)
        self.apertures_mean_frame = Frame.nans_like(self.frame)
        self.apertures_noise_frame = Frame.nans_like(self.frame)

        for i in range(len(aperture_centers)):

            center = aperture_centers[i]

            circle = Circle(center, aperture_radius)

            mask = Mask.from_shape(circle, self.frame.xsize, self.frame.ysize)

            self.apertures_frame[mask] = self.frame[mask]
            self.apertures_mean_frame[mask] = aperture_means[i]
            self.apertures_noise_frame[mask] = aperture_stddevs[i]

    # -----------------------------------------------------------------

    def fit_polynomial_to_apertures(self, aperture_centers, aperture_means):

        """
        This function ...
        :return:
        """

        x_values = [center.x for center in aperture_centers]
        y_values = [center.y for center in aperture_centers]

        # -- Fit polynomial --

        # Fit polynomial to aperture means
        degree = 4
        poly_init = models.Polynomial2D(degree=degree)
        fit_model = LevMarLSQFitter()
        polynomial = fit_model(poly_init, x_values, y_values, aperture_means)

        # Create x and y meshgrid for evaluating
        y_grid, x_grid = np.mgrid[:self.frame.ysize, :self.frame.xsize]

        # Evaluate the model
        data = polynomial(x_grid, y_grid)

        self.sky = Frame(data)

        # plotting.plot_box(data)

        # -- Fit spline --

        #f = interpolate.interp2d(x_values, y_values, aperture_means, kind='cubic')


        #x_grid = np.array(range(self.frame.xsize))
        #y_grid = np.array(range(self.frame.ysize))

        #data = f(x_grid, y_grid)

        # Set new sky frame
        #self.sky = Frame(data)

    # -----------------------------------------------------------------

    def interpolate_apertures(self, aperture_centers, aperture_means):

        """
        This function ...
        :param aperture_centers:
        :param aperture_means:
        :return:
        """

        # Inform the user
        log.info("Interpolating between the mean values of each aperture to fill the sky frame ...")

        x_values = np.array([center.x for center in aperture_centers])
        y_values = np.array([center.y for center in aperture_centers])

        x_ticks = np.arange(0, self.frame.xsize, 1)
        y_ticks = np.arange(0, self.frame.ysize, 1)
        z_grid = mlab.griddata(x_values, y_values, aperture_means, x_ticks, y_ticks)

        # Set the sky frame
        self.sky = Frame(z_grid)

    # -----------------------------------------------------------------

    def interpolate_apertures_try(self, aperture_centers, aperture_means):

        """
        This function ...
        :return:
        """

        return

        x_values = np.array([center.x for center in aperture_centers])
        y_values = np.array([center.y for center in aperture_centers])

        #X, Y = np.meshgrid(x_values, y_values)

        X = x_values
        Y = y_values
        Z = aperture_means

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

    def try_to_interpolate_smart(self, aperture_centers, aperture_means):

        """
        This function ...
        :param aperture_centers:
        :param aperture_means:
        :return:
        """

        model = Pipeline([('poly', PolynomialFeatures(degree=3)), ('linear', LinearRegression(fit_intercept=False))])

        # fit to an order-3 polynomial data
        x = np.arange(5)

        y = 3 - 2 * x + x ** 2 - x ** 3

        model = model.fit(x[:, np.newaxis], y)

    # -----------------------------------------------------------------

    def subtract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the sky from the frame ...")

        # Subtract the estimated sky from the image frame
        self.frame -= self.sky

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
        mask = Mask.from_shape(self.principal_ellipse * factor, self.frame.xsize, self.frame.ysize).inverse()

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
