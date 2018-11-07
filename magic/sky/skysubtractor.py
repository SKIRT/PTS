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
import math
import numpy as np
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D
from photutils.background.core import SExtractorBackground

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ..basics.mask import Mask
from ..region.composite import PixelCompositeRegion
from ..region.list import PixelRegionList, SkyRegionList
from ..tools import plotting, statistics, fitting, plotting
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ..core.mask import Mask as newMask
from ...core.basics.map import Map
from ...core.basics.configuration import save_mapping
from ...core.basics.distribution import Distribution
from ...core.plot.distribution import DistributionPlotter
from pts.core.tools.utils import lazyproperty
from ..region.list import load_as_pixel_region_list

# -----------------------------------------------------------------

estimation_methods = ["mean", "median", "polynomial", "photutils"]
photutils_interpolation_methods = ["mean", "median", "polynomial", "idw"]

# -----------------------------------------------------------------

class SkySubtractor(Configurable):

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
        super(SkySubtractor, self).__init__(*args, **kwargs)

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

        # The region of stars
        self.star_region = None

        # The output mask
        self.mask = None

        # The output mask, but without sigma clipping
        self.mask_not_clipped = None

        # The estimated sky (a single floating point value or a Frame, depending on the estimation method)
        self.sky = None

        # The estimated sky noise
        self.noise = None

        # The statistics
        self.statistics = None

        # The distributions
        self.distributions = Map()

        # Method-specific attributes
        self.poly_interpolated = None
        self.phot_background_mesh = None
        self.phot_background_rms_mesh = None
        self.phot_sky = None
        self.phot_rms = None
        self.phot_boundaries = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Create mask
        self.create_mask()

        # 3. Do an extra sigma-clipping step on the data
        if self.config.sigma_clip_mask: self.sigma_clip()

        # 4. Estimate the sky and sky noise
        self.estimate()

        # 5. Show
        self.show()

        # 6. Set statistics
        self.set_statistics()

        # 7. Create distributions
        if self.config.distributions: self.create_distributions()

        # 8. Write
        if self.config.write: self.write()

        # 9. Plot
        if self.config.plot: self.plot()

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
        if "sources_mask" in kwargs: self.sources_mask = kwargs.pop("sources_mask")
        elif self.config.sources_mask_plane is not None:
            if self.config.image is None: raise ValueError("The image path has to be specified")
            self.sources_mask = Mask.from_file(self.config.image, plane=self.config.sources_mask_plane)

        # NOT REQUIRED ANYMORE?
        if "principal_shape" in kwargs: self.principal_shape = kwargs.pop("principal_shape")
        elif self.config.principal_shape_region is not None: self.principal_shape = load_as_pixel_region_list(self.config.principal_shape_region, self.frame.wcs)[0]

        # Get optional input
        self.extra_mask = kwargs.pop("extra_mask", None)
        self.saturation_region = kwargs.pop("saturation_region", None)
        self.star_region = kwargs.pop("star_region", None)

    # -----------------------------------------------------------------

    @property
    def has_sources_mask(self):
        return self.sources_mask is not None

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sky mask ...")

        # Initialize a list to contain the mask contributions
        masks = []

        # Add sources mask
        if self.has_sources_mask: masks.append(self.sources_mask)

        # Add outside mask
        if self.outside_mask is not None: masks.append(self.outside_mask)

        # Interactive
        if self.config.interactive: plotting.plot_mask(self.outside_mask, title="outside mask")

        # Add saturation mask
        if self.saturation_mask is not None: masks.append(self.saturation_mask)

        # Add stars mask
        if self.stars_mask is not None: masks.append(self.stars_mask)

        # Add the extra mask (if specified)
        if self.extra_mask is not None: masks.append(self.extra_mask)

        # Add the eliminate mask
        if self.eliminate_mask is not None: masks.append(self.eliminate_mask)

        # NEW
        self.mask = newMask.union(*masks)

        # Save as the unclipped mask
        self.mask_not_clipped = self.mask

    # -----------------------------------------------------------------

    def sigma_clip(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing sigma-clipping on the pixel values ...")

        # Create the sigma-clipped mask
        self.mask = statistics.sigma_clip_mask(self.frame, self.config.sigma_clipping.sigma_level, self.mask, self.config.sigma_clipping.niterations)

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

        # Fit
        polynomial = fitting.fit_polynomial(self.frame, self.config.estimation.polynomial.order, mask=self.mask)

        # Evaluate the polynomial
        data = fitting.evaluate_model(polynomial, 0, self.frame.xsize, 0, self.frame.ysize)

        # Create sky map
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

    # -----------------------------------------------------------------

    def estimate_sky_photutils(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the sky and sky noise by using photutils ...")

        # Make cutout
        cutout, mask_cutout = self._make_cutout()

        # Run photutils
        background, background_rms = self._get_photutils_background(cutout, mask_cutout)

        # Make frames
        background, background_rms = self._make_background_frames(background, background_rms)

        # Interpolate the sky frame
        self._interpolate_sky(background)

        # Interpolate the noise frame
        self._interpolate_noise(background_rms)

    # -----------------------------------------------------------------

    def determine_aperture_radius(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the aperture radius ...")

        # Whether or not user has chosen this value
        if self.config.estimation.photutils.aperture_radius is not None: radius = self.config.estimation.photutils.aperture_radius
        else:

            # Check whether the FWHM is defined for the frame
            if self.frame.fwhm is None: raise RuntimeError("The FWHM of the frame is not defined: sky apertures cannot be generated")

            # Determine the radius for the sky apertures
            fwhm_pix = self.frame.fwhm_pix
            radius = self.config.estimation.photutils.aperture_fwhm_factor * fwhm_pix

        # Debugging
        log.debug("Using sky apertures with a radius of " + str(radius) + " pixels")

        # Set the aperture radius
        self.aperture_radius = radius

    # -----------------------------------------------------------------

    def _make_cutout(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Making cutout for the sky estimation ...")

        # Borders
        x_min = self.region_x_min
        x_max = self.region_x_max
        y_min = self.region_y_min
        y_max = self.region_y_max

        x_min = max(int(x_min), 0)
        x_max = min(int(x_max), self.xsize)

        y_min = max(int(y_min), 0)
        y_max = min(int(y_max), self.ysize)

        # Set the phot_boundaries
        self.phot_boundaries = dict()
        self.phot_boundaries["x_min"] = x_min
        self.phot_boundaries["x_max"] = x_max
        self.phot_boundaries["y_min"] = y_min
        self.phot_boundaries["y_max"] = y_max

        # MAKE CUTOUT
        cutout = self.frame[self.cutout_y_slice, self.cutout_x_slice]

        # CUTOUT MASK
        mask_cutout = self.mask[self.cutout_y_slice, self.cutout_x_slice]

        # Return
        return cutout, mask_cutout

    # -----------------------------------------------------------------

    def _make_background_frames(self, background, background_rms):

        """
        This function ...
        :return:
        """

        # Masked background
        background_frame = Frame.nans_like(self.frame)
        background_frame[self.cutout_y_slice, self.cutout_x_slice] = background

        # Masked background rms
        background_rms_frame = Frame.nans_like(self.frame)
        background_rms_frame[self.cutout_y_slice, self.cutout_x_slice] = background_rms

        # Return
        return background_frame, background_rms_frame
    # -----------------------------------------------------------------

    def _get_photutils_background(self, cutout, mask_cutout):

        """
        This function ...
        :return:
        """

        # Determine box size
        if self.config.estimation.photutils.fixed_width is not None:
            box_shape = (self.config.estimation.photutils.fixed_width, self.config.estimation.photutils.fixed_width)
        else:
            # Determine radius
            self.determine_aperture_radius()
            box_shape = (self.integer_aperture_diameter, self.integer_aperture_diameter)

        # Determine filter size
        filter_size = (self.config.estimation.photutils.filter_size, self.config.estimation.photutils.filter_size)

        # NO SIGMA CLIP BECAUSE WE HAVE ALREADY DONE THAT OURSELVES
        sigma_clip = None
        # bkg_estimator = MedianBackground()
        bkg_estimator = SExtractorBackground()
        try:
            bkg = Background2D(cutout, box_shape, filter_size=filter_size, sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator, mask=mask_cutout, filter_threshold=None, #exclude_mesh_method="threshold", (option does not exist anymore)
                               exclude_percentile=self.config.estimation.photutils.exclude_mesh_percentile) # used to be exclude_mesh_percentile
        except ValueError:

            plotting.plot_box(cutout)
            plotting.plot_mask(mask_cutout, title="mask")
            raise RuntimeError("Sky subtraction is not possible for this image")

        # Keep the background 2D object
        self.photutils_bkg = bkg

        #
        self.phot_background_mesh = Frame(bkg.background_mesh)
        self.phot_background_rms_mesh = Frame(bkg.background_rms_mesh)

        # NEW NEW
        self.phot_sky = Frame.nans_like(self.frame)
        self.phot_sky.name = "phot_sky"
        self.phot_sky.description = "photutils background"
        self.phot_sky.unit = self.frame.unit
        self.phot_sky.zero_point = self.frame.zero_point
        self.phot_sky.filter = self.frame.filter
        self.phot_sky.sky_subtracted = False
        self.phot_sky.fwhm = self.frame.fwhm

        self.phot_sky[self.cutout_y_slice, self.cutout_x_slice] = bkg.background

        self.phot_rms = Frame.nans_like(self.frame)
        self.phot_rms.name = "phot_rms"
        self.phot_rms.description = "photutils rms"
        self.phot_rms.unit = self.frame.unit
        self.phot_rms.zero_point = self.frame.zero_point
        self.phot_rms.filter = self.frame.filter
        self.phot_rms.sky_subtracted = False
        self.phot_rms.fwhm = self.frame.fwhm

        self.phot_rms[self.cutout_y_slice, self.cutout_x_slice] = bkg.background_rms

        # Return the background frame and the background rms
        return bkg.background, bkg.background_rms

    # -----------------------------------------------------------------

    def _interpolate_sky(self, background):

        """
        This function ...
        :param background:
        :return:
        """

        # Debugging
        log.debug("Interpolating the sky frame ...")

        # Mean
        if self.config.estimation.photutils.sky_interpolation_method == "mean": self.sky = self._interpolate_mean(background)

        # Median
        elif self.config.estimation.photutils.sky_interpolation_method == "median": self.sky = self._interpolate_median(background)

        # Polynomial
        elif self.config.estimation.photutils.sky_interpolation_method == "polynomial": self.sky = self._interpolate_polynomial(background)

        # IDW (photutils default)
        elif self.config.estimation.photutils.sky_interpolation_method == "idw": self.sky = self._interpolate_idw(background)

        # Invalid
        else: raise ValueError("Invalid interpolation method")

    # -----------------------------------------------------------------

    def _interpolate_noise(self, background_rms):

        """
        This function ...
        :param background_rms:
        :return:
        """

        # Debugging
        log.debug("Interpolating the noise frame ...")

        # Mean
        if self.config.estimation.photutils.noise_interpolation_method == "mean": self.noise = self._interpolate_mean(background_rms)

        # Median
        elif self.config.estimation.photutils.noise_interpolation_method == "median": self.noise = self._interpolate_median(background_rms)

        # Polynomial
        elif self.config.estimation.photutils.noise_interpolation_method == "polynomial": self.noise = self._interpolate_polynomial(background_rms)

        # IDW (photutils default)
        elif self.config.estimation.photutils.noise_interpolation_method == "idw": self.noise = background_rms

        # Invalid
        else: raise ValueError("Invalid interpolation method")

    # -----------------------------------------------------------------

    def _interpolate_mean(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Debugging
        log.debug("Interpolating using mean ...")

        mean_value = np.nanmean(frame.data)
        new = frame.copy()
        new[self.mask] = mean_value
        new.replace_nans(mean_value)
        return new

    # -----------------------------------------------------------------

    def _interpolate_median(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Debugging
        log.debug("Interpolating using median ...")

        median_value = np.nanmedian(frame.data)
        new = frame.copy()
        new[self.mask] = median_value
        new.replace_nans(median_value)
        return new

    # -----------------------------------------------------------------

    def _interpolate_polynomial(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Debugging
        log.debug("Interpolating using polynomial ...")

        masked_data = np.ma.array(frame.data, mask=self.mask.data)

        if self.config.interactive: plotting.plot_box(masked_data)

        mean_value = np.nanmean(masked_data)

        # Normalize the data to have a mean of 1. to avoid miniscule differences
        normalized = frame / mean_value

        # Fit
        polynomial = fitting.fit_polynomial(normalized, self.config.estimation.photutils.polynomial_order, mask=self.mask,
                                            show_warnings=True, fitter=self.config.estimation.photutils.polynomial_fitter)

        #print(polynomial)
        if fitting.all_zero_parameters(polynomial):

            # Debugging
            log.debug("Using the mean value of " + str(mean_value) + " in the estimated sky as the zero order polynomial term ...")

            polynomial = fitting.fit_polynomial(normalized, self.config.estimation.photutils.polynomial_order,
                                                mask=self.mask, show_warnings=True,
                                                fitter=self.config.estimation.photutils.polynomial_fitter, zero_order=1.)

            if fitting.all_zero_parameters(polynomial): raise RuntimeError("Could not fit a polynomial to the data")

        # Evaluate the polynomial
        data = fitting.evaluate_model(polynomial, 0, self.xsize, 0, self.ysize)
        data *= mean_value

        if self.config.interactive: plotting.plot_box(data, "fitted polynomial")

        # Make
        new = frame.copy()

        # Replace
        if self.config.estimation.photutils.replace_all_interpolated: new._data = data
        else: new[self.mask] = data[self.mask]

        # Replace NaNs
        nans = new.nans
        new[nans] = data[nans]

        # Return
        return new

    # -----------------------------------------------------------------

    def _interpolate_idw(self, frame):

        """
        This function ...
        :return:
        """

        mean = np.nanmean(frame)
        return frame.replace_nans(mean)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :param self:
        :return:
        """
        mean_sky = np.nanmean(self.sky)
        median_sky = np.nanmedian(self.sky)
        mean_noise = np.nanmean(self.noise)

        # Debugging
        log.debug("The mean of the sky frame is " + str(mean_sky))
        log.debug("The median of the sky frame is " + str(median_sky))
        log.debug("The mean of the sky noise frame is " + str(mean_noise))

        # Debugging
        log.debug("The mean value before subtraction is " + str(self.mean_frame))
        log.debug("The median value before subtraction is " + str(self.median_frame))
        log.debug("The standard deviation before subtration is " + str(self.stddev_frame))

        # Debugging
        log.debug("The mean value before subtraction without sigma-clipping is " + str(self.mean_frame_not_clipped))
        log.debug("The median value before subtraction without sigma-clipping is " + str(self.median_frame_not_clipped))
        log.debug("The standard deviation before subtraction without sigma-clipping is " + str(self.stddev_frame_not_clipped))

        # Debugging
        log.debug("The mean value after subtraction is " + str(self.mean_subtracted))
        log.debug("The median value after subtraction is " + str(self.median_subtracted))
        log.debug("The standard deviation after subtraction is " + str(self.stddev_subtracted))

        # Debugging
        #log.debug("The estimated sky value is " + str(self.sky))
        #log.debug("The mean value after subtraction is " + str(self.mean_subtracted))
        #log.debug("The median value after subtraction is " + str(self.median_subtracted))
        #log.debug("The standard deviation after subtraction is " + str(self.stddev_subtracted))
        #log.debug("The estimated noise is " + str(self.noise))

        # Debugging
        #log.debug("The mean sky value after subtraction is " + str(self.mean_subtracted))
        #log.debug("The median sky value after subtraction is " + str(self.median_subtracted))
        #log.debug("The global noise is " + str(self.noise))

    # -----------------------------------------------------------------

    def set_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting statistics ...")

        # Set
        self.statistics = Map()
        self.statistics.mean = self.mean_subtracted
        self.statistics.median = self.median_subtracted
        self.statistics.stddev = self.stddev_subtracted

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating distributions ...")

        # Create original distribution
        self.create_original_distribution()

        # Create no-clipping distribution
        self.create_noclipping_distribution()

        # Create subtracted distribution
        self.create_subtracted_distribution()

    # -----------------------------------------------------------------

    def create_original_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating distribution of original pixel values with sigma-clipping ...")

        # Get original pixel values
        original_1d = self.frame_masked_array.compressed()

        # Create distribution of original pixel values
        original = Distribution.from_values("Pixel values", original_1d)

        # Set distribution
        self.distributions.original = original

    # -----------------------------------------------------------------

    def create_noclipping_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating distribution of original pixel values without sigma-clipping ...")

        # Get original pixel values, but not sigma-clipped
        not_clipped_1d = self.not_clipped_masked_array.compressed()

        # Create distribution of pixel values
        not_clipped = Distribution.from_values("Pixel values", not_clipped_1d)

        # Set distribution
        self.distributions.not_clipped = not_clipped

    # -----------------------------------------------------------------

    def create_subtracted_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating distribution of subtracted pixel values ...")

        # Create distribution of subtracted pixel values
        subtracted = Distribution.from_values("Pixel values", self.subtracted_values)

        # Set distribution
        self.distributions.subtracted = subtracted

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write region
        self.write_region()

        # Write mask
        self.write_mask()

        # Write sky
        self.write_sky()

        # Write noise
        self.write_noise()

        # Write subtracted
        self.write_subtracted()

        # Write statistics
        self.write_statistics()

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

    def write_subtracted(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the subtracted frame ...")

        # Determine the path
        path = self.output_path_file("subtracted.fits")

        # Save
        self.subtracted.saveto(path)

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the statistics ...")

        # Determine the path
        path = self.output_path_file("statistics.dat")

        # Save
        save_mapping(path, self.statistics)

    # -----------------------------------------------------------------

    @property
    def photutils(self):
        return self.config.estimation.method == "photutils"

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot histograms
        if self.config.distributions: self.plot_histograms()

        # Plot photutils mesh
        if self.photutils: self.plot_photutils_mesh()

    # -----------------------------------------------------------------

    def plot_histograms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing histograms ...")

        # Write original
        self.write_original_histogram()

        # Write noclipping
        self.write_noclipping_histogram()

        # Write subtracted
        self.write_subtracted_histogram()

    # -----------------------------------------------------------------

    def write_original_histogram(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing original histogram ...")

        # Add distribution
        plotter = DistributionPlotter()
        plotter.add_distribution(self.distributions.original, "original")

        # Determine the path
        path = self.output_path_file("histogram_original.pdf")

        # Run the plotter
        plotter.run(output_path=path)

    # -----------------------------------------------------------------

    def write_noclipping_histogram(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Writing noclipping histogram ...")

        # Add distribution
        plotter = DistributionPlotter()
        plotter.add_distribution(self.distributions.not_clipped, "noclipping")

        # Determine the path
        path = self.output_path_file("histogram_notclipped.pdf")

        # Run the plotter
        plotter.run(output_path=path)

    # -----------------------------------------------------------------

    def write_subtracted_histogram(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing subtracted histogram ...")

        # Add distribution
        plotter = DistributionPlotter()
        plotter.add_distribution(self.distributions.subtracted, "subtracted")

        # Determine the path
        path = self.output_path_file("histogram_subtracted.pdf")

        # Run the plotter
        plotter.run(output_path=path)

    # -----------------------------------------------------------------

    def plot_photutils_mesh(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the photutils mesh ...")

        # Plot meshes
        plt.figure()
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(self.frame, origin='lower', cmap='Greys_r', norm=norm)
        self.photutils_bkg.plot_meshes(outlines=True, color='#1f77b4')

        # Determine path
        path = self.output_path_file("photutils_mesh.pdf")

        # Save
        plt.savefig(path)
        plt.close()

    # -----------------------------------------------------------------

    @lazyproperty
    def integer_aperture_diameter(self):

        """
        This function ...
        :return:
        """

        return 2 * self.integer_aperture_radius

    # -----------------------------------------------------------------

    @lazyproperty
    def integer_aperture_radius(self):

        """
        This function ...
        :return:
        """

        return int(math.ceil(self.aperture_radius))

    # -----------------------------------------------------------------

    @property
    def cutout_xmin(self):

        return self.phot_boundaries["x_min"]

    # -----------------------------------------------------------------

    @property
    def cutout_xmax(self):

        return self.phot_boundaries["x_max"]

    # -----------------------------------------------------------------

    @property
    def cutout_ymin(self):

        return self.phot_boundaries["y_min"]

    # -----------------------------------------------------------------

    @property
    def cutout_ymax(self):

        return self.phot_boundaries["y_max"]

    # -----------------------------------------------------------------

    @property
    def cutout_x_slice(self):

        """
        This function ...
        :return:
        """

        return slice(self.cutout_xmin, self.cutout_xmax)

    # -----------------------------------------------------------------

    @property
    def cutout_y_slice(self):

        """
        This function ...
        :return:
        """

        return slice(self.cutout_ymin, self.cutout_ymax)

    # -----------------------------------------------------------------

    @lazyproperty
    def inner_annulus_region(self):

        """
        This function ...
        :return:
        """

        return self.principal_shape * self.config.mask.annulus_inner_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def outer_annulus_region(self):

        """
        This function ...
        :return:
        """

        return self.principal_shape * self.config.mask.annulus_outer_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def inner_annulus_mask(self):

        """
        eg
        :return:
        """

        return self.inner_annulus_region.to_mask(self.frame.xsize, self.frame.ysize)

    # -----------------------------------------------------------------

    @lazyproperty
    def outer_annulus_mask(self):

        """
        geeh
        :return:
        """

        return self.outer_annulus_region.to_mask(self.frame.xsize, self.frame.ysize)

    # -----------------------------------------------------------------

    @property
    def subtraction_mask(self):

        """
        This property ..
        :return:
        """

        if self.config.sky_region is not None: return self.principal_mask
        else: return self.inner_annulus_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sky region ...")

        # If the sky region has to be loaded from file
        if self.config.sky_region is not None:

            sky_region = SkyRegionList.from_file(self.config.sky_region)
            return sky_region.to_pixel(self.frame.wcs)

        # If no region file is given by the user, create an annulus from the principal ellipse
        elif self.principal_shape is not None:

            inner = self.inner_annulus_region.copy()
            outer = self.outer_annulus_region.copy()
            inner.include = False

            # Create the annulus
            annulus = PixelCompositeRegion(outer, inner)

            # Create the sky region consisting of only the annulus
            region = PixelRegionList()
            region.append(annulus)
            return region

        else: raise ValueError("Cannot create region")

    # -----------------------------------------------------------------

    @property
    def region_bounding_box(self):

        """
        This function ...
        :return:
        """

        return self.region.bounding_box if self.region is not None else None

    # -----------------------------------------------------------------

    @property
    def region_x_min(self):

        """
        This function ...
        :return:
        """

        return self.region.x_min if self.region is not None else None

    # -----------------------------------------------------------------

    @property
    def region_x_max(self):

        """
        This function ...
        :return:
        """

        return self.region.x_max if self.region is not None else None

    # -----------------------------------------------------------------

    @property
    def region_y_min(self):

        """
        This function ...
        :return:
        """

        return self.region.y_min if self.region is not None else None

    # -----------------------------------------------------------------

    @property
    def region_y_max(self):

        """
        This function ...
        :return:
        """

        return self.region.y_max if self.region is not None else None

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        return self.frame.xsize

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        return self.frame.ysize

    # -----------------------------------------------------------------

    @lazyproperty
    def outside_mask(self):

        """
        This function ...
        :return:
        """

        # If region is defined
        if self.region is not None:

            # Create a mask from the pixels outside of the sky region
            outside_mask = self.region.to_mask(self.frame.xsize, self.frame.ysize).inverse()
            return outside_mask

        # Region not defined
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def principal_mask(self):

        """
        This function ...
        :return:
        """

        # If principal shape is provided
        if self.principal_shape is not None:

            # Create a mask from the principal shape
            principal_mask = self.principal_shape.to_mask(self.frame.xsize, self.frame.ysize)
            #masks.append(principal_mask)
            return principal_mask

        # Principal shape not defined
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def saturation_mask(self):

        """
        This function ...
        :return:
        """

        # Check whether saturation contours are defined
        if self.saturation_region is not None:

            # Expand all contours
            expanded_region = self.saturation_region * self.config.mask.saturation_expansion_factor

            # Create the saturation mask
            saturation_mask = expanded_region.to_mask(self.frame.xsize, self.frame.ysize)
            #self.mask += saturation_mask

            return saturation_mask

        # Saturation region is not provided
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def stars_mask(self):

        """
        This function ...
        :return:
        """

        # If star region is defined
        if self.star_region is not None:

            # Expand
            expanded_region = self.star_region * self.config.mask.stars_expansion_factor

            # Create the mask
            stars_mask = expanded_region.to_mask(self.frame.xsize, self.frame.ysize)

            # Return the mask
            return stars_mask

        # Star region not provided
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def eliminate_mask(self):

        """
        gege
        :return:
        """

        if self.config.eliminate is not None:

            reg = load_as_pixel_region_list(self.config.eliminate, wcs=self.frame.wcs)
            return reg.to_mask(self.frame.xsize, self.frame.ysize)

        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def sky_frame(self):

        """
        This function ...
        :return:
        """

        # Get the frame
        if isinstance(self.sky, Frame): result = self.sky
        else: result = Frame(np.full(self.frame.shape, self.sky))

        # Mask
        if self.extra_mask is not None and self.config.add_extra_mask:
            result[self.extra_mask] = self.config.mask_value

        # Return the result
        return result

    # -----------------------------------------------------------------

    @lazyproperty
    def nunmasked_pixels(self):

        """
        This function ....
        :return:
        """

        # Determine the number of unmasked pixels
        npixels = np.sum(self.mask.inverse())
        return npixels

    # -----------------------------------------------------------------

    @property
    def noise_frame(self):

        """
        This function ...
        :return:
        """

        # Get the frame
        if isinstance(self.noise, Frame): result = self.noise
        else: result = Frame(np.full(self.frame.shape, self.noise))

        # Mask
        if self.extra_mask is not None and self.config.add_extra_mask:
            result[self.extra_mask] = self.config.mask_value

        # Return the result
        return result

    # -----------------------------------------------------------------

    @property
    def mean(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped mean
        return np.ma.mean(np.ma.masked_array(self.frame.data, mask=self.mask.data))

    # -----------------------------------------------------------------

    @property
    def median(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped median
        return np.median(np.ma.masked_array(self.frame.data, mask=self.mask.data).compressed())

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

        # Subtract
        result = self.frame - self.sky

        # Set zero outside
        if self.config.set_zero_outside: result[self.subtraction_mask.inverse()] = 0.0

        # Eliminate negatives
        if self.config.eliminate_negatives: result[result <= 0.] = 0.0

        # mask extra
        if self.extra_mask is not None and self.config.add_extra_mask:

            # Set to zero
            result[self.extra_mask] = self.config.mask_value

        # Return the result
        return result

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted_compressed(self):
        return self.subtracted_masked_array.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_masked_array(self):
        return np.ma.masked_array(self.frame.data, mask=self.mask.data)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_compressed(self):
        return self.frame_masked_array.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_frame(self):
        return np.nanmean(self.frame_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def median_frame(self):
        return np.nanmedian(self.frame_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_frame(self):
        return np.nanstd(self.frame_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def not_clipped_masked_array(self):
        return np.ma.masked_array(self.frame.data, mask=self.mask_not_clipped.data)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_not_clipped_compressed(self):
        return self.not_clipped_masked_array.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_frame_not_clipped(self):
        return np.nanmean(self.frame_not_clipped_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def median_frame_not_clipped(self):
        return np.nanmedian(self.frame_not_clipped_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_frame_not_clipped(self):
        return np.nanstd(self.frame_not_clipped_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted_nans(self):
        return np.isnan(self.subtracted.data)

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted_nans_masked_array(self):
        return np.ma.masked_array(self.subtracted_nans, mask=self.mask.data)

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted_nans_compressed(self):
        return self.subtracted_nans_masked_array.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted_masked_array(self):
        return np.ma.masked_array(self.subtracted.data, mask=self.mask.data)

    # -----------------------------------------------------------------

    @lazyproperty
    def mean_subtracted(self):
        #return np.ma.mean(np.ma.masked_array(self.subtracted, mask=self.mask.data))
        return np.nanmean(self.subtracted_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def median_subtracted(self):
        #return np.median(np.ma.masked_array(self.subtracted, mask=self.mask.data).compressed())
        #masked_array = np.ma.masked_array(self.subtracted, mask=self.mask.data)
        return np.nanmedian(self.subtracted_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_subtracted(self):
        # return np.ma.masked_array(self.subtracted, mask=self.mask.data).std()
        return np.nanstd(self.subtracted_compressed)

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted_values(self):

        """
        This function ...
        :return: 
        """

        # Get the subtracted pixel values
        subtracted_1d = self.subtracted_masked_array.compressed()

        # Remove nans
        nans = self.subtracted_nans_compressed

        # Get
        subtracted_values = np.ma.masked_array(subtracted_1d, mask=nans).compressed()

        # Return
        return subtracted_values

# -----------------------------------------------------------------
