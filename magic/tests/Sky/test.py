#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from itertools import product

# Import astronomical modules
from photutils.datasets import make_100gaussians_image
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import biweight_location
from astropy.stats import biweight_midvariance, mad_std
from astropy.stats import sigma_clipped_stats
#from photutils.detection import detect_sources, detect_threshold # for version 0.2.X
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from photutils import Background2D, SigmaClip, MedianBackground # for 0.3.1
from photutils import SExtractorBackground
# from photutils.background import Background # for version 0.2.x
from astropy.modeling.models import Sersic2D
from astropy.coordinates import Angle
from photutils import make_source_mask
from photutils.datasets import make_random_gaussians, make_noise_image, make_gaussian_sources

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.magic.core.frame import Frame
from pts.magic.core.mask import Mask
from pts.magic.tools import plotting
from pts.magic.region.ellipse import PixelEllipseRegion
from pts.magic.basics.coordinate import PixelCoordinate
from pts.magic.basics.stretch import PixelStretch
from pts.magic.tools import statistics
from pts.core.basics.map import Map
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import save_mapping
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

# https://photutils.readthedocs.io/en/stable/photutils/background.html

# -----------------------------------------------------------------

description = "Test the sky subtraction"

# -----------------------------------------------------------------

class SkyTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SkyTest, self).__init__(*args, **kwargs)

        # FRAME COMPONENTS

        # The sources map
        self.sources = None

        # The noise map
        self.noise = None

        # The galaxy
        self.galaxy = None

        # Real sky frames
        self.constant_sky = None
        self.gradient_sky = None

        # TABLES
        self.source_table = None

        # MASKS

        # The sources mask
        self.sources_mask = None

        # The rotation mask
        self.rotation_mask = None

        # The galaxy
        self.galaxy_region = None

        # SKY ESTIMATION

        # Photutils Background2D object
        self.photutils_bkg = None

        # Sky reference estimation
        self.reference_sky = None

        # Path
        self.subtraction_path = None

        # The sky subtractor
        self.subtractor = None

        # STATISTICS

        # The statistics
        self.statistics = Map()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make rotation mask
        self.make_rotation_mask()

        # 3. Generate the sources
        self.make_sources()

        # 4. Make noise
        self.make_noise()

        # 5. Make galaxy
        self.make_galaxy()

        # 6. Make sky
        self.make_sky()

        # 7. Mask sources
        self.mask_sources()

        # 8. Statistics
        self.calculate_statistics()

        # 9. Reference
        self.reference()

        # 10. Subtract
        self.subtract()

        # 11. Write
        self.write()

        # 12. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SkyTest, self).setup(**kwargs)

        # Set the subtraction path
        self.subtraction_path = fs.create_directory_in(self.path, "subtraction")

    # -----------------------------------------------------------------

    def make_rotation_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making rotation mask ...")

        # Rotate
        if self.config.rotate:

            frame = Frame.zeros(self.config.shape)
            self.rotation_mask = frame.rotate(self.config.rotation_angle)

        else: self.rotation_mask = Mask.empty(self.config.shape[1], self.config.shape[0])

    # -----------------------------------------------------------------

    @property
    def effective_rotation_angle(self):

        """
        THis function ...
        :return:
        """

        if self.config.rotate: return self.config.rotation_angle
        else: return Angle(0.0, "deg")

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        THis function ...
        :return:
        """

        return self.rotation_mask.shape

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.rotation_mask.xsize

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.rotation_mask.ysize

    # -----------------------------------------------------------------

    @property
    def sources_sigma(self):

        """
        This function ...
        :return:
        """

        return self.config.fwhm * statistics.fwhm_to_sigma

    # -----------------------------------------------------------------

    def make_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the sources ...")

        flux_range = [self.config.flux_range.min, self.config.flux_range.max]
        xmean_range = [0, self.config.shape[1]]
        ymean_range = [0, self.config.shape[0]]

        # Ranges of sigma
        xstddev_range = [self.sources_sigma, self.sources_sigma]
        ystddev_range = [self.sources_sigma, self.sources_sigma]

        table = make_random_gaussians(self.config.nsources, flux_range, xmean_range,
                                      ymean_range, xstddev_range,
                                      ystddev_range, random_state=12345)
        self.source_table = table

        data = make_gaussian_sources(self.config.shape, table)
        self.sources = Frame(data)

        # mask
        self.sources[self.rotation_mask] = 0.0

        if self.config.plot: plotting.plot_box(self.sources, title="sources")

    # -----------------------------------------------------------------

    def make_noise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making noise map ...")

        # Make noise
        data = make_noise_image(self.config.shape, type='gaussian', mean=0., stddev=self.config.noise_stddev, random_state=12345)
        self.noise = Frame(data)

        # Mask
        self.noise[self.rotation_mask] = 0.0

        # Plot
        #if self.config.plot: plotting.plot_difference(self.frame, self.real_sky, title="original")
        if self.config.plot: plotting.plot_box(self.noise, title="noise")

    # -----------------------------------------------------------------

    def make_galaxy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding smooth galaxy source ...")

        effective_radius = self.config.galaxy_effective_radius
        effective_galaxy_angle = self.config.galaxy_angle + self.effective_rotation_angle

        axial_ratio = self.config.galaxy_axial_ratio
        angle_deg = effective_galaxy_angle.to("deg").value

        # Produce guess values
        initial_sersic_amplitude = self.config.galaxy_central_flux
        initial_sersic_r_eff = effective_radius
        initial_sersic_n = self.config.galaxy_sersic_index
        initial_sersic_x_0 = self.config.galaxy_position.x
        initial_sersic_y_0 = self.config.galaxy_position.y
        initial_sersic_ellip = (axial_ratio - 1.0) / axial_ratio
        initial_sersic_theta = np.deg2rad(angle_deg)

        # Produce sersic model from guess parameters, for time trials
        sersic_x, sersic_y = np.meshgrid(np.arange(self.xsize), np.arange(self.ysize))
        sersic_model = Sersic2D(amplitude=initial_sersic_amplitude, r_eff=initial_sersic_r_eff,
                                                        n=initial_sersic_n, x_0=initial_sersic_x_0,
                                                        y_0=initial_sersic_y_0, ellip=initial_sersic_ellip,
                                                        theta=initial_sersic_theta)
        sersic_map = sersic_model(sersic_x, sersic_y)

        # Set the galaxy frame
        self.galaxy = Frame(sersic_map)

        # Mask
        self.galaxy[self.rotation_mask] = 0.0

        limit_radius = self.config.galaxy_relative_asymptotic_radius * effective_radius

        # Create galaxy region
        galaxy_center = PixelCoordinate(initial_sersic_x_0, initial_sersic_y_0)
        galaxy_radius = PixelStretch(limit_radius, limit_radius / axial_ratio)
        self.galaxy_region = PixelEllipseRegion(galaxy_center, galaxy_radius, effective_galaxy_angle)

        # Set galaxy map zero outside certain radius
        self.galaxy[self.galaxy_mask.inverse()] = 0.0

        # Plot
        if self.config.plot: plotting.plot_box(self.galaxy, title="galaxy")

    # -----------------------------------------------------------------

    def make_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making sky ...")

        # make constant sky
        self.make_constant_sky()

        # Make gradient sky
        self.make_gradient_sky()

    # -----------------------------------------------------------------

    def make_constant_sky(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Making constant sky ...")

        self.constant_sky = Frame.filled_like(self.sources, self.config.constant_sky)

        # Mask
        self.constant_sky[self.rotation_mask] = 0.0

        # Plot

    # -----------------------------------------------------------------

    def make_gradient_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making gradient sky ...")

        y, x = np.mgrid[:self.ysize, :self.xsize]

        # Create gradient sky
        self.gradient_sky = Frame(x * y / 5000.)

        # Mask padded
        self.gradient_sky[self.rotation_mask] = 0.0

        # Plot
        #if self.config.plot: plotting.plot_difference(self.frame, self.real_sky, title="frame with background")
        if self.config.plot: plotting.plot_box(self.gradient_sky, title="gradient sky")

    # -----------------------------------------------------------------

    @lazyproperty
    def sky(self):

        """
        This function ...
        :return:
        """

        return self.constant_sky + self.gradient_sky

    # -----------------------------------------------------------------

    @lazyproperty
    def frame(self):

        """
        This fucntion ...
        :return:
        """

        return self.sources_with_galaxy_and_noise + self.sky

    # -----------------------------------------------------------------

    @lazyproperty
    def sources_with_noise(self):

        """
        This function ...
        :return:
        """

        return self.noise + self.sources

    # -----------------------------------------------------------------

    @lazyproperty
    def sources_with_galaxy_and_noise(self):

        """
        This function ...
        :return:
        """

        return self.sources_with_noise + self.galaxy

    # -----------------------------------------------------------------

    def mask_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Masking sources ...")

        # Create sources mask
        mask = make_source_mask(self.sources_with_noise.data, snr=2, npixels=5, dilate_size=11, mask=self.rotation_mask)
        self.sources_mask = Mask(mask)

        # Plot
        if self.config.plot: plotting.plot_mask(self.sources_mask, title="sources mask")

    # -----------------------------------------------------------------

    def calculate_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating statistics ...")

        # Calculate statistics no sigma clipping
        self.calculate_statistics_no_clipping()

        # Calculate statistics clipping
        self.calculate_statistics_clipping()

        # Calculate statistics masked sources
        self.calculate_statistics_masked()

    # -----------------------------------------------------------------

    def calculate_statistics_no_clipping(self):

        """
        This function ...
        :return:
        """

        # Compress (remove masked values)
        flattened = np.ma.array(self.sources_with_noise.data, mask=self.rotation_mask.data).compressed()

        median = np.median(flattened)
        biweight_loc = biweight_location(flattened)

        biweight_midvar = biweight_midvariance(flattened)
        median_absolute_deviation = mad_std(flattened)

        #print("median", median)
        #print("biweigth_loc", biweight_loc)
        #print("biweight_midvar", biweight_midvar)
        #print("median_absolute_deviation", median_absolute_deviation)

        self.statistics.no_clipping = Map()
        self.statistics.no_clipping.median = median
        self.statistics.no_clipping.biweight_loc = biweight_loc
        self.statistics.no_clipping.biweight_midvar = biweight_midvar
        self.statistics.no_clipping.median_absolute_deviation = median_absolute_deviation

        # SAME RESULTS:

        # median = np.median(self.original_frame)
        # biweight_loc = biweight_location(self.original_frame)
        # biweight_midvar = biweight_midvariance(self.original_frame)
        # median_absolute_deviation = mad_std(self.original_frame)

        # print("median", median)
        # print("biweigth_loc", biweight_loc)
        # print("biweight_midvar", biweight_midvar)
        # print("median_absolute_deviation", median_absolute_deviation)

    # -----------------------------------------------------------------

    def calculate_statistics_clipping(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Sigma-clipping ...")

        # Sigma clip
        mean, median, std = sigma_clipped_stats(self.sources_with_noise.data, sigma=3.0, iters=5,
                                                mask=self.rotation_mask)

        #print("sigma-clip mean:", mean)
        #print("sigma-clip median:", median)
        #print("sigma-clip std:", std)

        self.statistics.clipping = Map()
        self.statistics.clipping.mean = mean
        self.statistics.clipping.median = median
        self.statistics.clipping.std = std

    # -----------------------------------------------------------------

    def calculate_statistics_masked(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating statistics with sources masked ...")

        # Statistics
        mean, median, std = sigma_clipped_stats(self.sources_with_noise.data, sigma=3.0, mask=self.total_mask.data, iters=5)

        #print("sigma-clip mean after source masking:", mean)
        #print("sigma-clip median after source masking:", median)
        #print("sigma_clip std after source masking:", std)

        # Set the statistics
        self.statistics.masked = Map()
        self.statistics.masked.mean = mean
        self.statistics.masked.median = median
        self.statistics.masked.std = std

    # -----------------------------------------------------------------

    @lazyproperty
    def total_mask(self):

        """
        This function ...
        :return:
        """

        return self.sources_and_rotation_mask + self.galaxy_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def sources_and_rotation_mask(self):

        """
        This function ...
        :return:
        """

        return self.sources_mask + self.rotation_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_mask(self):

        """
        This function ...
        :return:
        """

        # Make mask
        return self.galaxy_region.to_mask(self.xsize, self.ysize)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_subtracted(self):

        """
        This function ...
        :return:
        """

        return self.frame - self.reference_sky

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted(self):

        """
        This function ...
        :return:
        """

        return self.frame - self.estimated_sky

    # -----------------------------------------------------------------

    @property
    def aperture_radius(self):

        """
        This function ...
        :return:
        """

        return self.config.fwhm * self.config.aperture_fwhm_factor

    # -----------------------------------------------------------------

    @property
    def aperture_diameter(self):

        """
        This function ...
        :return:
        """

        return 2.0 * self.aperture_radius

    # -----------------------------------------------------------------

    def reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating background with photutils ...")

        # Plot total mask
        if self.config.plot: plotting.plot_mask(self.total_mask, title="total mask")

        integer_aperture_radius = int(math.ceil(self.aperture_radius))
        box_shape = (integer_aperture_radius, integer_aperture_radius)
        filter_size = (3, 3)

        # Estimate the background
        sigma_clip = SigmaClip(sigma=3., iters=10)
        #bkg_estimator = MedianBackground()
        bkg_estimator = SExtractorBackground()

        bkg = Background2D(self.frame.data, box_shape, filter_size=filter_size, sigma_clip=sigma_clip,
                           bkg_estimator=bkg_estimator, mask=self.total_mask.data)

        # Statistics
        #print("median background", bkg.background_median)
        #print("rms background", bkg.background_rms_median)

        self.statistics.reference = Map()
        self.statistics.reference.median = bkg.background_median
        self.statistics.reference.rms = bkg.background_rms_median

        # Plot
        if self.config.plot: plotting.plot_box(bkg.background, title="background from photutils")

        # Set the sky
        self.reference_sky = Frame(bkg.background)

        # Set bkg object
        self.photutils_bkg = bkg

        # Plot
        if self.config.plot: plotting.plot_box(self.reference_sky, title="reference sky")

    # -----------------------------------------------------------------

    def subtract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the sky ...")

        # Settings
        settings = dict()
        settings["estimation"] = dict()
        settings["estimation"]["method"] = "photutils"
        settings["estimation"]["aperture_radius"] = self.aperture_radius
        settings["write"] = True
        #settings["estimation"]["finishing_step"] = "polynomial"
        #settings["estimation"]["polynomial_degree"] = self.config.polynomial_degree
        #settings["estimation"]["fill_method"] = "cubic"

        settings["plot"] = True

        # Input
        input_dict = dict()

        # Set input
        input_dict["frame"] = self.frame
        input_dict["principal_shape"] = self.galaxy_region
        input_dict["sources_mask"] = self.sources_mask
        input_dict["extra_mask"] = self.rotation_mask

        # Create command
        command = Command("subtract_sky", "subtract the sky from an artificially created image", settings, input_dict, cwd=self.subtraction_path)

        # Run the subtraction
        self.subtractor = self.run_command(command)

    # -----------------------------------------------------------------

    @property
    def estimated_sky(self):

        """
        This function ...
        :return:
        """

        return self.subtractor.sky_frame

    # -----------------------------------------------------------------

    @lazyproperty
    def subtracted(self):

        """
        This function ...
        :return:
        """

        return self.frame - self.estimated_sky

    # -----------------------------------------------------------------

    @lazyproperty
    def sky_residual(self):

        """
        This function ...
        :return:
        """

        return self.estimated_sky - self.sky

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the sources with noise
        self.write_sources()

        # Write the frame
        self.write_frame()

        # Write real sky map
        self.write_real_sky()

        # Write sources mask
        self.write_sources_mask()

        # Write galaxy mask
        self.write_galaxy_mask()

        # Write reference sky
        self.write_reference_sky()

        # Write estiamted sky
        self.write_estimated_sky()

        # Write residuals
        self.write_residual()

        # Write the statistics
        self.write_statistics()

    # -----------------------------------------------------------------

    def write_sources(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sources frame with noise ...")

        # Determine the path
        path = fs.join(self.path, "sources_noise.fits")

        # Save the frame
        self.sources_with_noise.saveto(path)

    # -----------------------------------------------------------------

    def write_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the frame ...")

        # Determine the path
        path = fs.join(self.path, "frame.fits")

        # SAve the frame
        self.frame.saveto(path)

    # -----------------------------------------------------------------

    def write_real_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the real sky map ...")

        # Determine the path
        path = fs.join(self.path, "real_sky.fits")

        # Save the map
        self.sky.saveto(path)

    # -----------------------------------------------------------------

    def write_sources_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sources mask ...")

        # Determine the path
        path = fs.join(self.path, "sources_mask.fits")

        # Save
        self.sources_mask.saveto(path)

    # -----------------------------------------------------------------

    def write_galaxy_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the galaxy mask ...")

        # Determine the path
        path = fs.join(self.path, "galaxy_mask.fits")

        # Save
        self.galaxy_mask.saveto(path)

    # -----------------------------------------------------------------

    def write_reference_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the reference sky map ...")

        # Determine the path
        path = fs.join(self.path, "reference_sky.fits")

        # Save
        self.reference_sky.saveto(path)

    # -----------------------------------------------------------------

    def write_estimated_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the estimated sky map ...")

        # Determine the path
        path = fs.join(self.path, "estimated_sky.fits")

        # Save
        self.estimated_sky.saveto(path)

    # -----------------------------------------------------------------

    def write_subtracted(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sky-subtracted frame ...")

        # Determine the path
        path = fs.join(self.path, "subtracted.fits")

        # Save
        self.subtracted.saveto(path)

    # -----------------------------------------------------------------

    def write_residual(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual map ...")

        # Determine the path
        path = fs.join(self.path, "residual.fits")

        # Save
        self.sky_residual.saveto(path)

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the statistics ...")

        # Determine the path
        path = fs.join(self.path, "statistics.dat")

        # Write
        save_mapping(path, self.statistics)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inofmrthe user
        log.info("Plotting ...")

        # Plot meshes
        if self.config.plotting.meshes: self.plot_meshes()

    # -----------------------------------------------------------------

    def plot_meshes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the meshes ...")

        # Plot meshes
        plt.figure()
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(self.frame, origin='lower', cmap='Greys_r', norm=norm)
        self.photutils_bkg.plot_meshes(outlines=True, color='#1f77b4')
        plt.show()

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    pass

# -----------------------------------------------------------------
