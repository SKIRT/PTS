#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.photometry.aperturenoise Contains the ApertureNoiseCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import gc
import pdb
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage
import scipy.stats
from scipy.optimize import curve_fit
from matplotlib import font_manager
from skimage.measure import block_reduce

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import time, tables
from ..tools import plotting
from ...core.basics.distribution import Distribution
from ...core.plot.distribution import DistributionPlotter
from ..core.mask import Mask
from ..dist_ellipse import distance_ellipse
from ..core.frame import Frame
from ..core.segmentationmap import SegmentationMap
from ..basics.vector import Position
from ..basics.coordinate import PixelCoordinate
from ..region.list import PixelRegionList
from ..region.circle import PixelCircleRegion
from ..region.composite import PixelCompositeRegion
from ..core.source import Source
from ..misc import chrisfuncs
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

# Define how many random aperture are desired/required/permitted
sky_success_min = 20  # the minimum number of apertures for the specified size
sky_gen_max = 200  # max number of attempts at generations of a coordinate for each aperture
sky_success_target = 50  # the desired number of apertures of the specified size

# -----------------------------------------------------------------

class ApertureNoiseCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ApertureNoiseCalculator, self).__init__(*args, **kwargs)

        # INPUT

        self.cutout = None
        self.band_name = None
        self.adj_semimaj_pix = None
        self.adj_axial_ratio = None
        self.adj_angle = None
        self.centre_i = None
        self.centre_j = None
        #self.mini_ap_rad_pix = None
        self.downsample_factor = None

        self.annulus_inner_factor = None
        self.annulus_outer_factor = None

        # OUTPUT

        # The aperture noise
        self.noise = None

        # The calculators
        self.exact_calculator = None
        self.extrapolation_calculator = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Try the exact method
        success = self.try_exact()

        # Try the extrapolation method
        if not success:

            # Debugging
            log.debug("Unable to estimate aperture noise using full-size randomly-placed sky apertures (only " + str(int(self.exact_calculator.sky_success_counter)) + " could be placed); switching to aperture extrapolation.")

            # Try the extrapolation method
            success = self.try_extrapolation()

        # If nothing was successful
        if success: log.info("Final aperture noise is " + str(chrisfuncs.ToPrecision(self.noise, 4)) + " (in map units).")
        else: raise RuntimeError("Could not determine the aperture noise")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ApertureNoiseCalculator, self).setup()

        # Get input
        self.cutout = kwargs.pop("cutout")
        self.band_name = kwargs.pop("band_name")
        self.adj_semimaj_pix = kwargs.pop("adj_semimaj_pix")
        self.adj_axial_ratio = kwargs.pop("adj_axial_ratio")
        self.adj_angle = kwargs.pop("adj_angle")
        self.centre_i = kwargs.pop("centre_i")
        self.centre_j = kwargs.pop("centre_j")
        self.downsample_factor = kwargs.pop("downsample_factor")

        self.annulus_inner_factor = kwargs.pop("annulus_inner_factor")
        self.annulus_outer_factor = kwargs.pop("annulus_outer_factor")

        # Set the calculators
        self.exact_calculator = ExactApertureNoiseCalculator()
        self.extrapolation_calculator = ExtrapolatingApertureNoiseCalculator()

        # Set config options
        self.exact_calculator.config.debug_plotting = self.config.debug_plotting
        self.exact_calculator.config.plot_path = self.config.plot_path
        self.extrapolation_calculator.config.debug_plotting = self.config.debug_plotting
        self.extrapolation_calculator.config.plot_path = self.config.plot_path

        self.exact_calculator.config.method = self.config.method # pts or caapr
        self.extrapolation_calculator.config.method = self.config.method # pts or caapr

    # -----------------------------------------------------------------

    def try_exact(self):

        """
        This function ...
        :return:
        """

        # Attempt to determine aperture noise the preferred way, using full-size randomly-placed apertures
        log.info("Estimating aperture noise using full-size randomly-placed sky apertures ...")

        # Set the input
        input_dict = dict()
        input_dict["cutout"] = self.cutout.copy()
        input_dict["band_name"] = self.band_name
        input_dict["adj_semimaj_pix"] = self.adj_semimaj_pix
        input_dict["adj_axial_ratio"] = self.adj_axial_ratio
        input_dict["adj_angle"] = self.adj_angle
        input_dict["centre_i"] = self.centre_i
        input_dict["centre_j"] = self.centre_j
        input_dict["downsample_factor"] = self.downsample_factor

        input_dict["annulus_inner_factor"] = self.annulus_inner_factor
        input_dict["annulus_outer_factor"] = self.annulus_outer_factor

        # Run the calculator
        self.exact_calculator.run(**input_dict)

        # If the function succeeded
        if self.exact_calculator.success:

            # Debugging
            log.debug("Aperture noise successfully estimated using full-size randomly-placed sky apertures")

            # Set the noise
            self.noise = self.exact_calculator.noise

            # Return that we have succeeded
            return True

        # Not succeeded
        else: return False

    # -----------------------------------------------------------------

    def try_extrapolation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating aperture noise using extrapolation of noise estimated from smaller sky apertures ...")

        # Set the input
        input_dict = dict()
        input_dict["cutout"] = self.cutout.copy()
        input_dict["band_name"] = self.band_name
        input_dict["adj_semimaj_pix"] = self.adj_semimaj_pix
        input_dict["adj_axial_ratio"] = self.adj_axial_ratio
        input_dict["adj_angle"] = self.adj_angle
        input_dict["centre_i"] = self.centre_i
        input_dict["centre_j"] = self.centre_j
        input_dict["downsample_factor"] = self.downsample_factor

        input_dict["annulus_inner_factor"] = self.annulus_inner_factor
        input_dict["annulus_outer_factor"] = self.annulus_outer_factor

        # Run the calculator
        self.extrapolation_calculator.run(**input_dict)

        # If even aperture extrapolation is unable to produce an aperture noise estimate, report null value
        if not self.extrapolation_calculator.success: return False
        else:
            self.noise = self.extrapolation_calculator.noise
            return True

# -----------------------------------------------------------------

class ExactApertureNoiseCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ExactApertureNoiseCalculator, self).__init__(*args, **kwargs)

        # INPUT

        self.cutout = None
        self.band_name = None
        self.adj_semimaj_pix = None
        self.adj_axial_ratio = None
        self.adj_angle = None
        self.centre_i = None
        self.centre_j = None
        self.mini_ap_rad_pix = None
        self.downsample_factor = None

        self.annulus_inner_factor = None
        self.annulus_outer_factor = None

        # Is set by the object itself
        self.adj_semimaj_pix_full = None


        # OUTPUT

        self.success = False
        self.noise = None
        self.napertures = None
        self.prior_mask = None
        self.flag_mask = None
        self.sky_success_counter = None
        self.cutout_inviolate = None


        ##

        self.apertures_frame = None
        self.apertures_sum_frame = None
        self.apertures_mean_frame = None
        self.apertures_noise_frame = None

        # Region for the sky aperture circles
        self.aperture_region = PixelRegionList()

        self.covering_apertures = None
        self.apertures_mask = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Function that attempts to estimate aperture noise using randomly-positioned sky apertures of given dimensions
        :param kwargs:
        :return:
        """

        # 2. Calculate the noise
        self.calculate()

    # -----------------------------------------------------------------

    def setup(self, **input_dict):

        """
        This function .
        cutout:                 Array upon which photometry is being perfomred upon
        adj_semimaj_pix:        Semi-major axis of photometric aperture, in pixels.
        adj_axial_ratio:        Axial ratio of photometryic aperture.
        adj_angle:              Position angle of photometric aperture, in degrees.
        adj_centre_i:           Zero-indexed, 0th-axis coordinate (equivalent to y-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.
        adj_centre_j:           Zero-indexed, 1st-axis coordinate (equivalent to x-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.
        :return:
        """

        # Call the setup function of the base class
        super(ExactApertureNoiseCalculator, self).setup()

        # Set input
        self.cutout = input_dict.pop("cutout")
        self.band_name = input_dict.pop("band_name")
        self.adj_semimaj_pix = input_dict.pop("adj_semimaj_pix")
        self.adj_axial_ratio = input_dict.pop("adj_axial_ratio")
        self.adj_angle = input_dict.pop("adj_angle")
        self.centre_i = input_dict.pop("centre_i")
        self.centre_j = input_dict.pop("centre_j")
        self.mini_ap_rad_pix = input_dict.pop("mini_ap_rad_pix", None)
        self.downsample_factor = input_dict.pop("downsample_factor")

        self.annulus_inner_factor = input_dict.pop("annulus_inner_factor")
        self.annulus_outer_factor = input_dict.pop("annulus_outer_factor")

    # -----------------------------------------------------------------

    @property
    def mini(self):

        """
        This function ...
        :return:
        """

        return self.mini_ap_rad_pix is not None

    # -----------------------------------------------------------------

    @property
    def downsample(self):

        """
        This function ...
        :return:
        """

        return int(self.downsample_factor) > 1

    # -----------------------------------------------------------------

    @property
    def mini_full_ratio(self):

        """
        This function ...
        :return:
        """

        ratio = self.adj_semimaj_pix / self.adj_semimaj_pix_full
        return ratio

    # -----------------------------------------------------------------

    def generate_positions(self, adj_semimin_pix, adj_semimin_pix_full, adj_semimaj_pix_full,
                           cutout_inviolate, sky_border):

        """
        This function ...
        :return:
        """

        # Generate random polar coordinates to draw from
        log.debug('Setup: Generating pool of random polar coordinates')

        random_size = sky_success_target * sky_gen_max * 10
        #random_size = 100
        random_failed = []
        random_theta_list = 360.0 * np.random.rand(random_size)
        random_r_list = adj_semimin_pix + np.abs(np.random.normal(loc=0.0, scale=5.0 * adj_semimaj_pix_full, size=random_size))

        # Distribution of radii
        initial_r_distribution = Distribution.from_values("Radius", random_r_list)
        initial_theta_distribution = Distribution.from_values("Angle", random_theta_list)

        # Locate contiguous map regions
        log.debug('Pruning: Locating contiguous coverage regions')

        cont_binary = np.zeros(self.cutout.shape)
        cont_binary[np.where(np.isnan(cutout_inviolate) == False)] = 1
        cont_structure = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        cont_label = scipy.ndimage.measurements.label(cont_binary, structure=cont_structure)[0]
        cont_search_mask = chrisfuncs.EllipseMask(cont_label, 3.0, 1.0, 0.0, self.centre_i, self.centre_j)
        cont_search_values = cont_label[np.where(cont_search_mask == 1)]

        #plotting.plot_mask(cont_binary, title="cont binary")
        #plotting.plot_box(cont_label, title="cont label")
        #plotting.plot_mask(cont_search_mask, title="cont search mask")
        #plotting.plot_mask(cont_search_values, title="cont search values")
        #print(cont_search_values)

        # Identify contiguous region associated with target source
        log.debug('Pruning: Identifying coverage region associated with source')

        if np.where(cont_search_values > 0)[0].shape[0] == 0: cont_target = 0
        else: cont_target = scipy.stats.mode(cont_search_values[np.where(cont_search_values > 0)])[0][0]
        cont_where_bad = np.where(cont_label != cont_target)

        # Remove random coordinates that are more distant than most distant part of the coverage region the target source lies in
        log.debug('Pruning: Removing random coords that definitely lie outside coverage region')

        cont_size_i, cont_size_j = self.cutout.shape
        cont_range_i = np.arange(cont_size_i) - self.centre_i
        cont_range_j = np.arange(cont_size_j) - self.centre_j
        cont_coord_i, cont_coord_j = np.meshgrid(cont_range_j, cont_range_i)  # Yes, i and j are supposed to be this way around inside meshgrid (for some reason)
        cont_coord_i[cont_where_bad] = np.NaN
        cont_coord_j[cont_where_bad] = np.NaN
        cont_dist = np.sqrt(cont_coord_i ** 2 + cont_coord_j ** 2)
        cont_dist_max = np.nanmax(cont_dist)
        random_r_coverage = np.where(random_r_list < (cont_dist_max - sky_border))
        random_r_list = random_r_list[random_r_coverage]
        random_theta_list = random_theta_list[random_r_coverage]

        final_r_distributution = Distribution.from_values("Radius", random_r_list)
        final_theta_distribution = Distribution.from_values("Angle", random_theta_list)

        plotter = DistributionPlotter()
        plotter.add_distribution(initial_r_distribution, "initial r")
        plotter.add_distribution(final_r_distributution, "final r")
        #plotter.run()

        plotter.clear()
        plotter.add_distribution(initial_theta_distribution, "initial theta")
        plotter.add_distribution(final_theta_distribution, "final theta")
        #plotter.run()

        # Convert random polar coordinates into cartesian coordinates
        log.debug('Pruning: Converting random polar coords to cartesian coords, and removing those that lie beyond map border')

        random_i_list = self.centre_i + (random_r_list * np.cos(np.radians(random_theta_list)))  # np.random.normal(loc=centre_i, scale=2.0*sky_ap_rad_pix)
        random_j_list = self.centre_j + (random_r_list * np.sin(np.radians(random_theta_list)))  # np.random.normal(loc=centre_j, scale=2.0*sky_ap_rad_pix)

        # Remove random coodinates that fall fall beyond border in i-coords
        random_not_i_border = np.where((random_i_list > sky_border) & (random_i_list < (self.cutout.shape[0] - sky_border)))
        random_i_list = random_i_list[random_not_i_border]
        random_j_list = random_j_list[random_not_i_border]

        # Remove random coodinates that fall fall beyond border in j-coords
        random_not_j_border = np.where((random_j_list > sky_border) & (random_j_list < (self.cutout.shape[1] - sky_border)))
        random_i_list = random_i_list[random_not_j_border]
        random_j_list = random_j_list[random_not_j_border]

        # Remove random coordinates that intersect source
        log.debug('Pruning: Removing random coords that intersect source')
        random_not_source = np.where(np.sqrt((np.abs(self.centre_i - random_i_list)) ** 2.0 + (
        np.abs(self.centre_j - random_j_list)) ** 2.0) > adj_semimin_pix_full)
        # random_not_source = np.where( (abs(centre_i-random_i_list)>adj_semimin_pix_full) & (abs(centre_j-random_j_list)>adj_semimin_pix_full) )
        random_i_list = random_i_list[random_not_source]
        random_j_list = random_j_list[random_not_source]

        # Remove random coordinates that correspond to NaN pixels
        log.debug('Pruning: Removing random coords that correspond to NaN pixels')
        random_i_list_pix = np.round(random_i_list).astype(int)
        random_j_list_pix = np.round(random_j_list).astype(int)
        random_ij_values = self.cutout[(random_i_list_pix, random_j_list_pix)]
        random_ij_pix_good = np.where(np.isnan(random_ij_values) == False)
        random_i_list = random_i_list[random_ij_pix_good]
        random_j_list = random_j_list[random_ij_pix_good]

        # Return the random coordinates
        return random_i_list, random_j_list, random_failed

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating ...")

        # Handle downsampling input
        if self.downsample: ds_request = int(self.downsample_factor)
        else: ds_request = 1

        # Standard downsampling target is for aperture diameter to correspond to 200 pixels
        ds_target = int( np.round( float(self.adj_semimaj_pix) ) / 100.0 )
        ds_factor = max([ ds_request, ds_target ])

        # If mini-apertures are being used, ensure downsampling isn't agressive (ie, apertures would be sub-Nyquist sampled)
        if self.mini:

            ds_mini = int(np.round(float(self.mini_ap_rad_pix - 1.0) ) / 2.355)

            if ds_mini >= 2: ds_factor = min([ds_factor, ds_mini])
            else: ds_factor = 1

        # If downsampling is makes sense, apply
        if ds_factor >= 2:

            log.debug('Setup: Downsampling map by factor of ' + str(ds_factor))

            self.cutout = block_reduce(self.cutout, block_size=(int(ds_factor),int(ds_factor)), func=np.mean, cval=np.NaN)
            self.cutout *= float(ds_factor) * float(ds_factor)
            self.centre_i /= float(ds_factor)
            self.centre_j /= float(ds_factor)
            self.adj_semimaj_pix /= float(ds_factor)

            if self.mini: self.mini_ap_rad_pix /= float(ds_factor)

        # Handle input variables if mini-apertures are required
        if self.mini:

            log.debug('Setup: Preparing inputs for mini-apertures')

            if isinstance(self.mini_ap_rad_pix, float) or isinstance(self.mini_ap_rad_pix, int):

                mini = float(self.mini_ap_rad_pix)
                self.adj_semimaj_pix_full = self.adj_semimaj_pix
                self.adj_semimaj_pix = mini

            else: pdb.set_trace()

        # Otherwise, adj_semimaj_pix_full = self.adj_semimaj_pix
        else: self.adj_semimaj_pix_full = self.adj_semimaj_pix

        adj_semimin_pix_full = self.adj_semimaj_pix_full / self.adj_axial_ratio

        # Define characteristics of circular aperture with same area as elliptical source aperture
        ap_area = np.pi * self.adj_semimaj_pix * (self.adj_semimaj_pix / self.adj_axial_ratio)
        sky_ap_rad_pix = ( ap_area / np.pi )**0.5
        sky_border = int(sky_ap_rad_pix + 1.0) #int( ( band_dict['annulus_outer'] * sky_ap_rad_pix ) + 1 )
        adj_semimin_pix = self.adj_semimaj_pix / self.adj_axial_ratio

        # ..
        #if self.mini:
        #    print("COMPARISON OF SELF.ADJ_SEMIMAJ_PIX and SKY_AP_RAD_PIX:", self.adj_semimaj_pix, sky_ap_rad_pix)
        #    #assert self.adj_semimaj_pix == sky_ap_rad_pix

        # Creating mask maps to describe no-go regions
        log.debug('Setup: Creating mask maps')

        exclude_mask = chrisfuncs.EllipseMask(self.cutout, self.adj_semimaj_pix_full, self.adj_axial_ratio, self.adj_angle, self.centre_i, self.centre_j)

        # Plot exclude mask
        if self.config.plot_path is not None:
            path = fs.join(self.config.plot_path, "exclude_mask.png")
            plotting.plot_mask(exclude_mask, path=path)

        total_mask = Mask.union(exclude_mask, np.isnan(self.cutout))

        # Determine the maximum theoretical number of circular apertures with this radius for the number of usable pixels in the cutout
        pixel_area = total_mask.nunmasked
        max_number_of_sky_apertures = number_of_apertures_for_radius(sky_ap_rad_pix, pixel_area)
        log.warning("The maximum number of sky apertures is " + str(max_number_of_sky_apertures))

        # Check
        if max_number_of_sky_apertures < sky_success_min:

            log.error("The theoretical maximum number of sky apertures for this image is " + str(max_number_of_sky_apertures) + " , but we need " + str(sky_success_min))

            self.success = False
            self.sky_success_counter = 0

            return

        # Set pixels in source aperture to all have NaN pixels, so they don't get sampled by sky annuli
        cutout_inviolate = self.cutout.copy()
        self.cutout[np.where(chrisfuncs.EllipseMask(self.cutout, self.adj_semimaj_pix_full, self.adj_axial_ratio, self.adj_angle, self.centre_i, self.centre_j) == 1)] = np.NaN


        # Segmentation map for the number of apertures covering each pixel
        self.covering_apertures = SegmentationMap.empty_like(self.cutout)

        # Create a mask that tags all pixels that have been covered by one of the apertures
        self.apertures_mask = Mask.empty_like(self.cutout)

        # Prior mask
        self.prior_mask = Mask.empty_like(self.cutout)

        # Flag mask
        self.flag_mask = np.zeros(self.cutout.shape)

        if self.config.method == "caapr":

            # Chris' method
            # sky_gen_max, adj_semimin_pix, adj_semimin_pix_full, cutout_inviolate, sky_border, sky_ap_rad_pix, exclude_mask
            self.generate_apertures_caapr(adj_semimin_pix, adj_semimin_pix_full, cutout_inviolate, sky_border, sky_ap_rad_pix, exclude_mask, ap_area)

        elif self.config.method == "pts":

            # PTS method
            self.generate_apertures_pts(total_mask, pixel_area, sky_ap_rad_pix, max_number_of_sky_apertures)

        else: raise ValueError("Invalid method (must be 'caapr' or 'pts')")

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def create_aperture_frames(self, aperture_centers, aperture_sums, aperture_means, aperture_stddevs, aperture_radius):

        """
        This function ...
        :param aperture_centers:
        :param aperture_sums:
        :param aperture_means:
        :param aperture_stddevs:
        :param aperture_radius:
        :return:
        """

        self.apertures_frame = Frame.nans_like(self.cutout)
        self.apertures_sum_frame = Frame.nans_like(self.cutout)
        self.apertures_mean_frame = Frame.nans_like(self.cutout)
        self.apertures_noise_frame = Frame.nans_like(self.cutout)

        for i in range(len(aperture_centers)):

            center = aperture_centers[i]

            circle = PixelCircleRegion(center, aperture_radius)

            mask = Mask.from_shape(circle, self.cutout.shape[1], self.cutout.shape[0])

            self.apertures_frame[mask] = self.cutout[mask.data]
            self.apertures_sum_frame[mask] = aperture_sums[i]
            self.apertures_mean_frame[mask] = aperture_means[i]
            self.apertures_noise_frame[mask] = aperture_stddevs[i]

    # -----------------------------------------------------------------

    def generate_apertures_pts(self, total_mask, pixel_area, sky_ap_rad_pix, max_number_of_sky_apertures):

        """
        This function ...
        :param total_mask:
        :param pixel_area:
        :param sky_ap_rad_pix:
        :param max_number_of_sky_apertures:
        :return:
        """

        # DETERMINE THE REQUIRED NUMBER OF APERTURES
        required_napertures = min(int(0.5 * max_number_of_sky_apertures), sky_success_target)

        # SHOW THE REQUIRED NUMBER OF APERTURES
        log.info("The required number of apertures for this radius is " + str(required_napertures))

        #shape = (self.cutout.shape[1], self.cutout.shape[0])
        center = Position(int(round(self.centre_j)), int(round(self.centre_i)))
        ratio = self.adj_axial_ratio
        angle = self.adj_angle * u("deg")

        distance_ell = Frame(distance_ellipse(self.cutout.shape, center, ratio, angle))

        path = fs.join(self.config.plot_path, "distance_ellipse.fits")
        distance_ell.saveto(path)

        # The maximum "major-relative" radius
        max_maj_distance = np.max(distance_ell)


        #circle_area = np.pi * radius ** 2

        # Get arrays of the coordinates of all pixels that are not masked
        #pixels_y, pixels_x = np.where(self.mask.inverse())

        # Get the number of pixels that are not masked (also the area of the frame not masked)
        #npixels = pixels_x.size


        # Counter to keep track of the number of 'succesful' apertures that have been used
        current_napertures = 0

        ngenerations_total = 0
        ngenerations_since_last_aperture = 0

        # Initialize lists to contain the mean sky levels and noise levels in each of the apertures
        aperture_centers = []
        aperture_sums = []
        aperture_means = []
        aperture_stddevs = []

        min_random_r = self.adj_semimaj_pix_full + sky_ap_rad_pix
        max_random_r = max_maj_distance - sky_ap_rad_pix

        #print("MIN RANDOM R:", min_random_r)
        #print("MAX RANDOM R:", max_random_r)

        #exclude_mask = chrisfuncs.Photom.EllipseMask(self.cutout, self.adj_semimaj_pix_full, self.adj_axial_ratio,
        #                                             self.adj_angle, self.centre_i, self.centre_j)
        #attempt_mask = np.zeros(self.cutout.shape)

        # Generate at least 20 apertures, ideally 50
        while True:

            # Draw a random pixel index
            #index = np.random.randint(npixels)

            # Get the x and y coordinate of the pixel
            #x = pixels_x[index]
            #y = pixels_y[index]

            # Inform the user
            log.info("Generating a coordinate for aperture " + str(current_napertures + 1) + " (" + str(required_napertures) + " are required)")

            # GENERATE A RANDOM THETA AND R
            #random_theta_list = 360.0 * np.random.rand(random_size)
            random_theta = 360.0 * np.random.random_sample()
            #random_r_list = adj_semimin_pix + np.abs(np.random.normal(loc=0.0, scale=5.0 * self.adj_semimaj_pix_full, size=random_size))
            random_normalized_r = np.random.uniform(min_random_r, max_random_r)

            unrotated_ellipse_angle = (random_theta - self.adj_angle) * u("deg")

            #print("AXIAL RATIO", self.adj_axial_ratio)
            #print("UNROTATED ELLIPSE ANGLE", unrotated_ellipse_angle)

            radius_at_angle = ellipse_radius_for_angle(1., 1./self.adj_axial_ratio, unrotated_ellipse_angle) # relative to major axis length
            random_real_r = radius_at_angle * random_normalized_r

            random_y = self.centre_i + (random_real_r * np.cos(np.radians(random_theta)))
            random_x = self.centre_j + (random_real_r * np.sin(np.radians(random_theta)))

            x = random_x
            y = random_y

            # Create a coordinate for the center of the aperture
            center = PixelCoordinate(x, y)

            # CHECK WHETHER THE COORDINATE LIES IN THE FRAME
            xsize = self.cutout.shape[1]
            ysize = self.cutout.shape[0]

            if not (0 < int(round(center.x)) < xsize): continue
            if not (0 < int(round(center.y)) < ysize): continue

            # CHECK WHETHER THE COORDINATE IS NOT MASKED

            if total_mask.masks(center): continue

            # COORDINATE IS ACCEPTED

            # we have generated a new coordinate
            ngenerations_total += 1
            ngenerations_since_last_aperture += 1

            # If more than a given number of unsuccessful sky apertures have been generated in a row, call it a day
            if ngenerations_since_last_aperture > sky_gen_max:
                sky_gen_fail = True
                log.debug('Unable to generate suitable random sky aperture after ' + str(ngenerations_since_last_aperture) + ' attempts')
                break

            # CREATE APERTURE

            # Create a circular aperture
            circle = PixelCircleRegion(center, sky_ap_rad_pix)

            # IS THIS APERTURE OK?

            # Create a Source from the frame
            source = Source.from_shape(self.cutout, circle, 1.3)

            # Get a mask of the pixels that overlap with the apertures mask
            apertures_mask_cutout = self.apertures_mask[source.y_slice, source.x_slice]
            overlapping = apertures_mask_cutout * source.mask

            # Calculate the overlap fraction with the apertures mask
            number_of_overlapping_pixels = np.sum(overlapping)
            overlap_fraction = number_of_overlapping_pixels / pixel_area

            # If the overlap fraction is larger than 10% for this aperture, skip it
            if overlap_fraction >= 0.1: log.debug("For this aperture, an overlap fraction of more than 10% was found with other apertures, skipping ...")


            ### CHRIS'S THINGS OF CHECKING THE APERTURE

            # CHRIS' WAY OF CREATING MASK FROM SHAPE
            # Create mask
            #ap_mask = chrisfuncs.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_y, random_x)
            ap_mask = circle.to_mask(self.cutout.shape[1], self.cutout.shape[0])

            #intersection = Mask.intersection(ap_mask, circle.to_mask(self.cutout.shape[1], self.cutout.shape[0]).inverse())
            #print("N INTERSECTION PIXELS:", np.sum(intersection.data.astype(float)))
            #plotting.plot_mask(intersection.data, title="Intersection of Chris and inverse PTS aperture mask (should be empty)")

            # PHOTOMETRY

            # Evaluate pixels in sky aperture
            log.debug('Checking: Evaluating pixels in sky aperture')
            ap_calc = chrisfuncs.EllipseSum(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_y, random_x)

            # Evaluate pixels in sky annulus
            log.debug('Checking: Evaluating pixels in sky annulus')

            #bg_inner_semimaj_pix = self.semimaj_pix_annulus_inner
            #bg_width = self.semimaj_pix_annulus_outer - bg_inner_semimaj_pix
            #print(bg_width)

            bg_inner_semimaj_pix = self.adj_semimaj_pix * self.annulus_inner_factor
            bg_width = (self.adj_semimaj_pix * self.annulus_outer_factor) - bg_inner_semimaj_pix

            bg_width = min(2.0, bg_width)

            #print("ANNULUS INNER FACTOR", self.annulus_inner_factor)
            #print("ANNULUS OUTER FACTOR", self.annulus_outer_factor)
            #print("SELF.ADJ_SEMIMAJ_PIX", self.adj_semimaj_pix)
            #print("BG_INNER_SEMIMAJ_PIX", bg_inner_semimaj_pix)
            #print("BG_WIDTH", bg_width)

            bg_calc = chrisfuncs.AnnulusSum(self.cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_y, random_x)

            # CONTINUED ...

            # Check if more than a given fraction of the pixels inside the source aperture are NaN; if so, reject
            if ap_calc[3][0].shape[0] == 0 or ap_calc[1] == 0: ap_nan_frac = 0.0
            if ap_calc[1] == 0: ap_nan_frac = 1.0
            else: ap_nan_frac = float(ap_calc[3][0].shape[0]) / float(ap_calc[1] + float(ap_calc[3][0].shape[0]))

            ap_nan_thresh = 0.10
            if ap_nan_frac > ap_nan_thresh:

                log.debug('Rejection: Aperture contains too many NaNs')
                #if self.config.debug_plotting.nans: plotting.plot_mask(ap_mask, "Rejection: aperture contains too many NaNs")
                if self.config.debug_plotting.nans: source.plot("Rejection: aperture contains too many NaNs")
                #random_failed.append(random_index)
                continue

            # Check if more than a given fraction of the pixels inside the sky annulus are NaN; if so, reject
            if bg_calc[3][0].shape[0] == 0: bg_nan_frac = 0.0
            if bg_calc[1] == 0: bg_nan_frac = 1.0
            else: bg_nan_frac = float(bg_calc[3][0].shape[0]) / float(bg_calc[1] + bg_calc[3][0].shape[0])

            bg_nan_thresh = 0.80
            if bg_nan_frac > bg_nan_thresh:

                log.debug('Rejection: Annulus contains too many NaNs')
                #if self.config.debug_plotting.annulus_nans: plotting.plot_mask(ap_mask, "Rejection: annulus contains too many NaNs")
                if self.config.debug_plotting.annulus_nans: source.plot("Rejection: annulus contains too many NaNs")
                #random_failed.append(random_index)
                continue

            # If coords have not been rejected for any reason, accept them and proceed
            #else:
            #    sky_success_counter += 1
            #    break


            # APERTURE IS ACCEPTED

            # Add the aperture area to the mask
            self.apertures_mask[source.y_slice, source.x_slice] += source.mask

            # Add to covering mask
            self.covering_apertures.add_shape(circle)

            # Create annulus
            base = PixelCircleRegion(center, bg_inner_semimaj_pix)
            exclude = PixelCircleRegion(center, bg_inner_semimaj_pix + bg_width)
            annulus = PixelCompositeRegion(base, exclude)

            # Add aperture circle to region
            self.aperture_region.append(circle)

            # Add annulus to region
            self.aperture_region.append(annulus)

            # Increment counters
            current_napertures += 1
            ngenerations_since_last_aperture = 0  # we have a new aperture, so set this counter to zero


            # Calculate actual flux in sky aperture, and record
            log.debug('Checking: Performing photometry with random sky aperture and annulus')
            bg_clip = chrisfuncs.SigmaClip(bg_calc[2], median=False, sigma_thresh=3.0)
            bg_avg = bg_clip[1]
            ap_sum = ap_calc[0] - (ap_calc[1] * bg_avg)
            #sky_sum_list.append(ap_sum)


            # Add center to list of center coordinates
            aperture_centers.append(center)

            # Calculate the mean sky value in this aperture
            masked_array_cutout = np.ma.MaskedArray(source.cutout, mask=source.background_mask)
            aperture_mean = np.ma.mean(masked_array_cutout)
            aperture_stddev = np.std(masked_array_cutout)

            # Add mean
            aperture_means.append(aperture_mean)

            # Add sum
            aperture_sums.append(ap_sum)

            # Add
            aperture_stddevs.append(aperture_stddev)


            if np.isnan(ap_sum): pdb.set_trace()

            # Add this aperture to the prior mask and flag mask
            #ap_mask = chrisfuncs.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_y, random_x)
            self.prior_mask += ap_mask
            #flag_mask[np.where(ap_mask == 1)] += 2.0 ** (current_napertures + 1.0)


            ## END

            ##

            # CHECK STOPPING CRITERIA

            # Stop when we have enough apertures
            if current_napertures == required_napertures: break

        # CALCULATE NOISE BASED ON THE APERTURE SUMS (THE SKY PHOTOMETRY APERTURES)

        sky_sum_list = np.array(aperture_sums)
        ap_noise = chrisfuncs.SigmaClip(sky_sum_list, tolerance=0.001, median=True, sigma_thresh=3.0)[0]

        ap_noise = abs(ap_noise)

        log.info("The aperture noise for this radius is " + str(ap_noise))

        if self.config.plot_path is not None:
            path = fs.join(self.config.plot_path, "prior.png")
            plotting.plot_mask(self.prior_mask, path=path)

        # Debugging
        log.debug('Aperture noise from current random apertures is ' + str(chrisfuncs.ToPrecision(ap_noise, 4)) + ' (in map units).')

        self.success = True
        self.noise = ap_noise
        # self.napertures = sky_success_counter

        ####

        # Create aperture frames
        aperture_radius = sky_ap_rad_pix
        self.create_aperture_frames(aperture_centers, aperture_sums, aperture_means, aperture_stddevs, aperture_radius)

        ####

    # -----------------------------------------------------------------

    def generate_apertures_caapr(self, adj_semimin_pix, adj_semimin_pix_full,
                                 cutout_inviolate, sky_border, sky_ap_rad_pix, exclude_mask, ap_area):

        """
        This function ...
        :return:
        """

        # Masks
        #prior_mask = np.zeros(self.cutout.shape)

        attempt_mask = np.zeros(self.cutout.shape)

        # Generate random positions
        random_i_list, random_j_list, random_failed = self.generate_positions(adj_semimin_pix, adj_semimin_pix_full, self.adj_semimaj_pix_full, cutout_inviolate, sky_border)

        # If none of the apertures are suitable, immediately report failure
        if random_i_list.shape[0] == 0:

            log.debug('Status: Pruning removed all generated coordinates')

            self.success = False
            self.sky_success_counter = 0

            self.cutout_inviolate = cutout_inviolate

            return

        # Plot the coordinates
        if self.config.plot_path is not None:
            path = fs.join(self.config.plot_path, "coordinates.png")
            plotting.plot_coordinates_on_image(self.cutout, random_j_list, random_i_list, path=path)

        log.debug("NUMBER OF COORDINATES: " + str(len(random_i_list)))

        # Commence creation of random sky apertures
        sky_success_counter = 0
        sky_sum_list = []
        sky_total_fail = False
        while True:

            # Repeatedly generate random sky apertures, until an acceptable aperture is generated
            sky_gen_counter = 0 # set counter to zero
            sky_gen_fail = False
            while True:

                sky_gen_counter += 1

                # If more than a given number of unsuccessful sky apertures have been generated in a row, call it a day
                if sky_gen_counter > sky_gen_max:

                    sky_gen_fail = True
                    log.debug('Status: Unable to generate suitable random sky aperture after ' + str(sky_gen_max) + ' attempts')
                    break

                # Select random coordinate for this iteration; if no un-used random coordinates can be found, reject
                log.debug('Checking: Selecting random coordinate')
                random_accept = False
                random_reject_count = 0
                while not random_accept:

                    random_index = int(np.floor( np.random.rand() * float(random_i_list.shape[0])))

                    if random_index not in random_failed: random_accept = True
                    else: random_reject_count += 1
                    if random_reject_count > (10 * random_i_list.shape[0]): break

                if random_reject_count > (10 * random_i_list.shape[0]):

                    log.debug('Rejection: Unable to find un-used random coodinates')
                    continue

                # Get the random coordinate
                random_i = random_i_list[random_index]
                random_j = random_j_list[random_index]

                # Create mask
                ap_mask = chrisfuncs.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)

                if log.is_debug:
                    #ap_mask = chrisfuncs.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                    #plotting.plot_mask(ap_mask, title="Aperture " + str(sky_success_counter + 1) + ", generation " + str(sky_gen_counter))
                    attempt_mask[np.where(ap_mask == 1)] = sky_success_counter
                    log.debug('Aperture: ' + str(sky_success_counter + 1) + ';   Generation: ' + str(sky_gen_counter) + ';   Pix Coords: [' + str(random_i) + ',' + str(random_j)+']')

                # Do sophisticated check that generated sky aperture does not intersect source; if it does, reject
                log.debug('Checking: Determining whether aperture intersects source (sophisticated check)')
                exclude_sum = chrisfuncs.EllipseSum(exclude_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)[0]
                if exclude_sum > 0:

                    log.debug('Rejection: Aperture intersects source (according to sophisticated check)')

                    # Plot
                    if self.config.debug_plotting.intersection: plotting.plot_mask(ap_mask, "Rejection: aperture intersects source")

                    random_failed.append(random_index)
                    continue

                # Do basic check that the majority of the pixels in the generated sky aperture have not already been
                # sampled by previous sky apertures; they have, reject
                log.debug('Checking: Determining if aperture over-sampled (basic check)')
                prior_calc = chrisfuncs.EllipseSum(self.prior_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                prior_calc[2][np.where(prior_calc[2] >= 1.0)] = 1.0
                prior_frac = np.sum(prior_calc[2]) / float(prior_calc[1])
                if prior_frac > 0.5:

                    log.debug('Rejection: Aperture over-sampled (according to basic check)')

                    # Plot
                    if self.config.debug_plotting.oversampled: plotting.plot_mask(ap_mask, "Rejection: aperture over-sampled (basic)")

                    random_failed.append(random_index)
                    continue

                # Do sophisticated check that the majority of the pixels in the generated sky aperture have not already been sampled by previous sky apertures; they have, reject
                log.debug('Checking: Determining if aperture oversampled (sophisticated check)')
                ap_mask_check = chrisfuncs.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                flag_mask_check = self.flag_mask.copy()
                flag_mask_check[np.where(ap_mask_check==1)] = int(2.0**(sky_success_counter+1.0))
                flag_tallies = np.array([ np.where(flag_mask_check == flag)[0].shape[0] for flag in (2.0**np.arange(0.0, sky_success_counter+2.0)).tolist() ])
                flag_check = np.where(flag_tallies < 0.5*ap_area)[0].shape[0]
                if flag_check > 1:

                    log.debug('Rejection: Aperture over-sampled (according to sophisticated check)')
                    if self.config.debug_plotting.oversampled: plotting.plot_mask(ap_mask, "Rejection: aperture over-sampled (sophisticated)")
                    random_failed.append(random_index)
                    continue

                # Evaluate pixels in sky aperture
                log.debug('Checking: Evaluating pixels in sky aperture')
                ap_calc = chrisfuncs.EllipseSum(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)

                # Evaluate pixels in sky annulus
                log.debug('Checking: Evaluating pixels in sky annulus')

                #bg_inner_semimaj_pix = self.semimaj_pix_annulus_inner
                #bg_width = self.semimaj_pix_annulus_outer - bg_inner_semimaj_pix

                bg_inner_semimaj_pix = self.adj_semimaj_pix * self.annulus_inner_factor
                bg_width = (self.adj_semimaj_pix * self.annulus_outer_factor) - bg_inner_semimaj_pix

                bg_width = min(2.0, bg_width)

                bg_calc = chrisfuncs.AnnulusSum(self.cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_i, random_j)

                # Check if more than a given fraction of the pixels inside the source aperture are NaN; if so, reject
                if ap_calc[3][0].shape[0] == 0 or ap_calc[1] == 0: ap_nan_frac = 0.0
                if ap_calc[1] == 0: ap_nan_frac = 1.0
                else: ap_nan_frac = float(ap_calc[3][0].shape[0]) / float(ap_calc[1]+float(ap_calc[3][0].shape[0]))
                ap_nan_thresh = 0.10

                if ap_nan_frac > ap_nan_thresh:

                    log.debug('Rejection: Aperture contains too many NaNs')
                    if self.config.debug_plotting.nans: plotting.plot_mask(ap_mask, "Rejection: aperture contains too many NaNs")
                    random_failed.append(random_index)
                    continue

                # Check if more than a given fraction of the pixels inside the sky annulus are NaN; if so, reject
                if bg_calc[3][0].shape[0] == 0: bg_nan_frac = 0.0
                if bg_calc[1] == 0: bg_nan_frac = 1.0
                else: bg_nan_frac = float(bg_calc[3][0].shape[0]) / float(bg_calc[1]+bg_calc[3][0].shape[0])

                bg_nan_thresh = 0.80
                if bg_nan_frac > bg_nan_thresh:

                    log.debug('Rejection: Annulus contains too many NaNs')
                    if self.config.debug_plotting.annulus_nans:
                        plotting.plot_mask(ap_mask, "Rejection: annulus contains too many NaNs")
                        #plotting.
                    random_failed.append(random_index)
                    continue

                # If coords have not been rejected for any reason, accept them and proceed
                else:

                    sky_success_counter += 1
                    break

            # If no suitable sky aperture could be generated on this iteration, decide how to proceed, based on how many had been successfully generated already
            if sky_gen_fail:

                if sky_success_counter < sky_success_min:

                    sky_total_fail = True
                    break

                else:

                    log.debug('Status: However, sufficient number of successful random apertures (' + str(int(sky_success_counter)) + ') already generated; proceeding')
                    break

            # Calculate actual flux in sky aperture, and record
            log.debug('Checking: Performing photometry with random sky aperture and annulus')
            bg_clip = chrisfuncs.SigmaClip(bg_calc[2], median=False, sigma_thresh=3.0)
            bg_avg = bg_clip[1]
            ap_sum = ap_calc[0] - (ap_calc[1] * bg_avg)
            sky_sum_list.append(ap_sum)
            if np.isnan(ap_sum):
                pdb.set_trace()

            # Add this aperture to the prior mask and flag mask
            ap_mask = chrisfuncs.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
            self.prior_mask += ap_mask
            self.flag_mask[np.where(ap_mask==1)] += 2.0**sky_success_counter

            # If target number of sky apertures have been processed, break out of loop
            if sky_success_counter >= sky_success_target:

                log.debug('Status: Target number of successful random apertures (' + str(int(sky_success_target)) + ') generated; proceeding')
                break

        # If total failure was encountered, end process and report now
        if sky_total_fail:

            self.success = False
            self.sky_success_counter = sky_success_counter

            self.cutout_inviolate = cutout_inviolate

            return

        # Otherwise, calculate aperture noise using returned aperture values, and return
        else:

            sky_sum_list = np.array(sky_sum_list)
            ap_noise = chrisfuncs.SigmaClip(sky_sum_list, tolerance=0.001, median=True, sigma_thresh=3.0)[0]

            #chrisfuncs.Cutout(prior_mask, '/home/saruman/spx7cjc/DustPedia/Prior.fits')

            if self.config.plot_path is not None:
                path = fs.join(self.config.plot_path, "prior.png")
                plotting.plot_mask(self.prior_mask, path=path)

            # Debugging
            log.debug('Aperture noise from current random apertures is ' + str(chrisfuncs.ToPrecision(ap_noise,4)) + ' (in map units).')

            self.success = True
            self.noise = ap_noise
            self.napertures = sky_success_counter

            self.cutout_inviolate = cutout_inviolate

            return

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write ...
        self.write_apertures_frame()

        # Write ...
        self.write_apertures_sum_frame()

        # Write ...
        self.write_apertures_mean_frame()

        # Write ...
        self.write_apertures_noise_frame()

        # Write ...
        self.write_aperture_region()

        # Write ...
        self.write_covering_apertures()

        # Write ...
        self.write_apertures_mask()

        # Write ...
        self.write_prior_mask()

        # Write flag mask
        self.write_flag_mask()

    # -----------------------------------------------------------------

    def write_apertures_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the apertures frame ...")

        # Save ...
        apertures_frame_path = fs.join(self.config.plot_path, "apertures.fits")
        self.apertures_frame.saveto(apertures_frame_path)

    # -----------------------------------------------------------------

    def write_apertures_sum_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the apertures sum frame ...")

        # Save ...
        apertures_sum_frame_path = fs.join(self.config.plot_path, "apertures_sum.fits")
        self.apertures_sum_frame.saveto(apertures_sum_frame_path)

    # -----------------------------------------------------------------

    def write_apertures_mean_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the apertures mean frame ...")

        # Save ...
        apertures_mean_frame_path = fs.join(self.config.plot_path, "apertures_mean.fits")
        self.apertures_mean_frame.saveto(apertures_mean_frame_path)

    # -----------------------------------------------------------------

    def write_apertures_noise_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the apertures noise frame ...")

        # Save ...
        apertures_noise_frame_path = fs.join(self.config.plot_path, "apertures_noise.fits")
        self.apertures_noise_frame.saveto(apertures_noise_frame_path)

    # -----------------------------------------------------------------

    def write_aperture_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the aperture region ...")

        # Save region
        region_path = fs.join(self.config.plot_path, "apertures.reg")
        self.aperture_region.saveto(region_path)

    # -----------------------------------------------------------------

    def write_covering_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the covering of the apertures ...")

        # Save covering map
        covering_path = fs.join(self.config.plot_path, "covering.fits")
        self.covering_apertures.saveto(covering_path)

    # -----------------------------------------------------------------

    def write_apertures_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the apertures mask ...")

        # Save aperture mask
        apertures_mask_path = fs.join(self.config.plot_path, "apertures_mask.fits")
        self.apertures_mask.saveto(apertures_mask_path)

    # -----------------------------------------------------------------

    def write_prior_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the prior mask ...")

        # Save prior mask
        prior_mask_path = fs.join(self.config.plot_path, "prior_mask.fits")
        self.prior_mask.saveto(prior_mask_path)

    # -----------------------------------------------------------------

    def write_flag_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the flag mask ...")

        # Save flag mask
        flag_mask_path = fs.join(self.config.plot_path, "flag_mask.fits")
        self.flag_mask.saveto(flag_mask_path)

# -----------------------------------------------------------------

class ExtrapolatingApertureNoiseCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The consturctor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ExtrapolatingApertureNoiseCalculator, self).__init__(*args, **kwargs)

        # INPUT

        self.cutout = None
        self.band_name = None
        self.adj_semimaj_pix = None
        self.adj_axial_ratio = None
        self.adj_angle = None
        self.centre_i = None
        self.centre_j = None
        self.downsample_factor = None

        self.annulus_inner_factor = None
        self.annulus_outer_factor = None

        # OUTPUT

        # The aperture noise
        self.noise = None

        # Success flag
        self.success = False

        # The table with the extrapolation values
        self.extrapolation_table = None

    # -----------------------------------------------------------------

    def run(self, **input_dict):

        """
        This function attempts to estimate the aperture noise using raomy-positioned sky aperturs of given dimensions
        :param input_dict:
        :return:
        """

        # 1. Call the setup function
        self.setup(**input_dict)

        # 2. Calculate
        self.calculate()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **input_dict):

        """
        This function .
        cutout:                 Array upon which photometry is being perfomred upon
        adj_semimaj_pix:        Semi-major axis of photometric aperture, in pixels.
        adj_axial_ratio:        Axial ratio of photometryic aperture.
        adj_angle:              Position angle of photometric aperture, in degrees.
        adj_centre_i:           Zero-indexed, 0th-axis coordinate (equivalent to y-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.
        adj_centre_j:           Zero-indexed, 1st-axis coordinate (equivalent to x-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.
        :return:
        """

        # Call the setup function of the base class
        super(ExtrapolatingApertureNoiseCalculator, self).setup()

        # INPUT
        self.cutout = input_dict.pop("cutout")
        self.band_name = input_dict.pop("band_name")
        self.adj_semimaj_pix = input_dict.pop("adj_semimaj_pix")
        self.adj_axial_ratio = input_dict.pop("adj_axial_ratio")
        self.adj_angle = input_dict.pop("adj_angle")
        self.centre_i = input_dict.pop("centre_i")
        self.centre_j = input_dict.pop("centre_j")
        self.downsample_factor = input_dict.pop("downsample_factor")

        self.annulus_inner_factor = input_dict.pop("annulus_inner_factor")
        self.annulus_outer_factor = input_dict.pop("annulus_outer_factor")

        # Initialize the extrapolation table
        #names = ["Aperture radius", "Aperture area", "Noise"]
        #dtypes = [float, float, float]
        #self.extrapolation_table = Table(names=names, dtype=dtypes, masked=True)

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the aperture noise ...")

        #source_id = source_dict['name']+'_'+band_dict['band_name']

        # Define charactaristics of circular aperture with same area as elliptical source aperture
        ap_area = np.pi * self.adj_semimaj_pix * ( self.adj_semimaj_pix / self.adj_axial_ratio )
        sky_ap_rad_pix = (ap_area / np.pi)**0.5

        # Generate list of mini-aperture sizes to use, and declare result lists
        #mini_ap_rad_base = 1.2 # NEW value used by Chris
        #mini_ap_rad_base = 2.0
        mini_ap_rad_base = 5.0
        mini_ap_rad_pix_input = mini_ap_rad_base**np.arange(1.0, np.ceil( math.log( sky_ap_rad_pix, mini_ap_rad_base)))[::-1]
        min_ap_rad_pix_output = []
        mini_ap_noise_output = []
        mini_ap_num_output = []

        # Loop over radii for mini-apertures
        for mini_ap_rad_pix in mini_ap_rad_pix_input:

            log.debug('Finding aperture noise for mini-apertures of radius ' + str(mini_ap_rad_pix)[:6] + ' pixels.')

            mini_calculator = ExactApertureNoiseCalculator()

            # Set configuration
            mini_calculator.config.debug_plotting = self.config.debug_plotting
            if self.config.plot_path is not None:
                mini_plot_path = fs.join(self.config.plot_path, "mini_" + str(mini_ap_rad_pix))
                fs.create_directory(mini_plot_path)
                mini_calculator.config.plot_path = mini_plot_path

            mini_calculator.config.method = self.config.method # pts or caapr

            input_dict_radius = dict()
            input_dict_radius["cutout"] = self.cutout.copy()
            input_dict_radius["band_name"] = self.band_name
            input_dict_radius["adj_semimaj_pix"] = self.adj_semimaj_pix
            input_dict_radius["adj_axial_ratio"] = self.adj_axial_ratio
            input_dict_radius["adj_angle"] = self.adj_angle
            input_dict_radius["centre_i"] = self.centre_i
            input_dict_radius["centre_j"] = self.centre_j
            input_dict_radius["mini_ap_rad_pix"] = mini_ap_rad_pix
            input_dict_radius["downsample_factor"] = self.downsample_factor

            input_dict_radius["annulus_inner_factor"] = self.annulus_inner_factor
            input_dict_radius["annulus_outer_factor"] = self.annulus_outer_factor

            mini_calculator.run(**input_dict_radius)

            # If mini-aperture succeeded, record and proceed; else, call it a day
            if mini_calculator.success:

                print("MINI CALCULATOR SUCCES, ADDING RADIUS ", mini_ap_rad_pix)

                min_ap_rad_pix_output.append(mini_ap_rad_pix)
                mini_ap_noise_output.append(mini_calculator.noise)
                mini_ap_num_output.append(mini_calculator.napertures)

            # No succes for this radius
            else: log.debug('Unable to place sufficient number of mini-apertures at this radius.')

            # Stop when we have 6 datapoints
            if len(mini_ap_noise_output) >= 6: break

        print("MINI APERTURE RADII WITH SUCCESS:", min_ap_rad_pix_output)

        # Convert output lists into arrays
        min_ap_rad_pix_output = np.array(min_ap_rad_pix_output)
        mini_ap_noise_output = np.array(mini_ap_noise_output)
        mini_ap_num_output = np.array(mini_ap_num_output).astype(float)

        # If insufficient points to make extrapolation, report failure; else proceed
        if min_ap_rad_pix_output.shape[0] < 2:

            self.success = False
            self.noise = None

            gc.collect()

            return

        else:

            mini_ap_areas = np.pi * min_ap_rad_pix_output**2.0

            # Calculate values to plot
            #log_mini_ap_area = np.log10(mini_ap_areas)
            #log_mini_ap_noise = np.log10(mini_ap_noise_output)
            #mini_ap_noise_err_rel = mini_ap_num_output**0.5 / mini_ap_num_output
            ##mini_ap_noise_err = np.abs( mini_ap_noise_output * mini_ap_noise_err_rel )
            #log_mini_ap_noise_err = mini_ap_noise_err_rel #chrisfuncs.LogError(mini_ap_noise_output, mini_ap_noise_err)

            ## NEW: CAAPR COMMIT 4515d19eed03a4e45469531c4fb671e1d5c90206

            # Calculate log of mini aperture area and noise
            log_mini_ap_area = np.log10(np.pi * min_ap_rad_pix_output**2.0)
            log_mini_ap_noise = np.log10(mini_ap_noise_output)

            # Calculate poisson uncertaity on calculated noise
            mini_ap_noise_err_rel = mini_ap_num_output ** 0.5 / mini_ap_num_output
            mini_ap_noise_err = np.abs(mini_ap_noise_output * mini_ap_noise_err_rel)

            # Weight points according to distance in log space from true aperture area
            mini_ap_noise_err *= (1.0 + (np.log10(ap_area) - log_mini_ap_area))

            # Translate uncertainties into log space
            log_mini_ap_noise_err = chrisfuncs.LogError(mini_ap_noise_output, mini_ap_noise_err)

            ##

            # Set the table columns
            names = ["Aperture radius", "Aperture area", "Noise"]
            dtypes = [float, float, float]
            data = [min_ap_rad_pix_output, mini_ap_areas, mini_ap_noise_output]
            self.extrapolation_table = Table(names=names, data=data, dtypes=dtypes)

            log.debug("log mini ap area " + str(log_mini_ap_area))
            log.debug("log mini ap noise " + str(log_mini_ap_noise))
            log.debug("log mini ap noise error " + str(log_mini_ap_noise_err))

            line_fit = curve_fit(Line, log_mini_ap_area, log_mini_ap_noise, sigma=log_mini_ap_noise_err)
            line_m, line_c = line_fit[0][0], line_fit[0][1]

            # Determine projected aperture noise value
            log_ap_area = np.log10(ap_area)
            log_ap_noise_proj = Line(log_ap_area, line_m, line_c) #10.0**( ( line_m * np.log10(ap_area) ) + line_c )
            ap_noise_proj = 10.0**(log_ap_noise_proj)

            # Generate points for best-fit line
            ax_y_min = np.floor( np.min([ np.min( log_mini_ap_noise - log_mini_ap_noise_err ), log_ap_noise_proj ]) )
            ax_y_max = np.ceil( np.max([ np.max( log_mini_ap_noise + log_mini_ap_noise_err ), log_ap_noise_proj ]) )
            ax_x_min = np.floor( np.min([ np.min(log_mini_ap_area), log_ap_area ]) )
            ax_x_max = np.ceil( np.max([ np.max(log_ap_area), log_ap_area ]) )
            line_x = np.linspace( ax_x_min, ax_x_max, num=10000 )
            line_y = Line(line_x, line_m, line_c)

            # Set up figure & axes
            fig = plt.figure(figsize=(8,6))
            ax_dims = [0.125, 0.125, 0.825, 0.825]
            ax = fig.add_axes(ax_dims)

            # Plot points and best-fit line
            ax.errorbar(log_mini_ap_area, log_mini_ap_noise, yerr=log_mini_ap_noise_err, ecolor='#4D78C9', elinewidth=1.15, capthick=1.15, marker='x', color='#0080FF', markersize=0.0, markeredgewidth=1.15, linewidth=0)
            ax.scatter(log_mini_ap_area, log_mini_ap_noise, c='#4D78C9', marker='o', s=75, linewidths=0, label='Mini-aperture noise values')
            ax.scatter(np.log10(ap_area), log_ap_noise_proj, c='#C03030', marker='H', s=150, linewidths=0, label='Extrapolated aperture noise')
            ax.plot(line_x, line_y, ls='--', lw=1.0, c='#4D78C9', label='Line of best fit')

            # Format axis limts and labels
            ax.set_xlabel(r'Aperture Area (log$_{10}$ pix)', fontsize=15)
            ax.set_ylabel(r'Aperture Noise (log$_{10}$ map units)', fontsize=15)
            ax.set_xlim(ax_x_min,ax_x_max)
            ax.set_ylim(ax_y_min,ax_y_max)

            for xlabel in ax.get_xticklabels(): xlabel.set_fontproperties(font_manager.FontProperties(size=15))
            for ylabel in ax.get_yticklabels(): ylabel.set_fontproperties(font_manager.FontProperties(size=15))

            # Plot axis 1 legend
            ax_handles, ax_labels = ax.get_legend_handles_labels()
            ax1_lgd = ax.legend(ax_handles, ax_labels, loc='best', scatterpoints=1, labelspacing=0.25, borderpad=0)
            ax1_lgd.draw_frame(False)
            plt.setp(plt.gca().get_legend().get_texts(), fontsize='12.5')

            if self.config.plot_path is not None:

                # Save figure, clean up, and report results
                fig.savefig(fs.join(self.config.plot_path, time.unique_name("aperture_Noise_Projection_" + self.band_name) + ".png"), dpi=100)
                gc.collect()
                fig.clear()
                plt.close('all')

            self.success = True
            self.noise = ap_noise_proj

            return

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write ...
        self.write_extrapolation_table()

    # -----------------------------------------------------------------

    def write_extrapolation_table(self):

        """
        This function ...
        :return:
        """

        # Determine the path
        path = fs.join(self.config.plot_path, "extrapolation.dat")

        # Write
        tables.write(self.extrapolation_table, path)

# -----------------------------------------------------------------

def number_of_apertures_for_radius(radius, pixel_area):

    """
    This function ...
    :param radius:
    :param pixel_area:
    :return:
    """

    # npixels = np.sum(self.mask.inverse())

    # Assuming optimal hexagonal packing, get an estimate of the maximum number of circles of given radius
    # can fit in the area covered by the pixels that are not masked. This is obviously a significant overestimation
    # especially in the case where the radius becomes of the same order of magnitude as the radius of the
    # galaxy annulus (the hexagonal packing assumes a rectangular area or at least rectangular-like edges)
    # With perfect hexagonal packing, the area of the rectangle that will be covered by the circles is π/(2√3),
    # which is approximately equal to 0.907
    # See: https://www.quora.com/How-many-3-75-inch-circles-will-fit-inside-a-17-inch-square
    coverable_area = 0.907 * pixel_area
    circle_area = np.pi * radius ** 2
    optimal_number_of_apertures = coverable_area / circle_area

    # Debugging
    log.debug("The upper limit to the number of apertures that fit in the part of the frame that is not masked (assuming hexagonal packing) is " + str(optimal_number_of_apertures))

    # Determine the number of apertures that are going to be used, take a third of the upper limit
    # napertures = int(optimal_number_of_apertures / 3.)

    # Debugging
    # log.debug("A total of " + str(napertures) + " apertures are going to be used to estimate the sky ...")

    # Return the number of apertures
    return int(optimal_number_of_apertures)

# -----------------------------------------------------------------

def ellipse_radius_for_angle(a, b, angle):

    """
    This function ...
    :return:
    """

    numerator = a * b
    denominator = math.sqrt( a**2 * math.sin(angle.to("radian").value)**2 + b**2 * math.cos(angle.to("radian").value)**2)

    return numerator / denominator

# -----------------------------------------------------------------

# Define straight-line function, and fit it to points
def Line(x,m,c): return m*x + c

# -----------------------------------------------------------------
