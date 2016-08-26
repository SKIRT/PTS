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
import ChrisFuncs

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ApertureNoiseCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(ApertureNoiseCalculator, self).__init__(config)

        # INPUT

        self.cutout = None
        self.source_dict = None
        self.band_dict = None
        self.kwargs_dict = None
        self.adj_semimaj_pix = None
        self.adj_axial_ratio = None
        self.adj_angle = None
        self.centre_i = None
        self.centre_j = None
        self.downsample_factor = None

        # OUTPUT

        # The aperture noise
        self.noise = None

        # The calculators
        self.exact_calculator = None
        self.extrapolation_calculator = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Try the exact method
        success = self.try_exact()

        # Try the interpolation method
        if not success:

            # Debugging
            log.debug("Unable to estiamte aperture noise using full-size randomly-placed sky apertures (only " + str(int(self.exact_calculator.sky_success_counter)) + " could be placed); switching to aperture extrapolation.")

            # Try the extrapolation method
            success = self.try_extrapolation()

        # If nothing was successful
        if success: log.info("Final aperture noise is " + str(ChrisFuncs.FromGitHub.randlet.ToPrecision(self.noise, 4)) + " (in map units).")
        else: raise RuntimeError("Could not determine the aperture noise")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param input_dict:
        :return:
        """

        # Call the setup function of the base class
        super(ApertureNoiseCalculator, self).setup()

        # Get input
        self.cutout = kwargs.pop("cutout")
        self.source_dict = kwargs.pop("source_dict")
        self.band_dict = kwargs.pop("band_dict")
        self.kwargs_dict = kwargs.pop("kwargs_dict")
        self.adj_semimaj_pix = kwargs.pop("adj_semimaj_pix")
        self.adj_axial_ratio = kwargs.pop("adj_axial_ratio")
        self.adj_angle = kwargs.pop("adj_angle")
        self.centre_i = kwargs.pop("centre_i")
        self.centre_j = kwargs.pop("centre_j")
        self.downsample_factor = kwargs.pop("downsample_factor")

        # Set the calculators
        self.exact_calculator = ExactApertureNoiseCalculator()
        self.extrapolation_calculator = ExtrapolatingApertureNoiseCalculator()

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
        input_dict["source_dict"] = self.source_dict
        input_dict["band_dict"] = self.band_dict
        input_dict["kwargs_dict"] = self.kwargs_dict
        input_dict["adj_semimaj_pix"] = self.adj_semimaj_pix
        input_dict["adj_axial_ratio"] = self.adj_axial_ratio
        input_dict["adj_angle"] = self.adj_angle
        input_dict["centre_i"] = self.centre_i
        input_dict["centre_j"] = self.centre_j
        #input_dict["downsample"] = int(band_dict['downsample_factor'])
        input_dict["downsample_factor"] = self.downsample_factor

        # Run the calculator
        self.exact_calculator.run(input_dict)

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
        input_dict["source_dict"] = self.source_dict
        input_dict["band_dict"] = self.band_dict
        input_dict["kwargs_dict"] = self.kwargs_dict
        input_dict["adj_semimaj_pix"] = self.adj_semimaj_pix
        input_dict["adj_axial_ratio"] = self.adj_axial_ratio
        input_dict["adj_angle"] = self.adj_angle
        input_dict["centre_i"] = self.centre_i
        input_dict["centre_j"] = self.centre_j

        # Run the calculator
        self.extrapolation_calculator.run(**input_dict)

        # If even aperture extrapolation is unable to produce an aperture noise estimate, report null value
        if not self.extrapolation_calculator.success: return False
        else:
            self.noise = self.extrapolation_calculator.noise
            return True

# -----------------------------------------------------------------

class ExactApertureNoiseCalculator(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # INPUT

        self.cutout = None
        self.source_dict = None
        self.band_dict = None
        self.kwargs_dict = None
        self.adj_semimaj_pix = None
        self.adj_axial_ratio = None
        self.adj_angle = None
        self.centre_i = None
        self.centre_j = None
        self.mini = None
        self.downsample_factor = None

        # OUTPUT

        self.success = False
        self.noise = None
        self.napertures = None
        self.prior_mask = None
        self.flag_mask = None
        self.sky_success_counter = None
        self.cutout_inviolate = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        Function that attempts to estimate aperture noise using randomly-positioned sky apertures of given dimensions
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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

        # cutout, source_dict, band_dict, kwargs_dict, adj_semimaj_pix, adj_axial_ratio, adj_angle, centre_i, centre_j, mini=False, downsample=False

        self.cutout = input_dict.pop("cutout")
        self.source_dict = input_dict.pop("source_dict")
        self.band_dict = input_dict.pop("band_dict")
        self.kwargs_dict = input_dict.pop("kwargs_dict")
        self.adj_semimaj_pix = input_dict.pop("adj_semimaj_pix")
        self.adj_axial_ratio = input_dict.pop("adj_axial_ratio")
        self.adj_angle = input_dict.pop("adj_angle")
        self.centre_i = input_dict.pop("centre_i")
        self.centre_j = input_dict.pop("centre_j")
        self.mini = input_dict.pop("mini")
        self.downsample_factor = input_dict.pop("downsample_factor")

    # -----------------------------------------------------------------

    @property
    def downsample(self):

        """
        This function ...
        :return:
        """

        return int(self.downsample_factor)

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        #source_id = source_dict['name'] + '_' + band_dict['band_name']
        ap_debug = False

        # Handle downsampling input
        if self.downsample != False: ds_request = int(self.downsample)
        else: ds_request = 1

        # Standard downsampling target is for aperture diameter to correspod to 200 pixels
        ds_target = int( np.round( float(self.adj_semimaj_pix) ) / 100.0 )
        ds_factor = max([ ds_request, ds_target ])

        # If mini-apertures are being used, ensure downsampling isn't agressive (ie, apertures would be sub-Nyquist sampled)
        if self.mini != False:

            ds_mini = int(np.round( float(self.mini - 1.0) ) / 2.355)

            if ds_mini >= 2: ds_factor = min([ ds_factor, ds_mini ])
            else: ds_factor = 1

        # If downsampling is makes sense, apply
        if ds_factor >= 2:

            if ap_debug: log.debug('Setup: Downsampling map by factor of ' + str(ds_factor))
            self.cutout = block_reduce(self.cutout, block_size=(int(ds_factor),int(ds_factor)), func=np.mean, cval=np.NaN)
            self.cutout *= float(ds_factor)*float(ds_factor)
            self.centre_i /= float(ds_factor)
            self.centre_j /= float(ds_factor)
            self.adj_semimaj_pix /= float(ds_factor)

            if self.mini != False: self.mini /= float(ds_factor)

        # Handle input variables if mini-apertures are required
        if self.mini != False:

            if ap_debug: log.debug('Setup: Preparing inputs for mini-apertures')

            if isinstance(self.mini, float) or isinstance(self.mini, int):

                mini = float(self.mini)
                adj_semimaj_pix_full = self.adj_semimaj_pix
                self.adj_semimaj_pix = mini

            else: pdb.set_trace()

        else: adj_semimaj_pix_full = self.adj_semimaj_pix

        adj_semimin_pix_full = adj_semimaj_pix_full / self.adj_axial_ratio

        # Define charactaristics of circular aperture with same area as elliptical source aperture
        ap_area = np.pi * self.adj_semimaj_pix * (self.adj_semimaj_pix / self.adj_axial_ratio)
        sky_ap_rad_pix = ( ap_area / np.pi )**0.5
        sky_border = int( sky_ap_rad_pix + 1.0 ) #int( ( band_dict['annulus_outer'] * sky_ap_rad_pix ) + 1 )
        adj_semimin_pix = self.adj_semimaj_pix / self.adj_axial_ratio

        # Creating mask maps to describe no-go regions
        if ap_debug: log.debug('Setup: Creating mask maps')

        prior_mask = np.zeros(self.cutout.shape)
        exclude_mask = ChrisFuncs.Photom.EllipseMask(self.cutout, adj_semimaj_pix_full, self.adj_axial_ratio, self.adj_angle, self.centre_i, self.centre_j)
        flag_mask = np.zeros(self.cutout.shape)
        attempt_mask = np.zeros(self.cutout.shape)

        # Set pixels in source aperture to all have NaN pixels, so they don't get sampled by sky annuli
        cutout_inviolate = self.cutout.copy()
        self.cutout[ np.where(ChrisFuncs.Photom.EllipseMask(self.cutout, adj_semimaj_pix_full, self.adj_axial_ratio, self.adj_angle, self.centre_i, self.centre_j) == 1) ] = np.NaN

        # Define how many random aperture are desired/required/permitted
        sky_success_target = 50
        sky_success_min = 20
        sky_gen_max = 100

        # Generate random polar coordinates to draw from
        if ap_debug: log.debug('Setup: Generating pool of random polar coordinates')

        random_size = sky_success_target * sky_gen_max * 10
        random_failed = []
        random_theta_list = 360.0 * np.random.rand(random_size)
        random_r_list = adj_semimin_pix + np.abs(np.random.normal(loc=0.0, scale=5.0*adj_semimaj_pix_full, size=random_size))

        # Locate contiguous map regions
        if ap_debug: log.debug('Pruning: Locating contiguous coverage regions')

        cont_binary = np.zeros(self.cutout.shape)
        cont_binary[ np.where( np.isnan(cutout_inviolate)==False ) ] = 1
        cont_structure = np.array([[1,1,1], [1,1,1], [1,1,1]])
        cont_label = scipy.ndimage.measurements.label(cont_binary, structure=cont_structure)[0]
        cont_search_mask = ChrisFuncs.EllipseMask(cont_label, 3.0, 1.0, 0.0, self.centre_i, self.centre_j)
        cont_search_values = cont_label[ np.where( cont_search_mask==1 ) ]

        # Identify contiguous region associated with target source
        if ap_debug: log.debug('Pruning: Identifying coverage region associated with source')

        if np.where(cont_search_values>0)[0].shape[0] == 0: cont_target = 0
        else: cont_target = scipy.stats.mode(cont_search_values[np.where(cont_search_values>0)])[0][0]
        cont_where_bad = np.where(cont_label != cont_target)

        # Remove random coordinates that are more distant than most distant part of the coverage region the target source lies in
        if ap_debug: log.debug('Pruning: Removing random coords that definetily lie outside coverage region')

        cont_size_i, cont_size_j = self.cutout.shape
        cont_range_i = np.arange(cont_size_i) - self.centre_i
        cont_range_j = np.arange(cont_size_j) - self.centre_j
        cont_coord_i, cont_coord_j = np.meshgrid(cont_range_j, cont_range_i)  # Yes, i and j are supposed to be this way around inside meshgrid (for some reason)
        cont_coord_i[cont_where_bad] = np.NaN
        cont_coord_j[cont_where_bad] = np.NaN
        cont_dist = np.sqrt(cont_coord_i**2 + cont_coord_j**2)
        cont_dist_max = np.nanmax(cont_dist)
        random_r_coverage = np.where( random_r_list < (cont_dist_max-sky_border) )
        random_r_list = random_r_list[random_r_coverage]
        random_theta_list = random_theta_list[random_r_coverage]

        # Convert random polar coordinates into cartesian coordinates
        if ap_debug: log.debug('Pruning: Converting random polar coords to cartesian coords, and removing those that lie beyond map border')
        random_i_list = self.centre_i + ( random_r_list * np.cos(np.radians(random_theta_list)) )#np.random.normal(loc=centre_i, scale=2.0*sky_ap_rad_pix)
        random_j_list = self.centre_j + ( random_r_list * np.sin(np.radians(random_theta_list)) )#np.random.normal(loc=centre_j, scale=2.0*sky_ap_rad_pix)

        # Remove random coodinates that fall fall beyond border in i-coords
        random_not_i_border = np.where( (random_i_list>sky_border) & (random_i_list < (self.cutout.shape[0]-sky_border)) )
        random_i_list = random_i_list[random_not_i_border]
        random_j_list = random_j_list[random_not_i_border]

        # Remove random coodinates that fall fall beyond border in j-coords
        random_not_j_border = np.where( (random_j_list>sky_border) & (random_j_list < (self.cutout.shape[1]-sky_border)) )
        random_i_list = random_i_list[random_not_j_border]
        random_j_list = random_j_list[random_not_j_border]

        # Remove random coordinates that intersect source
        if ap_debug: log.debug('Pruning: Removing random coords that intersect source')
        random_not_source = np.where( np.sqrt( (np.abs(self.centre_i-random_i_list))**2.0 + (np.abs(self.centre_j-random_j_list))**2.0 ) > adj_semimin_pix_full )
        #random_not_source = np.where( (abs(centre_i-random_i_list)>adj_semimin_pix_full) & (abs(centre_j-random_j_list)>adj_semimin_pix_full) )
        random_i_list = random_i_list[random_not_source]
        random_j_list = random_j_list[random_not_source]

        # Remove random coordinates that correspond to NaN pixels
        if ap_debug: log.debug('Pruning: Removing random coords that correspond to NaN pixels')
        random_i_list_pix = np.round(random_i_list).astype(int)
        random_j_list_pix = np.round(random_j_list).astype(int)
        random_ij_values = self.cutout[(random_i_list_pix, random_j_list_pix)]
        random_ij_pix_good = np.where(np.isnan(random_ij_values)==False)
        random_i_list = random_i_list[random_ij_pix_good]
        random_j_list = random_j_list[random_ij_pix_good]

        # If none of the apertures are suitable, immediately report failure
        if random_i_list.shape[0] == 0:

            if ap_debug: log.debug('Status: Pruning removed all generated coordinates')
            #ap_noise_dict = {'fail':True, 'prior_mask':prior_mask, 'flag_mask':flag_mask, 'sky_success_counter':0}

            self.success = False
            self.prior_mask = prior_mask
            self.flag_mask = flag_mask
            self.sky_success_counter = 0

            self.cutout_inviolate = cutout_inviolate

            return

        # Commence creation of random sky apertures
        sky_success_counter = 0
        sky_sum_list = []
        sky_total_fail = False
        while True:

            # Repeatedly generate random sky apertures, until an acceptable aperture is generated
            sky_gen_counter = 0
            sky_gen_fail = False
            while True:
                sky_gen_counter += 1

                # If more than a given number of unsuccessful sky apertures have been generated in a row, call it a day
                if sky_gen_counter>sky_gen_max:

                    sky_gen_fail = True
                    if ap_debug: log.debug('Status: Unable to generate suitable random sky aperture after ' + str(sky_gen_max) + ' attempts')
                    break

                # Select random coordinate set for this iteration; if no un-used random coordinates can be found, reject
                #if ap_debug: print 'Checking: Selecting random coordinates'
                random_accept = False
                random_reject_count = 0
                while random_accept == False:

                    random_index = int(np.floor( np.random.rand() * float(random_i_list.shape[0])))

                    if random_index not in random_failed: random_accept = True
                    else: random_reject_count += 1
                    if random_reject_count > (10 * random_i_list.shape[0]): break

                if random_reject_count > (10 * random_i_list.shape[0]):

                    if ap_debug: log.debug('Rejection: Unable to find un-used random coodinates')
                    continue

                random_i = random_i_list[random_index]
                random_j = random_j_list[random_index]

                if ap_debug:

                    ap_mask = ChrisFuncs.Photom.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                    attempt_mask[ np.where(ap_mask==1) ] = sky_success_counter
                    log.debug('Aperture: '+str(sky_success_counter + 1) + ';   Generation: '+str(sky_gen_counter) + ';   Pix Coords: [' + str(random_i) + ',' + str(random_j)+']')

                # Do sophisticated check that generated sky aperture does not intersect source; if it does, reject
                #if ap_debug: print 'Checking: Determining whether aperture intersects source (sophisticated check)'
                exclude_sum = ChrisFuncs.Photom.EllipseSum(exclude_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)[0]
                if exclude_sum > 0:

                    if ap_debug: log.debug('Rejection: Aperture intersects source (according to sophisticated check)')
                    random_failed.append(random_index)
                    continue

                # Do basic chrck that the majority of the pixels in the generated sky aperture have not already been sampled by previous sky apertures; they have, reject
                #if ap_debug: print 'Checking: Determining if aperture over-sampled (basic check)'
                prior_calc = ChrisFuncs.Photom.EllipseSum(prior_mask, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                prior_calc[2][ np.where(prior_calc[2]>=1.0) ] = 1.0
                prior_frac = np.sum(prior_calc[2]) / float(prior_calc[1])
                if prior_frac>0.5:
                    if ap_debug: log.debug('Rejection: Aperture over-sampled (according to basic check)')
                    random_failed.append(random_index)
                    continue

                # Do sophisticated check that the majority of the pixels in the generated sky aperture have not already been sampled by previous sky apertures; they have, reject
                #if ap_debug: print 'Checking: Determinging if aperture oversampled (sophisticated check)'
                ap_mask_check = ChrisFuncs.Photom.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
                flag_mask_check = flag_mask.copy()
                flag_mask_check[np.where(ap_mask_check==1)] = int(2.0**(sky_success_counter+1.0))
                flag_tallies = np.array([ np.where(flag_mask_check==flag)[0].shape[0] for flag in (2.0**np.arange(0.0,sky_success_counter+2.0)).tolist() ])
                flag_check = np.where(flag_tallies<(0.5*ap_area))[0].shape[0]
                if flag_check > 1:

                    if ap_debug: log.debug('Rejection: Aperture over-sampled (according to sophisticated check)')
                    random_failed.append(random_index)
                    continue

                # Evaluate pixels in sky aperture
                #if ap_debug: print 'Checking: Evaluating pixels in sky aperture'
                ap_calc = ChrisFuncs.Photom.EllipseSum(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)

                # Evaluate pixels in sky annulus
                #if ap_debug: print 'Checking: Evaluating pixels in sky anmnulus'
                bg_inner_semimaj_pix = self.adj_semimaj_pix * band_dict['annulus_inner']
                bg_width = (self.adj_semimaj_pix * band_dict['annulus_outer']) - bg_inner_semimaj_pix
                bg_calc = ChrisFuncs.Photom.AnnulusSum(self.cutout, bg_inner_semimaj_pix, bg_width, 1.0, 0.0, random_i, random_j)

                # Check if more than a givn fraction of the pixels inside the source aperture are NaN; if so, reject
                if ap_calc[3][0].shape[0] == 0 or ap_calc[1] == 0: ap_nan_frac = 0.0
                if ap_calc[1] == 0: ap_nan_frac = 1.0
                else: ap_nan_frac = float(ap_calc[3][0].shape[0]) / float(ap_calc[1]+float(ap_calc[3][0].shape[0]))
                ap_nan_thresh = 0.10

                if ap_nan_frac > ap_nan_thresh:

                    if ap_debug: log.debug('Rejection: Aperture contains too many NaNs')
                    random_failed.append(random_index)
                    continue

                # Check if more than a given fraction of the pixels inside the sky annulus are NaN; if so, reject
                if bg_calc[3][0].shape[0] == 0: bg_nan_frac = 0.0
                if bg_calc[1] == 0: bg_nan_frac = 1.0
                else: bg_nan_frac = float(bg_calc[3][0].shape[0]) / float(bg_calc[1]+bg_calc[3][0].shape[0])

                bg_nan_thresh = 0.80
                if bg_nan_frac > bg_nan_thresh:

                    if ap_debug: log.debug('Rejection: Annulus contains too many NaNs')
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

                    if ap_debug: log.debug('Status: However, sufficient number of successful random apertures ('+str(int(sky_success_counter))+') already generated; proceeding')
                    break

            # Calculate actual flux in sky aperture, and record
            #if ap_debug: print 'Checking: Performing photometry with random sky aperture and annulus'
            bg_clip = ChrisFuncs.SigmaClip(bg_calc[2], median=False, sigma_thresh=3.0)
            bg_avg = bg_clip[1]
            ap_sum = ap_calc[0] - (ap_calc[1] * bg_avg)
            sky_sum_list.append(ap_sum)
            if np.isnan(ap_sum):
                pdb.set_trace()

            # Add this aperture to the prior mask and flag mask
            if not ap_debug: ap_mask = ChrisFuncs.Photom.EllipseMask(self.cutout, sky_ap_rad_pix, 1.0, 0.0, random_i, random_j)
            prior_mask += ap_mask
            flag_mask[np.where(ap_mask==1)] += 2.0**sky_success_counter

            # If target number of sky apertures have been processed, break out of loop
            if sky_success_counter>=sky_success_target:

                if ap_debug: log.debug('Status: Target number of successful random apertures (' + str(int(sky_success_target)) + ') generated; proceeding')
                break

        # If total failure was encountered, end process and report now
        if sky_total_fail:

            #ap_noise_dict = {'fail': True, 'prior_mask': prior_mask, 'flag_mask': flag_mask, 'sky_success_counter':sky_success_counter}

            #return ap_noise_dict

            self.success = False
            self.prior_mask = prior_mask
            self.flag_mask = flag_mask
            self.sky_success_counter = sky_success_counter

            self.cutout_inviolate = cutout_inviolate

            return

        # Otherwise, calculate aperture noise using returned aperture values, and return
        else:

            sky_sum_list = np.array(sky_sum_list)
            ap_noise = ChrisFuncs.SigmaClip(sky_sum_list, tolerance=0.001, median=True, sigma_thresh=3.0)[0]

            #ap_noise_dict = {'fail': False, 'ap_noise':ap_noise, 'ap_num':sky_success_counter, 'prior_mask':prior_mask, 'flag_mask':flag_mask}
            #ChrisFuncs.Cutout(prior_mask, '/home/saruman/spx7cjc/DustPedia/Prior.fits')

            # Debugging
            log.debug('Aperture noise from current random apertures is '+str(ChrisFuncs.FromGitHub.randlet.ToPrecision(ap_noise,4))+' (in map units).')

            self.success = True
            self.noise = ap_noise
            self.napertures = sky_success_counter
            self.prior_mask = prior_mask
            self.flag_mask = flag_mask

            self.cutout_inviolate = cutout_inviolate

            return

# -----------------------------------------------------------------

class ExtrapolatingApertureNoiseCalculator(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The consturctor ...
        """

        # INPUT

        self.cutout = None
        self.source_dict = None
        self.band_dict = None
        self.kwargs_dict = None
        self.adj_semimaj_pix = None
        self.adj_axial_ratio = None
        self.adj_angle = None
        self.centre_i = None
        self.centre_j = None

        # OUTPUT

        # The aperture noise
        self.noise = None

        # Success flag
        self.success = False

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

        self.cutout = input_dict.pop("cutout")
        self.source_dict = input_dict.pop("source_dict")
        self.band_dict = input_dict.pop("band_dict")
        self.kwargs_dict = input_dict.pop("kwargs_dict")
        self.adj_semimaj_pix = input_dict.pop("adj_semimaj_pix")
        self.adj_axial_ratio = input_dict.pop("adj_axial_ratio")
        self.adj_angle = input_dict.pop("adj_angle")
        self.centre_i = input_dict.pop("centre_i")
        self.centre_j = input_dict.pop("centre_j")

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
        mini_ap_rad_base = 2.0
        mini_ap_rad_pix_input = mini_ap_rad_base**np.arange(1.0, np.ceil( math.log( sky_ap_rad_pix, mini_ap_rad_base)))[::-1]
        min_ap_rad_pix_output = []
        mini_ap_noise_output = []
        mini_ap_num_output = []

        # Loop over radii for mini-apertures
        for mini_ap_rad_pix in mini_ap_rad_pix_input:

            log.debug('Finding aperture noise for mini-apertures of radius '+str(mini_ap_rad_pix)[:6] + ' pixels.')

            mini_calculator = ExactApertureNoiseCalculator()


            input_dict_radius = dict()
            input_dict_radius["cutout"] = self.cutout.copy()
            input_dict_radius["source_dict"] = self.source_dict
            input_dict_radius["band_dict"] = self.band_dict
            input_dict_radius["kwargs_dict"] = self.kwargs_dict
            input_dict_radius["adj_semimaj_pix"] = self.adj_semimaj_pix
            input_dict_radius["adj_axial_ratio"] = self.adj_axial_ratio
            input_dict_radius["adj_angle"] = self.adj_angle
            input_dict_radius["centre_i"] = self.centre_i
            input_dict_radius["centre_j"] = self.centre_j
            input_dict_radius["mini_ap_rad_pix"] = mini_ap_rad_pix
            input_dict_radius["downsample_factor"] = self.downsample_factor

            mini_calculator.run(**input_dict_radius)

            #mini_ap_noise_dict = ApNoise(cutout.copy(), self.source_dict, self.band_dict, self.kwargs_dict, self.adj_semimaj_pix, self.adj_axial_ratio, self.adj_angle, self.centre_i, self.centre_j, mini=mini_ap_rad_pix, downsample=int(band_dict['downsample_factor']))

            # If mini-aperture succeeded, record and proceed; else, call it a day
            if mini_calculator.success:

                min_ap_rad_pix_output.append(mini_ap_rad_pix)
                mini_ap_noise_output.append(mini_calculator.noise)
                mini_ap_num_output.append(mini_calculator.napertures)

            # No succes for this radius
            else: log.debug('Unable to place sufficient number of mini-apertures at this radius.')

            # Stop when we have 6 datapoints
            if len(mini_ap_noise_output) >= 6: break

        # Convert output lists into arrays
        min_ap_rad_pix_output = np.array(min_ap_rad_pix_output)
        mini_ap_noise_output = np.array(mini_ap_noise_output)
        mini_ap_num_output = np.array(mini_ap_num_output).astype(float)

        # If insufficient points to make extrapolation, report failure; else proceed
        if min_ap_rad_pix_output.shape[0] < 2:

            #ap_noise_dict = {'fail':True, 'ap_noise':np.NaN}

            self.success = False
            self.noise = None

            gc.collect()
            #return ap_noise_dict

            return

        else:

            # Calculate values to plot
            log_mini_ap_area = np.log10(np.pi*min_ap_rad_pix_output**2.0)
            log_mini_ap_noise = np.log10(mini_ap_noise_output)
            mini_ap_noise_err_rel = mini_ap_num_output**0.5 / mini_ap_num_output
            #mini_ap_noise_err = np.abs( mini_ap_noise_output * mini_ap_noise_err_rel )
            log_mini_ap_noise_err = mini_ap_noise_err_rel#ChrisFuncs.LogError(mini_ap_noise_output, mini_ap_noise_err)

            # Define straight-line function, and fit it to points
            def Line(x,m,c):
                return (m*x)+c

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
            for xlabel in ax.get_xticklabels():
                xlabel.set_fontproperties(font_manager.FontProperties(size=15))
            for ylabel in ax.get_yticklabels():
                ylabel.set_fontproperties(font_manager.FontProperties(size=15))

            # Plot axis 1 legend
            ax_handles, ax_labels = ax.get_legend_handles_labels()
            ax1_lgd = ax.legend(ax_handles, ax_labels, loc='best', scatterpoints=1, labelspacing=0.25, borderpad=0)
            ax1_lgd.draw_frame(False)
            plt.setp(plt.gca().get_legend().get_texts(), fontsize='12.5')

            # Save figure, clean up, and report results
            fig.savefig(fs.join( kwargs_dict['temp_dir_path'], "Aperture_Noise_Projection.png"), dpi=100)
            gc.collect()
            fig.clear()
            plt.close('all')

            #ap_noise_dict = {'fail': False, 'ap_noise': ap_noise_proj}
            #return ap_noise_dict

            self.success = True
            self.noise = ap_noise_proj

            return

# -----------------------------------------------------------------

def CircularApertureStandardDeviationFinder(fits, area, ann=True, ann_inner=1.5, ann_width=1.0, angle=0.0, axial_ratio=1.0, apertures=100, clipping=True, mask=False, maskap=[1,1]):

    """
    Got this function from Pieter De Vis
    :param fits:
    :param area:
    :param ann:
    :param ann_inner:
    :param ann_width:
    :param angle:
    :param axial_ratio:
    :param apertures:
    :param clipping:
    :param mask:
    :param maskap:
    :return:
    """

    # Calculate aperture's radius from given circular area
    rad = np.sqrt(area/np.pi)

    # Define width of border
    if ann==True:
        factor = ann_inner + ann_width
    elif ann==False:
        factor = 1.0

    # Define exclusion zone
    semi_maj = ( area / (axial_ratio*np.pi) )**0.5
    source_i_centre, source_j_centre = int(round(float(fits.shape[0])/2.0)), int(round(float(fits.shape[1])/2.0))
    exclusion_mask = ChrisFuncs.EllipseMask(fits, 1.2*semi_maj, axial_ratio, angle, source_i_centre, source_j_centre)

    # Definine initial values and initiate loop
    if rad < 2:
        rad = 2

    # Define limits for fixed sigma-clip
    fits_clip = ChrisFuncs.SigmaClip(fits[ where( (exclusion_mask==0)) ], sigma_thresh=2.0, median=True, no_zeros=True)
    clip_upper = fits_clip[1] + (3.0 * fits_clip[0])
    clip_lower = fits_clip[1] - (3.0 * fits_clip[0])

    # Loop over apertures
    ap_sum_array = []
    rejects_array = []
    hist=[]
    a = 0
    b = 0
    sum_ratio_array = []
    while a<apertures:
        print a
        # Generate box coordinates, excluding cutout border, and target source
        ap_border = int((factor*rad)+1)
        neverending = 1
        while neverending > 0:
            neverending = 1
            i_centre = random.randint(ap_border, fits.shape[0]-ap_border)
            j_centre = random.randint(ap_border, fits.shape[1]-ap_border)

            # Do crude check that generated coords do not intersect source
            if ( abs(0.5*fits.shape[0] - i_centre) < rad ) and ( abs(0.5*fits.shape[1] - j_centre) < rad ):
                continue

            # Do sophisticated check that generated coords do not intersect source
            exclusion_cutout = exclusion_mask[ i_centre-(rad+1):i_centre+(rad+1) , j_centre-(rad+1):j_centre+(rad+1) ]
            exclusion_sum = np.sum(exclusion_cutout)
            if exclusion_sum>0:
                continue
            else:
                break

        # Extract slice around aperture; if slice has either dimension 0, then continue
        cutout = fits[ i_centre-(rad*factor+1):i_centre+(rad*factor+1) , j_centre-(rad*factor+1):j_centre+(rad*factor+1) ]
        box_centre = int(round(cutout.shape[0]/2.0))
        if cutout.shape[0]==0 or cutout.shape[1]==0:
            continue

        # If background-subtraction required, measure sky annulus
        if ann==True:
            bg_phot = ChrisFuncs.AnnulusSum(cutout, ann_inner*rad, ann_width*rad, 1.0, 0.0, box_centre, box_centre)
        bg_mean = np.mean(bg_phot[2])
        if 0 in bg_phot[2]:
            continue
            if clipping==True:
            bg_clipped = bg_phot[2][ np.where( (bg_phot[2]<clip_upper) & (bg_phot[2]>clip_lower) ) ]
            bg_clipped = bg_clipped[ np.where( bg_clipped<clip_upper ) ]
            bg_clipped = bg_clipped[ np.where( bg_clipped>clip_lower ) ]
            bg_mean=np.mean(bg_clipped)
        elif ann==False:
            bg_mean = 0.0

        # Measure sky aperture, and reject if it contains nothing but NaN pixels
        ap_phot = ChrisFuncs.EllipseSum(cutout, rad, 1.0, 0.0, box_centre, box_centre)
    if mask:
        cutoutmask = maskap[ i_centre-(rad*factor+1):i_centre+(rad*factor+1) , j_centre-(rad*factor+1):j_centre+(rad*factor+1) ]
        ap_mask = ChrisFuncs.EllipseSum(cutoutmask, rad, 1.0, 0.0, box_centre, box_centre)
		if 0 in ap_mask[2]:
			print "Steve mask"
			continue
        ap_sum  = np.sum(ap_phot[2]) - bg_mean*ap_phot[1]
        if ap_phot[2].shape[0] == 0:
            continue
        if np.where( ap_phot[2]==0 )[0].shape[0] == ap_phot[2].shape[0]:
            continue
    if (0 in ap_phot[2]) and mask: continue
    b += 1
	if b > 1000:
        panic
        # Perform fixed sigma-clip of sky aperture; reject if >20% of pixels are clipped, otherwise scale flux for removed pixels
	if clipping:
        rejects_array.append(ap_phot[0])
        ap_clipped = ap_phot[2][ np.where( (ap_phot[2]<clip_upper) & (ap_phot[2]>clip_lower) ) ]
        ap_clipped = ap_clipped[ np.where( ap_clipped<clip_upper ) ]
		ap_clipped = ap_clipped[ np.where( ap_clipped>clip_lower ) ]
		clip_frac = float(ap_clipped.shape[0]) / float(ap_phot[2].shape[0])
		if clip_frac<0.8: continue
		ap_sum = np.sum(ap_clipped) / clip_frac - bg_mean*ap_phot[1]
		sum_ratio_array.append(ap_sum/ap_phot[0])

		# Store calculated sum and tick up counter
		ap_clipped = ap_phot[2][ np.where( ap_phot[2]<clip_upper ) ]
		clip_frac = float(ap_clipped.shape[0]) / float(ap_phot[2].shape[0])
		hist.append(clip_frac)
	a += 1
	ap_sum_array.append(ap_sum)	

    #bina=pl.hist(ap_sum_array)
    #pl.hist(ap_sum_array2,bins=bina[1])
    #pl.savefig("histograms/histap_"+str(time.time())+".png")
    #pl.clf()
    
    # Now take standard deviation of output array to get a noise value as pure as the driven snow
    ap_sum_array = np.array(ap_sum_array)
    ap_sum_array = ap_sum_array[ np.where( np.isnan(ap_sum_array)==False ) ]
    ap_sum_array = ap_sum_array[ np.where( ap_sum_array!=0.0 ) ]
    ap_sigma = np.std(ap_sum_array)
    
    rejects_array = np.array(rejects_array)
    rejects_array = rejects_array[0:apertures]
    ap_mean=np.mean(np.array(ap_sum_array2)/float(ap_phot[1]))

    # Report calculated standard deviation
    return ap_sigma, ap_sum_array, ap_mean
    
# -----------------------------------------------------------------
