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
import numpy as np
import lmfit
import time

# Import astronomical modules
import astropy.convolution
import astropy.io.fits
import astropy.wcs
import astropy.modeling

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ..misc import chrisfuncs

# -----------------------------------------------------------------

class ApertureCorrector(Configurable):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(ApertureCorrector, self).__init__(*args, **kwargs)

        # The input
        self.psf = None
        self.cutout = None

        # The aperture correction factor
        self.factor = None

    # -----------------------------------------------------------------

    def run(self, **input_dict):

        """
        This function ...
        :param input_dict:
        :return:
        """

        # 1. Call the setup function
        self.setup(**input_dict)

        # 2. Calculate the correction factor
        self.calculate_aperture_correction()

    # -----------------------------------------------------------------

    def setup(self, **input_dict):

        """
        This function ...
        :param input_dict:
        :return:
        """

        # Call the setup function of the base class
        super(ApertureCorrector, self).setup()

        # Set the input
        self.psf = input_dict.pop("psf")
        self.cutout = input_dict.pop("cutout")

    # -----------------------------------------------------------------

    def calculate_aperture_correction(self):

        """
        # Define function that uses provided beam profile to aperture-correct photometry

        INPUT_DICT:
        ----------

        psf:    the PSF, or a False boolean
        cutout: Array upon which photometry is being perfomred upon.

        CONFIG:
        -------

        pix_arcsec: The width, in arscec, of the pixels in the map photometry is being performed upon
                    (this is needed in case there is a pixel size mismatch with PSF).
        semimaj_pix: Semi-major axis of photometric aperture, in pixels.
        axial_ratio: Axial ratio of photometryic aperture.
        angle: Position angle of photometric aperture, in degrees.

        centre_i: Zero-indexed, 0th-axis coordinate (equivalent to y-axis one-indexed coordinates in FITS terms)
                  of centre position of photometric aperture.
        centre_j: Zero-indexed, 1st-axis coordinate (equivalent to x-axis one-indexed coordinates in FITS terms) of
                  centre position of photometric aperture.

        ... and other ... (see configuration definition)

        """

        ### INPUT

        psf = self.psf.data  # PSF MUST BE PREPARED

        # Cutout
        cutout = self.cutout.data

        #####

        # Produce mask for pixels we care about for fitting (ie, are inside photometric aperture and background annulus)
        mask = chrisfuncs.EllipseMask(cutout, self.config.semimaj_pix, self.config.axial_ratio,
                                             self.config.angle, self.config.centre_i,
                                             self.config.centre_j)  # *band_dict['annulus_outer']

        ##
        # from pts.magic.tools import plotting
        # plotting.plot_mask(mask)
        ##

        # Produce guess values
        initial_sersic_amplitude = cutout[int(round(self.config.centre_i)), int(round(self.config.centre_j))]
        initial_sersic_r_eff = self.config.semimaj_pix / 10.0
        initial_sersic_n = 1.0
        initial_sersic_x_0 = self.config.centre_j
        initial_sersic_y_0 = self.config.centre_i
        initial_sersic_ellip = (self.config.axial_ratio - 1.0) / self.config.axial_ratio
        initial_sersic_theta = np.deg2rad(self.config.angle)

        # Produce sersic model from guess parameters, for time trials
        sersic_x, sersic_y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))
        sersic_model = astropy.modeling.models.Sersic2D(amplitude=initial_sersic_amplitude, r_eff=initial_sersic_r_eff,
                                                        n=initial_sersic_n, x_0=initial_sersic_x_0,
                                                        y_0=initial_sersic_y_0, ellip=initial_sersic_ellip,
                                                        theta=initial_sersic_theta)
        sersic_map = sersic_model(sersic_x, sersic_y)

        # Make sure that PSF array is smaller than sersic model array (as required for convolution); if not, remove its edges such that it is
        if psf.shape[0] > sersic_map.shape[0] or psf.shape[1] > sersic_map.shape[1]:
            excess = max(psf.shape[0] - sersic_map.shape[0], psf.shape[1] - sersic_map.shape[1])
            border = max(2, int(np.round(np.ceil(float(excess) / 2.0) - 1.0)))
            psf = psf[border:, border:]
            psf = psf[:-border, :-border]

        # Determine whether FFT convolution or direct convolution is faster for this kernel,
        # using sersic model produced with guess parameters
        time_fft = time.time()
        conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
        time_fft = time.time() - time_fft
        time_direct = time.time()
        conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)
        time_direct = time.time() - time_direct

        if time_fft < time_direct: use_fft = True
        else: use_fft = False

        # Set up parameters to fit galaxy with 2-dimensional sersic profile
        params = lmfit.Parameters()
        params.add('sersic_amplitide', value=initial_sersic_amplitude, vary=True)
        params.add('sersic_r_eff', value=initial_sersic_r_eff, vary=True, min=0.0, max=self.config.semimaj_pix)
        params.add('sersic_n', value=initial_sersic_n, vary=True, min=0.1, max=10)
        params.add('sersic_x_0', value=initial_sersic_x_0, vary=False)
        params.add('sersic_y_0', value=initial_sersic_y_0, vary=False)
        params.add('sersic_ellip', value=initial_sersic_ellip, vary=True, min=0.5 * initial_sersic_ellip, max=0.5 * (1.0 - initial_sersic_ellip) + initial_sersic_ellip)
        params.add('sersic_theta', value=initial_sersic_theta, vary=False)

        # Solve with LMfit to find parameters of best-fit sersic profile
        result = lmfit.minimize(chi_squared_sersic, params, args=(cutout, psf, mask, use_fft), method='leastsq', ftol=1E-5, xtol=1E-5, maxfev=200)

        # Extract best-fit results
        sersic_amplitide = result.params['sersic_amplitide'].value
        sersic_r_eff = result.params['sersic_r_eff'].value
        sersic_n = result.params['sersic_n'].value
        sersic_x_0 = result.params['sersic_x_0'].value
        sersic_y_0 = result.params['sersic_y_0'].value
        sersic_ellip = result.params['sersic_ellip'].value
        sersic_theta = result.params['sersic_theta'].value

        # Construct underlying sersic map and convolved sersic map, using best-fit parameters
        sersic_model = astropy.modeling.models.Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
        sersic_map = sersic_model(sersic_x, sersic_y)

        if use_fft: conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
        else: conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)

        # Determine annulus properties before proceeding with photometry
        # bg_inner_semimaj_pix = input_dict['semimaj_pix'] * input_dict['annulus_inner'] # number of pixels of semimajor axis of inner annulus ellipse
        # bg_width = (input_dict['semimaj_pix'] * input_dict['annulus_outer']) - bg_inner_semimaj_pix # number of pixels of difference between outer major axis and minor major axis

        bg_inner_semimaj_pix = self.config.semimaj_pix_annulus_inner
        bg_width = self.config.semimaj_pix_annulus_outer - bg_inner_semimaj_pix
        axial_ratio_annulus = self.config.axial_ratio_annulus
        angle_annulus = self.config.annulus_angle
        centre_i_annulus = self.config.annulus_centre_i
        centre_j_annulus = self.config.annulus_centre_j

        # Evaluate pixels in source aperture and background annulus in UNCONVOLVED sersic map
        if self.config.subpixel_factor == 1.0:
            sersic_ap_calc = chrisfuncs.EllipseSum(sersic_map, self.config.semimaj_pix, self.config.axial_ratio, self.config.angle, self.config.centre_i, self.config.centre_j)
            sersic_bg_calc = chrisfuncs.AnnulusSum(sersic_map, bg_inner_semimaj_pix, bg_width, axial_ratio_annulus, angle_annulus, centre_i_annulus, centre_j_annulus)
        elif self.config.subpixel_factor > 1.0:
            sersic_ap_calc = chrisfuncs.EllipseSumUpscale(sersic_map, self.config.semimaj_pix, self.config.axial_ratio, self.config.angle, self.config.centre_i, self.config.centre_j, upscale=self.config.subpixel_factor)
            sersic_bg_calc = chrisfuncs.AnnulusSumUpscale(sersic_map, bg_inner_semimaj_pix, bg_width, axial_ratio_annulus, angle_annulus, centre_i_annulus, centre_j_annulus, upscale=self.config.subpixel_factor)
        else: raise ValueError("Invalid subpixel factor: " + str(self.config.subpixel_factor))

        # Background-subtract and measure UNCONVOLVED sersic source flux
        sersic_bg_clip = chrisfuncs.SigmaClip(sersic_bg_calc[2], median=False, sigma_thresh=3.0)
        sersic_bg_avg = sersic_bg_clip[1] * self.config.subpixel_factor**2.0
        sersic_ap_sum = sersic_ap_calc[0] - (sersic_ap_calc[1] * sersic_bg_avg)  # sersic_ap_calc[1] = number of pixels counted for calculating sum (total flux in ellipse)

        # Evaluate pixels in source aperture and background annulus in CONVOLVED sersic map
        if self.config.subpixel_factor == 1.0:
            conv_ap_calc = chrisfuncs.EllipseSum(conv_map, self.config.semimaj_pix, self.config.axial_ratio, self.config.angle, self.config.centre_i, self.config.centre_j)
            conv_bg_calc = chrisfuncs.AnnulusSum(conv_map, bg_inner_semimaj_pix, bg_width, axial_ratio_annulus, angle_annulus, centre_i_annulus, centre_j_annulus)
        elif self.config.subpixel_factor > 1.0:
            conv_ap_calc = chrisfuncs.EllipseSumUpscale(conv_map, self.config.semimaj_pix, self.config.axial_ratio, self.config.angle, self.config.centre_i, self.config.centre_j, upscale=self.config.subpixel_factor)
            conv_bg_calc = chrisfuncs.AnnulusSumUpscale(conv_map, bg_inner_semimaj_pix, bg_width, axial_ratio_annulus, angle_annulus, centre_i_annulus, centre_j_annulus, upscale=self.config.subpixel_factor)
        else: raise ValueError("Invalid subpixel factor: " + str(self.config.subpixel_factor))

        # Background-subtract and measure CONVOLVED sersic source flux
        conv_bg_clip = chrisfuncs.SigmaClip(conv_bg_calc[2], median=False, sigma_thresh=3.0)
        conv_bg_avg = conv_bg_clip[1] * self.config.subpixel_factor**2.0
        conv_ap_sum = conv_ap_calc[0] - (conv_ap_calc[1] * conv_bg_avg)  # conv_ap_calc[1] = number of pixels counted for calculating sum (total flux in ellipse)

        # Find difference between flux measured on convoled and unconvoled sersic maps
        ap_correction = np.nanmax([1.0, (sersic_ap_sum / conv_ap_sum)])

        # Return aperture correction
        #return ap_correction

        # Set the factor
        self.factor = ap_correction

# -----------------------------------------------------------------

def chi_squared_sersic(params, cutout, psf, mask, use_fft, lmfit=True):

    """
    This function defines LMfit convolved-sersic function
    :param params:
    :param cutout:
    :param psf:
    :param mask:
    :param use_fft:
    :param lmfit:
    :return:
    """

    # Extract variable parameters
    sersic_amplitide = params['sersic_amplitide'].value
    sersic_r_eff = params['sersic_r_eff'].value
    sersic_n = params['sersic_n'].value
    sersic_x_0 = params['sersic_x_0'].value
    sersic_y_0 = params['sersic_y_0'].value
    sersic_ellip = params['sersic_ellip'].value
    sersic_theta = params['sersic_theta'].value

    # Generate sersic model given current parameters
    sersic_x, sersic_y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n,
                                                    x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip,
                                                    theta=sersic_theta)
    sersic_map = sersic_model(sersic_x, sersic_y)

    # Convolve sersic model with PSF
    if use_fft: conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    else: conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)

    # Calculate residuals, filtered by mask
    residuals = cutout - conv_map
    mask[np.where(mask == 0.0)] = np.nan
    residuals *= mask
    residuals.flatten()
    residuals = residuals[np.where(np.isnan(residuals) == False)]

    # Return residuals
    if lmfit: return residuals ** 2.0
    else: return residuals, cutout - conv_map

# -----------------------------------------------------------------
