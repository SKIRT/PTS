#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.aperturecorrection

# -----------------------------------------------------------------

# Import standard modules
import time
import numpy as np
import scipy.ndimage
import lmfit

# Import astronomical modules
import astropy.io.fits
import astropy.wcs
import astropy.convolution
import astropy.modeling

# Other
import ChrisFuncs

# -----------------------------------------------------------------

def StandaloneApCorrect(input_dict):

    """
    Define function that uses provided beam profile to aperture-correct photometry

    Entries in input_dict:

    psf_path: Either a string giving the path to FITS file that contains the PSF, or a False boolean (in which case an airy disc PSF will be assumed).

    cutout: Array upon whcih photometry is being perfomred upon.

    pix_arcsec: The width, in arscec, of the pixels in the map photometry is being performed upon (this is needed in case there is a pixel size mismatch with PSF).

    semimaj_pix: Semi-major axis of photometric aperture, in pixels.

    axial_ratio: Axial ratio of photometryic aperture.

    angle: Position angle of photometric aperture, in degrees.

    centre_i: Zero-indexed, 0th-axis coordinate (equivalent to y-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.

    centre_j: Zero-indexed, 1st-axis coordinate (equivalent to x-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture.

    annulus_inner: The semi-major axis of the inner edge of the background annulus, in units of the semi-major axis of the source ellipse.

    annulus_outer: The semi-major axis of the outer edge of the background annulus, in units of the semi-major axis of the source ellipse.
    """

    # If no PSF given, assume Airy disc; else extract PSF from provided file
    if (str(input_dict['psf_path'])==False) or (input_dict['psf_path'] is None):
        psf = astropy.convolution.kernels.AiryDisk2DKernel(input_dict['psf_path']).array
    else:

        # Read in PSF, and establish pixel size
        psf_in, psf_header = astropy.io.fits.getdata(input_dict['psf_path'], header=True)
        psf_wcs = astropy.wcs.WCS(psf_header)
        if psf_wcs.wcs.has_cd():
            psf_cdelt = psf_wcs.wcs.cd.max()
        else:
            psf_cdelt = psf_wcs.wcs.cdelt.max()
        psf_cdelt_arcsec = abs( psf_cdelt * 3600.0 )

        # If PSF pixel size is different to map pixel size, rescale PSF accordingly
        if (input_dict['pix_arcsec']/psf_cdelt_arcsec)>1.001 or (input_dict['pix_arcsec']/psf_cdelt_arcsec)<0.999:

            # 'Trim' the edges of the input PSF until the rescaled PSF has odd dimensions
            psf_even = True
            while psf_even:
                zoom_factor = float(psf_cdelt_arcsec) / float(input_dict['pix_arcsec'])
                psf = scipy.ndimage.zoom(psf_in, (zoom_factor,zoom_factor), mode='nearest')
                if (psf.shape[0]%2!=0) and (psf.shape[1]%2!=0):
                    psf_even = False
                else:
                    psf_in = psf_in[1:,1:]
                    psf_in = psf_in[:-1,:-1]

        # Else, if pixel sizes are already the same, leave as-is
        elif ((input_dict['pix_arcsec']/psf_cdelt_arcsec)>=0.999) and ((input_dict['pix_arcsec']/psf_cdelt_arcsec)<=0.001):
            psf = psf_in.copy()

    # Normalise PSF
    psf /= np.nansum(psf)

    # Extract cutout
    cutout = input_dict['cutout']

    # Produce mask for pixels we care about for fitting (ie, are inside photometric aperture and background annulus)
    mask = ChrisFuncs.Photom.EllipseMask(cutout, input_dict['semimaj_pix'], input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j']) #*band_dict['annulus_outer']

    # Produce guess values
    initial_sersic_amplitide = cutout[ input_dict['centre_i'], input_dict['centre_j'] ]
    initial_sersic_r_eff = input_dict['semimaj_pix'] / 10.0
    initial_sersic_n = 1.0
    initial_sersic_x_0 = input_dict['centre_j']
    initial_sersic_y_0 = input_dict['centre_i']
    initial_sersic_ellip = ( input_dict['axial_ratio'] - 1.0 ) / input_dict['axial_ratio']
    initial_sersic_theta = np.deg2rad( input_dict['angle'] )

    # Produce sersic model from guess parameters, for time trials
    sersic_x, sersic_y = np.meshgrid( np.arange(cutout.shape[1]), np.arange(cutout.shape[0]) )
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=initial_sersic_amplitide, r_eff=initial_sersic_r_eff, n=initial_sersic_n, x_0=initial_sersic_x_0, y_0=initial_sersic_y_0, ellip=initial_sersic_ellip, theta=initial_sersic_theta)
    sersic_map = sersic_model(sersic_x,sersic_y)

    # Make sure that PSF array is smaller than sersic model array (as required for convolution); if not, remove its edges such that it is
    if psf.shape[0]>sersic_map.shape[0] or psf.shape[1]>sersic_map.shape[1]:
        excess = max( psf.shape[0]-sersic_map.shape[0], psf.shape[1]-sersic_map.shape[1] )
        border = int( np.round( np.ceil( float(excess) / 2.0 ) - 1.0 ) )
        psf = psf[border:,border:]
        psf = psf[:-border,:-border]

    # Determine wither FFT convolution or direct convolution is faster for this kernel, using sersic model produced with guess parameters
    time_fft = time.time()
    conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    time_fft = time.time() - time_fft
    time_direct = time.time()
    conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)
    time_direct = time.time() - time_direct
    if time_fft<time_direct:
        use_fft = True
    else:
        use_fft = False

    # Set up parameters to fit galaxy with 2-dimensional sersic profile
    params = lmfit.Parameters()
    params.add('sersic_amplitide', value=initial_sersic_amplitide, vary=True)
    params.add('sersic_r_eff', value=initial_sersic_r_eff, vary=True, min=0.0, max=input_dict['semimaj_pix'])
    params.add('sersic_n', value=initial_sersic_n, vary=True, min=0.1, max=10)
    params.add('sersic_x_0', value=initial_sersic_x_0, vary=False)
    params.add('sersic_y_0', value=initial_sersic_y_0, vary=False)
    params.add('sersic_ellip', value=initial_sersic_ellip, vary=True, min=0.5*initial_sersic_ellip, max=0.5*(1.0-initial_sersic_ellip)+initial_sersic_ellip)
    params.add('sersic_theta', value=initial_sersic_theta, vary=False)

    # Solve with LMfit to find parameters of best-fit sersic profile
    result = lmfit.minimize(Sersic_LMfit, params, args=(cutout, psf, mask, use_fft), method='leastsq', ftol=1E-5, xtol=1E-5)

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
    if use_fft==True:
        conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    elif use_fft==False:
        conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)

    # Determine annulus properties before proceeding with photometry
    bg_inner_semimaj_pix = input_dict['semimaj_pix'] * input_dict['annulus_inner']
    bg_width = (input_dict['semimaj_pix'] * input_dict['annulus_outer']) - bg_inner_semimaj_pix

    # Evaluate pixels in source aperture and background annulus  unconvoled sersic map
    sersic_ap_calc = ChrisFuncs.Photom.EllipseSum(sersic_map, input_dict['semimaj_pix'], input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])
    sersic_bg_calc = ChrisFuncs.Photom.AnnulusSum(sersic_map, bg_inner_semimaj_pix, bg_width, input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])

    # Background-subtract and measure unconvoled sersic source flux
    sersic_bg_clip = ChrisFuncs.SigmaClip(sersic_bg_calc[2], median=False, sigma_thresh=3.0)
    sersic_bg_avg = sersic_bg_clip[1]
    sersic_ap_sum = sersic_ap_calc[0] - (sersic_ap_calc[1] * sersic_bg_avg)

    # Evaluate pixels in source aperture and background annulus in convolved sersic map
    conv_ap_calc = ChrisFuncs.Photom.EllipseSum(conv_map, input_dict['semimaj_pix'], input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])
    conv_bg_calc = ChrisFuncs.Photom.AnnulusSum(conv_map, bg_inner_semimaj_pix, bg_width, input_dict['axial_ratio'], input_dict['angle'], input_dict['centre_i'], input_dict['centre_j'])

    # Background-subtract and measure convolved sersic source flux
    conv_bg_clip = ChrisFuncs.SigmaClip(conv_bg_calc[2], median=False, sigma_thresh=3.0)
    conv_bg_avg = conv_bg_clip[1]
    conv_ap_sum = conv_ap_calc[0] - (conv_ap_calc[1] * conv_bg_avg)

    # Find difference between flux measued on convoled and unconvoled sersic maps
    ap_correction = np.nanmax([ 1.0, (sersic_ap_sum/conv_ap_sum) ])

    # Return aperture correction
    return ap_correction

# -----------------------------------------------------------------

def Sersic_LMfit(params, cutout, psf, mask, use_fft, lmfit=True):

    """
    Define LMfit convolved-sersic function
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

    # Generate sersic model given current paramters
    sersic_x, sersic_y = np.meshgrid( np.arange(cutout.shape[1]), np.arange(cutout.shape[0]) )
    sersic_model = astropy.modeling.models.Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
    sersic_map = sersic_model( sersic_x, sersic_y )

    # Convolve sersic model with PSF
    if use_fft==True:
        conv_map = astropy.convolution.convolve_fft(sersic_map, psf, normalize_kernel=True)
    elif use_fft==False:
        conv_map = astropy.convolution.convolve(sersic_map, psf, normalize_kernel=True)

    # Calculate residuals, filtered by mask
    residuals = cutout - conv_map
    mask[ np.where(mask==0.0) ] = np.nan
    residuals *= mask
    residuals.flatten()
    residuals = residuals[ np.where( np.isnan(residuals)==False ) ]

    # Return residuals
    if lmfit==True:
        return residuals**2.0
    elif lmfit==False:
        return residuals, cutout-conv_map

# -----------------------------------------------------------------

# EXAMPLE USE:
"""
input_dict = {}
input_dict['psf_path'] = '/home/herdata/spx7cjc/Beams/SPIRE_250.fits'
input_dict['cutout'] = astropy.io.fits.getdata('/home/saruman/spx7cjc/DustPedia/SPIRE/Cutouts/DustPedia/NGC4030_SPIRE_250.fits')
input_dict['pix_arcsec'] = 6.0
input_dict['semimaj_pix'] = 41.0
input_dict['axial_ratio'] = 1.1795263352195566
input_dict['angle'] = 115.16660752050387
input_dict['centre_i'] = 300.0
input_dict['centre_j'] = 300.0
input_dict['annulus_inner'] = 1.25
input_dict['annulus_outer'] = 1.601

ap_corr = StandaloneApCorrect(input_dict)
"""

# -----------------------------------------------------------------