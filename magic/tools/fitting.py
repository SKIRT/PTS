#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.fitting Contains functions used for fitting models to two-dimensional data.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import warnings
import numpy as np
from scipy import stats

# Import astronomical modules
from astropy.modeling import models, fitting

# Import the relevant PTS classes and modules
from . import general, statistics
from ..basics.vector import Position, Extent
from ..basics.coordinate import PixelCoordinate

# -----------------------------------------------------------------

def linear_regression(x, y):

    """
    This function ...
    :param x:
    :param y:
    :return:
    """

    # Fit
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    # Return
    return slope, intercept

# -----------------------------------------------------------------

def get_linear_fit_parameters(x, y, xlog=False, ylog=False):

    """
    This function ...
    :param x:
    :param y:
    :param xlog:
    :param ylog:
    :return:
    """

    # Convert?
    if xlog: x = np.log10(x)
    if ylog: y = np.log10(y)

    # Keep only finite values
    valid = np.isfinite(x) * np.isfinite(y)
    x = x[valid]
    y = y[valid]

    # Perform fit
    slope, intercept = linear_regression(x, y)

    # Return
    return slope, intercept

# -----------------------------------------------------------------

def get_linear_values(x, slope, intercept, xlog=False, ylog=False):

    """
    Thins function ...
    :param x:
    :param slope:
    :param intercept:
    :param xlog:
    :param ylog:
    :return:
    """

    # Calculate new values
    if xlog: x = np.log10(x)
    y = slope * x + intercept

    # Convert?
    if ylog: y = 10 ** y

    # Return
    return y

# -----------------------------------------------------------------

def get_linear_fitted_values(x, y, new_x, xlog=False, ylog=False, return_parameters=False):

    """
    This function ...
    :param x:
    :param y:
    :param new_x:
    :param xlog:
    :param ylog:
    :param return_parameters:
    :return:
    """

    # Get the parameters
    slope, intercept = get_linear_fit_parameters(x, y, xlog=xlog, ylog=ylog)

    # Get the values
    new = get_linear_values(new_x, slope, intercept, xlog=xlog, ylog=ylog)

    # Return
    if return_parameters: return new, slope, intercept
    else: return new

# -----------------------------------------------------------------

def fit_polynomial(data, degree, mask=None, sigma_clip_background=False, show_warnings=False,
                   fitter="levenberg-marquardt", zero_order=None):

    """
    This function ...
    :param data:
    :param degree:
    :param mask:
    :param sigma_clip_background:
    :param show_warnings:
    :param fitter:
    :return:
    """

    if sigma_clip_background: mask = statistics.sigma_clip_mask(data, sigma_level=3.0, mask=mask)

    # Fit the data using astropy.modeling
    poly_init = models.Polynomial2D(degree=degree)

    # Set initial values
    if zero_order is not None: poly_init.c0_0 = zero_order
    #poly_init.c1_0 = 1.0
    #poly_init.c2_0 = 2.0
    #poly_init.c3_0 = 3.0
    #poly_init.c0_1 = 4.0
    #poly_init.c0_2 = 5.0
    #poly_init.c0_3 = 6.0
    #poly_init.c1_1 = 7.0
    #poly_init.c1_2 = 8.0
    #poly_init.c2_1 = 9.0

    # Get the fitter
    if fitter == "levenberg-marquardt": fit_model = fitting.LevMarLSQFitter()
    elif fitter == "linear":
        if degree > 1: raise ValueError("Cannot use linear fitter for polynomials with degree > 1")
        fit_model = fitting.LinearLSQFitter()
    elif fitter == "simplex": fit_model = fitting.SimplexLSQFitter()
    elif fitter == "sequential_least_squares": fit_model = fitting.SLSQPLSQFitter()
    else: raise ValueError("Invalid vlaue for 'fitter'")

    # Split x, y and z values that are not masked
    x_values, y_values, z_values = general.split_xyz(data, mask=mask, arrays=True)

    # Ignore model linearity warning from the fitter
    if show_warnings: poly = fit_model(poly_init, x_values, y_values, z_values)
    else:
        with warnings.catch_warnings():

            warnings.simplefilter('ignore')
            poly = fit_model(poly_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Return the polynomial model and the new mask
    if sigma_clip_background: return poly, mask

    # Return the polynomial model
    else: return poly

# -----------------------------------------------------------------

def all_zero_parameters(polynomial):

    """
    This function ...
    :param polynomial:
    :return:
    """

    from ...core.tools import sequences
    return sequences.all_equal_to(polynomial.parameters, 0.0)

# -----------------------------------------------------------------

def evaluate_model(model, x_min, x_max, y_min, y_max, x_delta=1, y_delta=1):

    """
    This function ...
    :param model:
    :param x_min:
    :param x_max:
    :param y_min:
    :param y_max:
    :param x_delta:
    :param y_delta:
    :return:
    """

    # Create x and y meshgrid for evaluating the model
    y_plotvalues, x_plotvalues = np.mgrid[y_min:y_max:y_delta, x_min:x_max:x_delta]

    # Evaluate the model
    evaluated_model = model(x_plotvalues, y_plotvalues)

    # Return the evaluated data
    return evaluated_model

# -----------------------------------------------------------------

# NOT USED CURRENTLY
def fit_two_2D_Gaussians(box, x_shift=0.0, y_shift=0.0, zoom_factor=1.0, mask=None):

    """
    This function ...
    :param box:
    :param x_shift:
    :param y_shift:
    :param zoom_factor:
    :param mask:
    :return:
    """

    # Get the dimensions of the box
    box_ysize = box.shape[0]
    box_xsize = box.shape[1]

    upperright_x = 0.75*box_xsize
    upperright_y = 0.75*box_ysize

    lowerleft_x = 0.25*box_xsize
    lowerleft_y = 0.25*box_ysize

    init_x_stddev = 0.2*box_xsize
    init_y_stddev = 0.2*box_ysize

    two_gaussians_init = models.Gaussian2D(amplitude=1., x_mean=upperright_x, y_mean=upperright_y, x_stddev=init_x_stddev,
                                           y_stddev=init_y_stddev) + \
                         models.Gaussian2D(amplitude=1., x_mean=lowerleft_x, y_mean=lowerleft_y, x_stddev=init_x_stddev,
                                           y_stddev=init_y_stddev)

    fit_model = fitting.LevMarLSQFitter()

    x_values = []
    y_values = []
    z_values = []

    for x in range(box_xsize):
        for y in range(box_ysize):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:

                x_values.append(x)
                y_values.append(y)
                z_values.append(box[y,x])

    # Ignore model linearity warning from the fitter
    with warnings.catch_warnings():

        warnings.simplefilter('ignore')
        two_gaussians = fit_model(two_gaussians_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Adjust the position of the model to a different coordinate frame
    if zoom_factor > 1.0:

        two_gaussians.x_mean_0.value = two_gaussians.x_mean_0.value/zoom_factor + x_shift
        two_gaussians.y_mean_0.value = two_gaussians.y_mean_0.value/zoom_factor + y_shift
        two_gaussians.x_stddev_0.value /= zoom_factor
        two_gaussians.y_stddev_0.value /= zoom_factor

        two_gaussians.x_mean_1.value = two_gaussians.x_mean_1.value/zoom_factor + x_shift
        two_gaussians.y_mean_1.value = two_gaussians.y_mean_1.value/zoom_factor + y_shift
        two_gaussians.x_stddev_1.value /= zoom_factor
        two_gaussians.y_stddev_1.value /= zoom_factor

    else:

        two_gaussians.x_mean_0.value += x_shift
        two_gaussians.y_mean_0.value += y_shift
        two_gaussians.x_mean_1.value += x_shift
        two_gaussians.y_mean_1.value += y_shift

    # Return the model
    return two_gaussians

# -----------------------------------------------------------------

def fit_2D_ShiftedGaussian(box, center=None, fixed_center=False, max_center_offset=None, sigma=None, zoom_factor=1.0, mask=None, amplitude=None):

    """
    This function ...
    :param box:
    :param center:
    :param fixed_center:
    :param max_center_offset:
    :param sigma:
    :param zoom_factor:
    :param mask:
    :param amplitude:
    :return:
    """

    # Compound model class that represent a Gaussian function that can be shifted up and down
    ShiftedGaussian = models.Gaussian2D + models.Const2D

    # Parameters are: amplitude_0, x_mean_0, y_mean_0, x_stddev_0, y_stddev_0, theta_0, amplitude_1

# -----------------------------------------------------------------

def fit_2D_Gaussian(box, center=None, fixed_center=False, max_center_offset=None, sigma=None, zoom_factor=1.0, mask=None, amplitude=None, max_sigma_offset=None):

    """
    This function ...
    :param box:
    :param center:
    :param fixed_center:
    :param deviation_center:
    :param radius:
    :param x_shift:
    :param y_shift:
    :param zoom_factor:
    :param mask:
    :return:
    """

    # Get the dimensions of the box
    box_xsize = box.xsize
    box_ysize = box.ysize

    # Set the initial guess for the center of the model (the one that is specified, otherwise the center of the box)
    init_xmean = center.x if center is not None else 0.5*(box.xsize-1)
    init_ymean = center.y if center is not None else 0.5*(box.ysize-1)

    # Set the initial guess for the width of the model (the one that is specified, otherwise one tenth of the size of the box)
    init_x_stddev = sigma if sigma is not None else 0.1*box.xsize
    init_y_stddev = sigma if sigma is not None else 0.1*box.ysize

    # Initialize an empty dictionary to specify fixed parameters
    fixed_parameters = {'theta': 0.0}

    if fixed_center:

        fixed_parameters['x_mean'] = True
        fixed_parameters['y_mean'] = True

    # Initialize an empty dictionary to specify bounds
    bounds = {}

    if max_center_offset is not None:

        bounds['x_mean'] = [init_xmean-max_center_offset, init_xmean+max_center_offset]
        bounds['y_mean'] = [init_ymean-max_center_offset, init_ymean+max_center_offset]

    if max_sigma_offset is not None:

        bounds["x_stddev"] = [init_x_stddev - max_sigma_offset, init_x_stddev + max_sigma_offset]
        bounds["y_stddev"] = [init_y_stddev - max_sigma_offset, init_y_stddev + max_sigma_offset]

    # Define the 'tied' dictionary to specify that the y_stddev should vary along with x_stddev
    tied = {'y_stddev': (lambda model: model.x_stddev)}

    # Fit the data using astropy.modeling
    gaussian_init = models.Gaussian2D(amplitude=1., x_mean=init_xmean, y_mean=init_ymean, x_stddev=init_x_stddev,
                                      y_stddev=init_y_stddev, fixed=fixed_parameters, bounds=bounds, tied=tied)
    fit_model = fitting.LevMarLSQFitter()

    x_values = []
    y_values = []
    z_values = []

    for x in range(box_xsize):
        for y in range(box_ysize):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:

                x_values.append(x)
                y_values.append(y)
                z_values.append(box[y,x])

    # Ignore model linearity warning from the fitter
    with warnings.catch_warnings():

        warnings.simplefilter('ignore')
        gaussian = fit_model(gaussian_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Fix negative sigmas
    if gaussian.x_stddev.value < 0: gaussian.x_stddev.value = -gaussian.x_stddev.value
    if gaussian.y_stddev.value < 0: gaussian.y_stddev.value = -gaussian.y_stddev.value

    # Adjust the position of the model to a different coordinate frame
    if zoom_factor > 1.0:

        gaussian.x_mean.value = gaussian.x_mean.value/zoom_factor
        gaussian.y_mean.value = gaussian.y_mean.value/zoom_factor
        gaussian.x_stddev.value /= zoom_factor
        gaussian.y_stddev.value /= zoom_factor

    # Return the Gaussian model
    return gaussian

# -----------------------------------------------------------------

def fit_2D_Airy(box, center=None, fixed_center=False, max_center_offset=None, radius=None, zoom_factor=1.0, mask=None, amplitude=None):

    """
    This function ...
    :param box:
    :param center:
    :param fixed_center:
    :param deviation_center:
    :param radius:
    :param x_shift:
    :param y_shift:
    :param zoom_factor:
    :param mask:
    :return:
    """

    # Get the dimensions of the box
    box_ysize = box.shape[0]
    box_xsize = box.shape[1]

    # Set the initial guess for the center of the model (the one that is specified, otherwise the center of the box)
    init_x0 = center.x if center is not None else 0.5*(box_xsize-1)
    init_y0 = center.y if center is not None else 0.5*(box_ysize-1)

    # Set the initial radius for the model (the one that is specified, otherwise one tenth of the width of the box)
    init_radius = radius if radius is not None else 0.1*box_xsize

    # Initialize an empty dictionary to specify fixed parameters
    fixed_parameters = {}

    if fixed_center:

        fixed_parameters['x_0'] = True
        fixed_parameters['y_0'] = True

    # Initialize an empty dictionary to specify bounds
    bounds = {}

    if max_center_offset is not None:

        bounds['x_mean'] = [init_x0-max_center_offset, init_x0+max_center_offset]
        bounds['y_mean'] = [init_y0-max_center_offset, init_y0+max_center_offset]

    # Fit the data using astropy.modeling
    airy_init = models.AiryDisk2D(amplitude=1., x_0=init_x0, y_0=init_y0, radius=init_radius, fixed=fixed_parameters, bounds=bounds)
    fit_model = fitting.LevMarLSQFitter()

    x_values = []
    y_values = []
    z_values = []

    for x in range(box_xsize):
        for y in range(box_ysize):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:

                x_values.append(x)
                y_values.append(y)
                z_values.append(box[y,x])

    # Ignore model linearity warning from the fitter
    with warnings.catch_warnings():

        warnings.simplefilter('ignore')
        airy = fit_model(airy_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Adjust the position of the model to a different coordinate frame
    if zoom_factor > 1.0:

        airy.x_0.value = airy.x_0.value/zoom_factor
        airy.y_0.value = airy.y_0.value/zoom_factor
        airy.radius /= zoom_factor

    # Return the fitted two-dimensional Airy Disk model
    return airy

# -----------------------------------------------------------------

def fit_2D_Moffat(box, center=None, fixed_center=False, deviation_center=None, x_shift=0.0, y_shift=0.0,
                  zoom_factor=1.0, mask=None):

    """
    This function ...
    :param box:
    :param center:
    :param fixed_center:
    :param deviation_center:
    :param x_shift:
    :param y_shift:
    :param zoom_factor:
    :param mask:
    :return:
    """

    # Get the dimensions of the box
    box_ysize = box.shape[0]
    box_xsize = box.shape[1]

    # Set the initial guess for the center of the model (the one that is specified, otherwise the center of the box)
    init_x0 = center[0] if center is not None else 0.5*(box_xsize-1)
    init_y0 = center[1] if center is not None else 0.5*(box_ysize-1)

    # Initialize an empty dictionary to specify fixed parameters
    fixed_parameters = {}

    if fixed_center:

        fixed_parameters['x_0'] = True
        fixed_parameters['y_0'] = True

    # Initialize an empty dictionary to specify bounds
    bounds = {}

    if deviation_center is not None:

        bounds['x_mean'] = [init_x0-deviation_center, init_x0+deviation_center]
        bounds['y_mean'] = [init_y0-deviation_center, init_y0+deviation_center]

    # Fit the data using astropy.modeling
    moffat_init = models.Moffat2D(amplitude=1., x_0=init_x0, y_0=init_y0, gamma=1.0, alpha=1.0, fixed=fixed_parameters, bounds=bounds)
    fit_model = fitting.LevMarLSQFitter()

    x_values = []
    y_values = []
    z_values = []

    for x in range(box_xsize):
        for y in range(box_ysize):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:

                x_values.append(x)
                y_values.append(y)
                z_values.append(box[y,x])

    # Ignore model linearity warning from the fitter
    with warnings.catch_warnings():

        warnings.simplefilter('ignore')
        moffat = fit_model(moffat_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Adjust the position of the model to a different coordinate frame
    if zoom_factor > 1.0:

        moffat.x_0.value = moffat.x_0.value/zoom_factor + x_shift
        moffat.y_0.value = moffat.y_0.value/zoom_factor + y_shift

    else:

        moffat.x_0.value += x_shift
        moffat.y_0.value += y_shift

    # Return the fitted two-dimensional Moffat model
    return moffat

# -----------------------------------------------------------------

def fit_2D_MexicanHat(box, center=None, fixed_center=False, deviation_center=None, radius=None, x_shift=0.0,
                      y_shift=0.0, zoom_factor=1.0, mask=None):

    """
    This function ...
    :param box:
    :param center:
    :param fixed_center:
    :param deviation_center:
    :param radius:
    :param x_shift:
    :param y_shift:
    :param zoom_factor:
    :param mask:
    :return:
    """

    # Get the dimensions of the box
    box_ysize = box.shape[0]
    box_xsize = box.shape[1]

    # Set the initial guess for the center of the model (the one that is specified, otherwise the center of the box)
    init_x0 = center[0] if center is not None else 0.5*(box_xsize-1)
    init_y0 = center[1] if center is not None else 0.5*(box_ysize-1)

    # Set the initial guess for the radius of the model (the one that is specified, otherwise one tenth of the width of the box)
    init_sigma = radius if radius is not None else 0.1*box_xsize

    # Initialize an empty dictionary to specify fixed parameters
    fixed_parameters = {}

    if fixed_center:

        fixed_parameters['x_0'] = True
        fixed_parameters['y_0'] = True

    # Initialize an empty dictionary to specify bounds
    bounds = {}

    if deviation_center is not None:

        bounds['x_mean'] = [init_x0-deviation_center, init_x0+deviation_center]
        bounds['y_mean'] = [init_y0-deviation_center, init_y0+deviation_center]

    # Fit the data using astropy.modeling
    mexicanhat_init = models.MexicanHat2D(amplitude=1., x_0=init_x0, y_0=init_y0, sigma=init_sigma, fixed=fixed_parameters, bounds=bounds)
    fit_model = fitting.LevMarLSQFitter()

    x_values = []
    y_values = []
    z_values = []

    for x in range(box_xsize):
        for y in range(box_ysize):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:

                x_values.append(x)
                y_values.append(y)
                z_values.append(box[y,x])

    # Ignore model linearity warning from the fitter
    with warnings.catch_warnings():

        warnings.simplefilter('ignore')
        mexicanhat = fit_model(mexicanhat_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Adjust the position of the model to a different coordinate frame
    if zoom_factor > 1.0:

        mexicanhat.x_0.value = mexicanhat.x_0.value/zoom_factor + x_shift
        mexicanhat.y_0.value = mexicanhat.y_0.value/zoom_factor + y_shift
        mexicanhat.sigma.value /= zoom_factor

    else:

        mexicanhat.x_0.value += x_shift
        mexicanhat.y_0.value += y_shift

    # Return the fitted two-dimensional Mexican Hat model
    return mexicanhat

# -----------------------------------------------------------------

def center(model):

    """
    This function ...
    :param model:
    :return:
    """

    if isinstance(model, models.Gaussian2D): return PixelCoordinate(x=model.x_mean.value, y=model.y_mean.value)
    elif isinstance(model, models.AiryDisk2D): return PixelCoordinate(x=model.x_0.value, y=model.y_0.value)
    else: raise ValueError("Unsupported model type: " + str(type(model)))

# -----------------------------------------------------------------

def sigma(model):

    """
    This function ...
    :param model:
    :return:
    """

    if isinstance(model, models.Gaussian2D): return Extent(x=model.x_stddev.value, y=model.y_stddev.value).norm
    elif isinstance(model, models.AiryDisk2D): return airy_radius_to_gaussian_sigma(model.radius)
    else: raise ValueError("Unsupported model type: " + str(type(model)))

# -----------------------------------------------------------------

def sigma_symmetric(model):

    """
    This function ...
    :param model: 
    :return: 
    """

    stddev = sigma(model)
    if stddev.x != stddev.y: raise ValueError("x and y stddev are not equal")
    return stddev.x

# -----------------------------------------------------------------

def airy_radius_to_gaussian_sigma(radius):

    """
    This function ...
    :param radius:
    :return:
    """

    return 0.42 * radius * 0.81989397882

# -----------------------------------------------------------------

def gaussian_sigma_to_airy_radius(sigma):

    """
    This function ...
    :param sigma:
    :return:
    """

    return sigma / (0.42 * 0.81989397882)

# -----------------------------------------------------------------

def fwhm_to_airy_radius(fwhm):

    """
    This function ...
    :param fwhm:
    :return:
    """

    sigma = statistics.fwhm_to_sigma * fwhm
    return gaussian_sigma_to_airy_radius(sigma)

# -----------------------------------------------------------------

def fwhm(model):

    """
    This function ...
    :param model:
    :return:
    """

    return statistics.sigma_to_fwhm * sigma(model)

# -----------------------------------------------------------------

def fwhm_symmetric(model):

    """
    This function ...
    :param model: 
    :return: 
    """

    return statistics.sigma_to_fwhm * sigma_symmetric(model)

# -----------------------------------------------------------------

def shift_model(model, x_shift, y_shift):

    """
    This function ...
    :param model:
    :param x_shift:
    :param y_shift:
    :return:
    """

    # If the model is a 2D Gaussian function
    if isinstance(model, models.Gaussian2D):

        model.x_mean += x_shift
        model.y_mean += y_shift

    # If the model is a 2D Airy Disk function
    elif isinstance(model, models.AiryDisk2D):

        model.x_0 += x_shift
        model.y_0 += y_shift

    # Unsupported models
    else: raise ValueError("Unsupported model (should be 'Gaussian2D' or 'AiryDisk2D'")

# -----------------------------------------------------------------

def shifted_model(model, x_shift, y_shift):

    """
    This function ...
    :param model:
    :return:
    """

    # Make a copy of the original model
    new_model = copy.deepcopy(model)

    # Shift the new model
    shift_model(new_model, x_shift, y_shift)

    # Return the new model
    return new_model

# -----------------------------------------------------------------
