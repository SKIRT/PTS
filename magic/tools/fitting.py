#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import warnings
import numpy as np
from scipy import ndimage

# Import astronomical modules
from astropy.modeling import models, fitting

# Import Astromagic modules
from . import general
from . import statistics
from ..core.vector import Position, Extent

# *****************************************************************

def fit_2D_model(data, mask, background, model='Gaussian', x_center=None, y_center=None, radius=None, x_shift=0.0, y_shift=0.0, pixel_deviation=0.5, upsample_factor=1.0):

    """
    This function fits a 2D model to ...
    :param data:
    :param mask:
    :param background:
    :param model:
    :param x_center:
    :param y_center:
    :param radius:
    :param x_shift:
    :param y_shift:
    :param pixel_deviation:
    :param upsample_factor:
    :return:
    """

    # If the box around the star has to be upsampled for fitting
    if upsample_factor > 1.0:

        data = ndimage.interpolation.zoom(data, zoom=upsample_factor)
        mask = ndimage.interpolation.zoom(mask, zoom=upsample_factor)
        background = ndimage.interpolation.zoom(background, zoom=upsample_factor)

        x_center *= upsample_factor
        y_center *= upsample_factor
        radius *= upsample_factor
        pixel_deviation *= upsample_factor

    # Subtract the background from the box
    data_without_background = data - background

    # Don't keep the center fixed; this did not work
    fixed_center = False

    # Define the fitting functions
    fitting_functions = {'Gaussian': fit_2D_Gaussian, 'Airy': fit_2D_Airy, 'Moffat': fit_2D_Moffat, 'MexicanHat': fit_2D_MexicanHat}

    # Do the fitting, obtain a model function
    model_function = fitting_functions[model](data_without_background, center=(x_center, y_center),
                                              fixed_center=fixed_center, deviation_center=pixel_deviation,
                                              radius=radius, x_shift=x_shift, y_shift=y_shift,
                                              zoom_factor=upsample_factor, mask=mask)

    # Return the model
    return model_function

# *****************************************************************

def fit_polynomial(data, degree, mask=None, sigma_clip_background=False):

    """
    This function ...
    :param box:
    :param degree:
    :param x_shift:
    :param y_shift:
    :param mask:
    :return:
    """

    if sigma_clip_background: mask = statistics.sigma_clip_mask(data, sigma_level=3.0, mask=mask)

    # Fit the data using astropy.modeling
    poly_init = models.Polynomial2D(degree=degree)
    fit_model = fitting.LevMarLSQFitter()

    # Split x, y and z values that are not masked
    x_values, y_values, z_values = general.split_xyz(data, mask=mask, arrays=True)

    # Ignore model linearity warning from the fitter
    with warnings.catch_warnings():

        warnings.simplefilter('ignore')
        poly = fit_model(poly_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Return the polynomial model and the new mask
    if sigma_clip_background: return poly, mask

    # Return the polynomial model
    else: return poly

# *****************************************************************

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

# *****************************************************************

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

# *****************************************************************

def fit_2D_Gaussian(box, center=None, fixed_center=False, max_center_offset=None, sigma=None, zoom_factor=1.0, mask=None, amplitude=None):

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

# *****************************************************************

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

# *****************************************************************

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

# *****************************************************************

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

# *****************************************************************

def center(model):

    """
    This function ...
    :param model:
    :return:
    """

    if isinstance(model, models.Gaussian2D): return Position(x=model.x_mean.value, y=model.y_mean.value)
    elif isinstance(model, models.AiryDisk2D): return Position(x=model.x_0.value, y=model.y_0.value)
    else: raise ValueError("Unsupported model type")

# *****************************************************************

def sigma(model):

    """
    This function ...
    :param model:
    :return:
    """

    if isinstance(model, models.Gaussian2D): return Extent(x=model.x_stddev.value, y=model.y_stddev.value).norm
    elif isinstance(model, models.AiryDisk2D): return 0.42 * model.radius * 0.81989397882
    else: raise ValueError("Unsupported model type")

# *****************************************************************

def fwhm(model):

    """
    This function ...
    :param model:
    :return:
    """

    return 2.355 * sigma(model)

# *****************************************************************
