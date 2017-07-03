#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.analysis.stars Contains functions for fitting stars.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import ndimage

# Import the relevant PTS classes and modules
from ..tools import fitting, plotting, coordinates, cropping

# -----------------------------------------------------------------

# THIS FUNCTION IS NOT USED ANYMORE, BUT CONTAINS UPSAMPLING CODE POTENTIALLY USEFUL LATER
def make_star_model(shape, data, annuli_mask, fit_mask, background_outer_sigmas, fit_sigmas,
                    model_name, upsample_factor=1.0, interpolate_background=True, sigma_clip_background=True, plot=False):

    """
    This function ...
    :param shape:
    :param data:
    :param annuli_mask:
    :param fit_mask:
    :param background_inner_sigmas:
    :param background_outer_sigmas:
    :param fit_sigmas:
    :param upsample_factor:
    :param interpolate_background:
    :param sigma_clip_background:
    :param plot:
    :return:
    """

    # Get the shape's parameters
    x_center, y_center, x_radius, y_radius, _ = regions.ellipse_parameters(shape)

    # Set the radii for cutting out the background box
    radius = 0.5*(x_radius + y_radius)
    x_radius_outer = background_outer_sigmas*x_radius
    y_radius_outer = background_outer_sigmas*y_radius

    # Cut out the background
    background, x_min_back, x_max_back, y_min_back, y_max_back = cropping.crop(data, x_center, y_center, x_radius_outer, y_radius_outer)

    # Cut out the mask for the background
    background_mask = cropping.crop_check(annuli_mask, x_min_back, x_max_back, y_min_back, y_max_back)

    # Set the radii for cutting out the box for fitting
    x_radius_fitting = fit_sigmas*x_radius
    y_radius_fitting = fit_sigmas*y_radius

    # Cut out a box of selected frame around the star
    star, x_min, x_max, y_min, y_max = cropping.crop(data, x_center, y_center, x_radius_fitting, y_radius_fitting)

    # If the cropped region contains only one pixel row or column, a star model cannot be made
    if star.shape[0] == 1 or star.shape[1] == 1: return False, shape, None, None

    # Cut out the mask for fitting
    star_mask = fit_mask[y_min:y_max, x_min:x_max]

    # Estimate the background
    background_mask_beforeclipping = np.copy(background_mask)
    est_background, background_mask = estimate_background(background, background_mask, interpolate=interpolate_background, sigma_clip=sigma_clip_background)

    # Crop the interpolated background to the frame of the box
    star_background = cropping.crop_check(est_background, x_min-x_min_back, x_max-x_min_back, y_min-y_min_back, y_max-y_min_back)

    # Calculate the relative coordinates of the center
    x_center_rel, y_center_rel = coordinates.relative_coordinate(x_center, y_center, x_min, y_min)

    # Fit the star
    model_function = fitting.fit_2D_model(star, star_mask, star_background, model=model_name, x_center=x_center_rel,
                                          y_center=y_center_rel, radius=radius, x_shift=x_min, y_shift=y_min,
                                          upsample_factor=upsample_factor, pixel_deviation=0.5)

    # Evaluate the model
    evaluated_model = fitting.evaluate_model(model_function, x_min, x_max, y_min, y_max, x_delta=1.0/upsample_factor, y_delta=1.0/upsample_factor)

    # Check for succesful fit
    success = (np.isclose(model_function.x_stddev.value, x_radius, rtol=0.2) and np.isclose(model_function.y_stddev.value, y_radius, rtol=0.2))

    if success:

        if upsample_factor > 1.0: evaluated_model = ndimage.interpolation.zoom(evaluated_model, zoom=1.0/upsample_factor)

        # Plot
        if plot: plotting.plot_star_model(background=np.ma.masked_array(background,mask=background_mask_beforeclipping),
                                          background_clipped=np.ma.masked_array(background,mask=background_mask),
                                          est_background=est_background,
                                          star=np.ma.masked_array(star,mask=star_mask),
                                          est_background_star= star_background,
                                          fitted_star=evaluated_model)

        # Adjust the parameters of the shape to the model of this star
        shape.coord_list[0] = model_function.x_mean.value
        shape.coord_list[1] = model_function.y_mean.value
        shape.coord_list[2] = model_function.x_stddev.value
        shape.coord_list[3] = model_function.y_stddev.value

    # Return ...
    return success, shape, evaluated_model, (x_min, x_max, y_min, y_max)

# -----------------------------------------------------------------
