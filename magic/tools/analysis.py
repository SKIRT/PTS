#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import copy

# Import Astromagic modules
from . import fitting
from . import plotting
from ..core import regions
from . import statistics
from . import coordinates
from . import cropping
from . import interpolation
from ..core.box import Box
from ..core.source import Source
from ..core.vector import Position, Extent

# Import astronomical modules
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.table import Table
from photutils import find_peaks
from find_galaxy import find_galaxy
from astropy import log
from photutils import daofind
from astropy.stats import sigma_clipped_stats

# *****************************************************************

def find_source_daofind(frame, center, radius, angle, config, track_record):

    """
    This function ...
    :param data:
    :return:
    """

    # TODO: FIX THIS FUNCTION

    threshold_sigmas = 5.0

    # Calculate the sigma-clipped statistics of the data
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    result_table = daofind(data - median, fwhm=3.0, threshold=threshold_sigmas*std)

    result_table.rename_column('xcentroid', 'x_peak')
    result_table.rename_column('ycentroid', 'y_peak')

    # If requested, make a plot with the source(s) indicated
    if plot: plotting.plot_peaks(data, result_table['x_peak'], result_table['y_peak'], radius=4.0)

    # Return the list of source positions
    #return result_table, median

    source = []
    return source

# *****************************************************************

def find_source_iraf(frame, center, radius, angle, config, track_record):

    """
    This function ...
    :param data:
    :return:
    """

    # TODO: FIX THIS FUNCTION

# *****************************************************************

def fit_model_to_source(source, config, track_record=None, level=0):

    """
    This function searches for sources ...
    :param box:
    :param x_peak:
    :param y_peak:
    :param x_min:
    :param y_min:
    :param plot:
    :return:
    """

    # Find source
    if config.use_center_or_peak == "center": position = source.center
    elif config.use_center_or_peak == "peak": position = source.peak
    else: raise ValueError("Invalid option (should be 'center' or 'peak')")

    # If the box is too small, don't bother looking for stars (anymore)
    if source.cutout.xsize < config.minimum_pixels or source.cutout.ysize < config.minimum_pixels: return None, None

    # Estimate and subtract the background of the source
    if not source.is_subtracted:

        # median, mean, polynomial or interpolation
        source.estimate_background(config.background_est_method, config.sigma_clip_background)
        source.subtract_background()

    # Get the model name
    model_name = config.model_names[level]

    # Fit the model to the background-subtracted box
    model = source.subtracted.fit_model(position, model_name, amplitude=source.cutout.value(source.peak))

    # If the amplitude is negative, the model is invalid
    if model.amplitude < 0:

        log.warning("Source profile has negative amplitude")
        return None, None

    # Calculate the difference between the mean position of the model and the position of the center / peak
    difference = fitting.center(model) - position

    # If ...
    if difference.norm > config.max_model_offset:

        #log.warning("Center of source profile lies outside of the region where it was expected (distance = "+str(difference.norm))

        # Show a plot for debugging
        if config.debug.model_offset:
            rel_peak = source.cutout.rel_position(source.peak)
            plotting.plot_peak_model(source.cutout, rel_peak.x, rel_peak.y, model)

        # Create a new zoomed-in source
        source = source.zoom(config.zoom_factor)

        # Estimate and subtract the background
        source.estimate_background(config.background_est_method, config.sigma_clip_background)
        source.subtract_background()

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Try again (iterative procedure of zooming in, stops if the size of the cutout becomes too small)
        return fit_model_to_source(source, config, track_record, level)

    # The fit succeeded
    else:

        # Show a plot for debugging
        if config.debug.success:
            rel_peak = source.cutout.rel_position(source.peak)
            plotting.plot_peak_model(source.cutout, rel_peak.x, rel_peak.y, model)

        # Return the model
        return source, model

# *****************************************************************

def estimate_background(data, mask, interpolate=True, sigma_clip=True):

    """
    This function ...
    :param data:
    :param mask:
    :return:
    """

    # Sigma clipping
    if sigma_clip: mask = statistics.sigma_clip_mask(data, sigma_level=3.0, mask=mask)

    # Decide whether to interpolate the background or to calculate a single median background value
    if interpolate: background = interpolation.in_paint(data, mask)

    else:

        # Calculate the median
        median = np.ma.median(np.ma.masked_array(data, mask=mask))

        # Create the background array
        background = np.full(data.shape, median)

    # Return the background
    return background, mask

# *****************************************************************

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
    x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

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

# *****************************************************************

def find_galaxy_orientation(data, region, plot=False):

    """
    This function ...
    :param data:
    :param region:
    :param plot:
    :return:
    """

    # TODO: improve this function, check the documentation of the find_galaxy class to improve the fit, instead of using cropping to 'solve' fitting problems

    # Verify that the region contains only one shape
    assert len(region) == 1, "The region can only contain one shape"
    shape = region[0]

    x_position = shape.coord_list[0]
    y_position = shape.coord_list[1]

    # Look for the galaxy orientation
    orientation = find_galaxy(data[::-1,:], quiet=True, plot=plot)
    if plot: plt.show()

    # The length of the major axis of the ellipse
    major = orientation.majoraxis

    # The width and heigth of the ellips
    width = major
    height = major * (1 - orientation.eps)

    if not np.isclose(x_position, orientation.ypeak, rtol=0.02) or not np.isclose(y_position, orientation.xpeak, rtol=0.02):

        x_size = data.shape[1]
        y_size = data.shape[0]
        size = max(x_size, y_size)

        smaller_data, x_min, x_max, y_min, y_max = cropping.crop(data, x_position, y_position, size/4.0, size/4.0)

        # Again look for the galaxy orientation
        orientation = find_galaxy(smaller_data[::-1,:], quiet=True, plot=plot)
        if plot: plt.show()

        # The length of the major axis of the ellipse
        major = orientation.majoraxis

        # The width and heigth of the ellips
        width = major
        height = major * (1 - orientation.eps)

        # Correct the center coordinate for the cropping
        orientation.ypeak += x_min
        orientation.xpeak += y_min

        if not np.isclose(x_position, orientation.ypeak, rtol=0.02) or not np.isclose(y_position, orientation.xpeak, rtol=0.02):
            log.warning("Could not find a galaxy at the specified position")

    # Return the parameters of the galaxy
    return (orientation.ypeak, orientation.xpeak, width, height, orientation.theta)

# *****************************************************************

def find_source(frame, center, radius, angle, config, track_record=None, special=False):

    """
    This function ...
    :return:
    """

    # Segmentation method
    if config.detection_method == "segmentation": return find_source_segmentation(frame, center, radius, angle, config, track_record, special=special)

    # Peaks method
    elif config.detection_method == "peaks": return find_source_peaks(frame, center, radius, angle, config, track_record, special=special)

    # DAOFIND source detection
    elif config.detection_method == "daofind": return find_source_daofind(frame, center, radius, angle, config, track_record)

    # IRAF's starfind algorithm
    elif config.detection_method == "iraf": return find_source_iraf(frame, center, radius, angle, config, track_record)

    # Unknown detection method
    else: raise ValueError("Unknown source detection method")

# *****************************************************************

def find_source_segmentation(frame, center, radius, angle, config, track_record=None, expansion_level=1, special=False):

    """
    This function ...
    :return:
    """

    # Create a source object
    source = Source(frame, center, radius, angle, config.background_outer_factor)

    # If always subtract background is enabled
    if config.always_subtract_background:

        # Subtract the background from the source
        source.estimate_background(config.background_est_method, sigma_clip=config.sigma_clip_background)
        source.subtract_background()

    # Create a kernel
    sigma = config.kernel.fwhm * gaussian_fwhm_to_sigma
    kernel_size = int(round(4.0 * config.kernel.cutoff_level))
    kernel = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size)

    # Create a mask for the center segment found for the source
    source.find_center_segment(config.threshold_sigmas, kernel=kernel, min_pixels=config.min_pixels)

    # If no center segment was found, subtract the background first
    if not np.any(source.mask) and not config.always_subtract_background:

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.no_segment_before: source.plot(title="No segment found, gradient background will be removed")

        # Subtract the background from the source
        source.estimate_background(config.background_est_method, sigma_clip=config.sigma_clip_background)
        source.subtract_background()

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.no_segment_after: source.plot(title="After removing gradient background")

    # If still no center segment was found, return without source
    if not np.any(source.mask):

        # Show a plot for debugging
        if config.debug.no_segment: source.plot(title="No center segment was found")

        return None

    # If the mask extents to the boundary of the cutout box and if enabled, expand the ellipse and repeat the procedure
    if source.mask.hits_boundary() and config.expand:

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.expand: plotting.plot_box(np.ma.masked_array(source.cutout, mask=source.mask), title="Masked segment hits boundary")

        # If the maximum expansion level has been reached, no source could be found
        if expansion_level >= config.max_expansion_level: return None

        else:

            # Calculate the expanded parameters
            radius *= config.expansion_factor
            expansion_level += 1

            # Repeat the procedure for the expanded ellipse
            return find_source_segmentation(frame, center, radius, angle, config, track_record=track_record, expansion_level=expansion_level)

    else:

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.success: plotting.plot_box(np.ma.masked_array(source.cutout, mask=source.mask), title="Masked segment doesn't hit boundary")

        # Dilate the mask if requested
        if config.dilate: source.mask = source.mask.dilated(config.connectivity, config.iterations)

        # Show a plot for debugging
        if config.debug.dilated: plotting.plot_box(np.ma.masked_array(source.cutout, mask=source.mask), title="With dilated mask")

        # Return the source
        return source

# *****************************************************************

def find_source_peaks(frame, center, radius, angle, config, track_record=None, level=0, special=False):

    """
    This function ...
    :param frame:
    :param center:
    :param radius:
    :param angle:
    :param config:
    :param level:
    :return:
    """

    # If the maximum or minimum level is reached, return without source
    if level < config.min_level or level > config.max_level: return None

    # Create a source object
    source = Source(frame, center, radius, angle, config.background_outer_factor)

    # If the frame is zero in this box, continue to the next object
    if not np.any(source.cutout): return None

    # If the box is too small, skip this object
    if source.cutout.xsize < config.minimum_pixels or source.cutout.ysize < config.minimum_pixels: return None

    # If always subtract background is enabled
    if config.always_subtract_background:

        source.estimate_background(config.background_est_method, config.sigma_clip_background)
        source.subtract_background()

    # Find the location of peaks in the box (do not remove gradient yet for performance reasons)
    peaks = source.locate_peaks(config.threshold_sigmas)

    # If no peaks could be detected, remove a potential background gradient from the box before detection
    if len(peaks) == 0 and not config.always_subtract_background:

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.zero_peaks_before or special: source.plot(title="0 peaks, gradient background will be removed")

        # Estimate and subtract the background (remove the background gradient)
        source.estimate_background(config.background_est_method, config.sigma_clip_background)
        source.subtract_background()

        # Find the location of peaks in the box
        peaks = source.locate_peaks(config.threshold_sigmas)

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.zero_peaks_after or special: source.plot(title=str(len(peaks)) + " peak(s) found after removing gradient background", peaks=peaks)

    # If no sources were detected
    if len(peaks) == 0:

        # Add a snapshot of the source to the track record
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.zero_peaks or special: source.plot(title="0 peaks")

        # If the level was negative, no source can be found
        if level < 0: return None

        # Scale the ellipse in which to look for a source
        radius *= config.scale_factor

        # Find a source in the zoomed-out region
        return find_source_peaks(frame, center, radius, angle, config, track_record=track_record, level=level+1)

    # If more than one source was detected
    elif len(peaks) > 1:

        # Add a snapshot of the source to the track record
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.more_peaks or special: source.plot(title="More peaks", peaks=peaks)

        # If the level was positive, no source can be found
        if level > 0: return None

        # Scale the ellipse in which to look for a source
        radius /= config.scale_factor

        # Find a source in the zoomed-in region
        return find_source_peaks(frame, center, radius, angle, config, track_record=track_record, level=level-1)

    # If one source was detected
    elif len(peaks) == 1:

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.one_peak or special: source.plot(title="1 peak")

        # Get the x and y coordinate of the peak
        x_peak = peaks[0].x
        y_peak = peaks[0].y

        # Check whether peak position corresponds to the center of the cutout
        if not (np.isclose(x_peak, center.x, atol=config.peak_offset_tolerance) and np.isclose(y_peak, center.y, atol=config.peak_offset_tolerance)):

            # Show a plot for debugging
            if config.debug.off_center or special: source.plot(title="Peak and center position do not match")

            # No source was found
            return None

        # Else, return the source
        else: return source

# *****************************************************************
