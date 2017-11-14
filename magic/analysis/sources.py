#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.analysis.sources Contains functions for finding sources etc.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
from scipy import ndimage
import copy

# Import astronomical modules
from astropy.convolution import Gaussian2DKernel
from photutils import daofind
from astropy.stats import sigma_clipped_stats

# Import astronomical modules
import astropy.units as u
from astropy.coordinates import Angle
from photutils import source_properties, properties_table

# Import the relevant PTS classes and modules
from ..tools import fitting, plotting, statistics, coordinates, cropping, interpolation, masks
from ..core.detection import Detection
from ..region.ellipse import PixelEllipseRegion
from ...core.basics.log import log
from ..region import tools as regions
from ..basics.coordinate import PixelCoordinate
from ..basics.stretch import PixelStretch

# -----------------------------------------------------------------

def find_contours(data, segments, sigma_level):

    """
    This function ...
    :param data:
    :param segments:
    :param sigma_level:
    :return:
    """

    # Initialize a list for the contours
    contours = []

    # Get the segment properties
    # Since there is only one segment in the source.mask (the center segment), the props
    # list contains only one entry (one galaxy)
    properties_list = source_properties(data, segments)

    for properties in properties_list:

        # Obtain the position, orientation and extent
        position = PixelCoordinate(properties.xcentroid.value, properties.ycentroid.value)
        a = properties.semimajor_axis_sigma.value * sigma_level
        b = properties.semiminor_axis_sigma.value * sigma_level
        angle = properties.orientation.value # in radians
        angle = Angle(angle, u.rad)

        radius = PixelStretch(a, b)

        meta = {"text": str(properties.label)}

        # Create the contour
        contours.append(PixelEllipseRegion(position, radius, angle, meta=meta))

    # Return the contours
    return contours

# -----------------------------------------------------------------

def find_contour_mask(mask):

    """
    This function ...
    :param mask:
    :return:
    """

    x_min = 0
    x_max = mask.xsize - 1
    y_min = 0
    y_max = mask.ysize - 1

    from ..core.cutout import Cutout
    maskdata = mask.data
    fake_data = Cutout(maskdata.astype(float), x_min, x_max, y_min, y_max)
    contour = find_contour(fake_data, maskdata, sigma_level=3)
    return contour

# -----------------------------------------------------------------

def find_contour(box, mask, sigma_level):

    """
    This function ...
    :param box:
    :param mask:
    :param sigma_level:
    :return:
    """

    props = source_properties(box, mask)
    #tbl = properties_table(props)

    x_shift = box.x_min
    y_shift = box.y_min

    # Since there is only one segment in the self.source.mask (the center segment), the props
    # list contains only one entry (one galaxy)
    if len(props) == 0: return None
    properties = props[0]

    # Obtain the position, orientation and extent
    position = PixelCoordinate(properties.xcentroid.value + x_shift, properties.ycentroid.value + y_shift)
    a = properties.semimajor_axis_sigma.value * sigma_level
    b = properties.semiminor_axis_sigma.value * sigma_level
    angle = properties.orientation.value # in radians
    angle = Angle(angle, u.rad)

    #print("a", a)
    #print("b", b)

    # Get radius
    radius = PixelStretch(a, b)

    # Create and return the elliptical contour
    return PixelEllipseRegion(position, radius, angle)

# -----------------------------------------------------------------

def find_source_daofind(frame, ellipse, config, track_record, special=False):

    """
    This function ...
    :param data:
    :return:
    """

    # TODO: FIX THIS FUNCTION

    sigma_level = 5.0

    # Calculate the sigma-clipped statistics of the data
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    result_table = daofind(data - median, fwhm=3.0, threshold=sigma_level*std)

    result_table.rename_column('xcentroid', 'x_peak')
    result_table.rename_column('ycentroid', 'y_peak')

    # If requested, make a plot with the source(s) indicated
    if plot: plotting.plot_peaks(data, result_table['x_peak'], result_table['y_peak'], radius=4.0)

    # Return the list of source positions
    #return result_table, median

    source = []
    return source

# -----------------------------------------------------------------

def find_source_iraf(frame, ellipse, config, track_record, special=False):

    """
    This function ...
    :param data:
    :return:
    """

    # TODO: FIX THIS FUNCTION

# -----------------------------------------------------------------

def fit_model_to_source(source, config, track_record=None, level=0, special=False, stop_if_fail=False):

    """
    This function searches for sources ...
    :param source:
    :param config:
    :param track_record:
    :param level:
    :param special:
    :param stop_if_fail:
    :return:
    """

    # Find source
    if config.use_center_or_peak == "center": position = source.center
    elif config.use_center_or_peak == "peak": position = source.peak
    else: raise ValueError("Invalid option (should be 'center' or 'peak')")

    # If the box is too small, don't bother looking for stars (anymore)
    if source.cutout.xsize < config.minimum_pixels or source.cutout.ysize < config.minimum_pixels: return None, None

    # Estimate and subtract the background of the source
    if not source.has_background: source.estimate_background(config.background_est_method, config.sigma_clip_background)

    # Get the model name
    model_name = config.model_names[level]

    # Fit the model to the background-subtracted box
    try:
        model = source.subtracted.fit_model(position, model_name, amplitude=source.cutout.value(source.peak))
    except IndexError:
        log.debug("Index error occurred while fitting ...")
        log.debug("PEAK= (" + str(source.peak.x) + "," + str(source.peak.y) + ")")
        log.debug("source.cutout.x_min,y_min = " + str(source.cutout.x_min) + "," + str(source.cutout.y_min))
        log.debug("rel_peak = " + str(source.peak.x - source.cutout.x_min) + "," + str(source.peak.y - source.cutout.y_min))

        # TODO: NO SOLUTION YET AS TO WHY SOMETIMES THE PEAK POSITIONS ARE OUTSIDE OF THE SOURCE.CUTOUT
        return None, None

    # If the amplitude is negative, the model is invalid
    if model.amplitude < 0:

        log.warning("Model fitted to source has negative amplitude")
        return None, None

    # Calculate the difference between the mean position of the model and the position of the center / peak
    difference = fitting.center(model) - position

    # Failed fit
    if difference.norm > config.max_model_offset:

        # If stop if fail
        if stop_if_fail: return None, None

        # Show a plot for debugging
        if config.debug.model_offset or special:

            rel_peak = source.cutout.rel_position(source.peak)
            rel_model = fitting.shifted_model(model, -source.cutout.x_min, -source.cutout.y_min)
            plotting.plot_peak_model(source.cutout, rel_peak.x, rel_peak.y, rel_model, title="Center of source and peak do not match")

        # Set min pixels
        min_pixels = 4

        # Create a new zoomed-in source
        source = source.zoom(config.zoom_factor, min_xpixels=min_pixels, min_ypixels=min_pixels)

        # Check if we cannot zoom further
        if source.cutout.xsize <= min_pixels: stop_if_fail = True
        elif source.cutout.ysize <= min_pixels: stop_if_fail = True
        else: stop_if_fail = False

        # Estimate and subtract the background
        source.estimate_background(config.background_est_method, config.sigma_clip_background)

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Try again (iterative procedure of zooming in, stops if the size of the cutout becomes too small)
        return fit_model_to_source(source, config, track_record, level, special=special, stop_if_fail=stop_if_fail)

    # The fit succeeded
    else:

        # Show a plot for debugging
        if config.debug.success or special:

            rel_peak = source.cutout.rel_position(source.peak)
            rel_model = fitting.shifted_model(model, -source.cutout.x_min, -source.cutout.y_min)
            plotting.plot_peak_model(source.cutout, rel_peak.x, rel_peak.y, rel_model, title="Found a model that corresponds to the peak position")

        # Return the model
        return source, model

# -----------------------------------------------------------------

def estimate_background(data, mask, interpolate=True, sigma_clip=True):

    """
    This function ...
    :param data:
    :param mask:
    :param interpolate:
    :param sigma_clip:
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

# -----------------------------------------------------------------

def make_star_model(shape, data, annuli_mask, fit_mask, background_outer_sigmas, fit_sigmas,
                    model_name, upsample_factor=1.0, interpolate_background=True, sigma_clip_background=True, plot=False):

    """
    This function ...
    :param shape:
    :param data:
    :param annuli_mask:
    :param fit_mask:
    :param background_outer_sigmas:
    :param fit_sigmas:
    :param model_name:
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

def find_source(frame, ellipse, config, track_record=None, special=False):

    """
    This function ...
    :param frame:
    :param ellipse:
    :param track_record:
    :param special:
    :return:
    """

    # Segmentation method
    if config.detection_method == "segmentation": return find_source_segmentation(frame, ellipse, config, track_record, special=special)

    # Peaks method
    elif config.detection_method == "peaks": return find_source_peaks(frame, ellipse, config, track_record, special=special)

    # DAOFIND source detection
    elif config.detection_method == "daofind": return find_source_daofind(frame, ellipse, config, track_record, special=special)

    # IRAF's starfind algorithm
    elif config.detection_method == "iraf": return find_source_iraf(frame, ellipse, config, track_record, special=special)

    # Unknown detection method
    else: raise ValueError("Unknown source detection method")

# -----------------------------------------------------------------

def find_source_segmentation(frame, ellipse, config, track_record=None, expansion_level=1, special=False, sigma_level=None):

    """
    This function ...
    :param frame:
    :param ellipse:
    :param config:
    :param track_record:
    :param expansion_level:
    :param special:
    :param sigma_level:
    :return:
    """

    if special: log.debug("finding segmentation source, expansion level = " + str(expansion_level))

    # Allow for a custom sigma level
    sigma_level = config.sigma_level if sigma_level is None else sigma_level

    # Create a source object
    source = Detection.from_ellipse(frame, ellipse, config.background_outer_factor)

    # If the source cutout is zero or nan everywhere, return None (no source can be found here)
    if np.all(np.isnan(source.cutout)) or not np.any(source.cutout):
        if special: log.debug("no source can be found (cutout is zero or nan everywhere)")
        return None

    # If there are any nans, return None ??? yes, do we want this ? (temporary fix)
    if np.any(np.isnan(source.cutout)):

        #import os
        #from ..core import Frame
        #frame = Frame(source.cutout)
        #frame.saveto(os.path.join(os.getcwd(), "lalalalalal-nans.fits"))
        #plotting.plot_box(source.cutout)

        if special:
            log.debug("nans present in source cutout, setting corresponding pixels to zero")
            plotting.plot_box(source.cutout, title="cutout with nans")

        # Set nans zero
        source.cutout[np.isnan(source.cutout)] = 0.0

        if special: plotting.plot_box(source.cutout, title="nans replaced by 0.0")

    # If always subtract background is enabled
    if config.always_subtract_background:

        # Subtract the background from the source
        try: # weird error coming out for example with M81 GALEX FUV image (saturation detection)
            source.estimate_background(config.background_est_method, sigma_clip=config.sigma_clip_background)
        except:
            if special: log.debug("no source can be found (exception encountered while estimating background)")
            return None

    # Create a kernel
    sigma = config.kernel.fwhm * statistics.fwhm_to_sigma
    kernel_size = int(round(4.0 * config.kernel.cutoff_level))
    kernel = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size)
    kernel.normalize() # to suppress warning

    if special: log.debug("looking for center segment")

    # Create a mask for the center segment found for the source
    mask = source.find_center_segment(sigma_level, kernel=kernel, min_pixels=config.min_pixels)

    # If no center segment was found, subtract the background first
    if not np.any(mask) and not config.always_subtract_background:

        if special: log.debug("no center segment found")

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.no_segment_before or special: source.plot(title="No segment found, gradient background will be removed")

        # Subtract the background from the source
        try: # weird error coming out for example with M81 GALEX FUV image (saturation detection)
            source.estimate_background(config.background_est_method, sigma_clip=config.sigma_clip_background)
        except:
            if special: log.debug("no source can be found (exception encountered while estimating background)")
            return None

        if special: log.debug("looking for center segment again")

        # Search for a center segment again
        mask = source.find_center_segment(sigma_level, kernel=kernel, min_pixels=config.min_pixels)

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.no_segment_after or special: source.plot(title="After removing gradient background")

    # If still no center segment was found, return without source
    if not np.any(mask):

        if special: log.debug("still no center segment found")

        # Show a plot for debugging
        if config.debug.no_segment or special: source.plot(title="No center segment was found")

        # No source was found
        return None

    mask_without_appendages = mask.copy()

    # If overlapping is not allowed, see whether removing appendages helps by making it not overlap
    if not config.allow_overlap:

        # If the mask extents to the boundary of the cutout box en if enabled, apply binary opening to the mask to
        if masks.overlap(source.background_mask, mask) and config.remove_appendages:

            if special: log.debug("mask overlaps the background mask")

            # Show a plot for debugging
            if config.debug.overlap_before or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=mask), title="Overlapping mask before appendage removal")

            if special: log.debug("removing appendages")

            # Remove appendages from the mask
            mask_without_appendages = mask.remove_appendages()

            # Show a plot for debugging
            if config.debug.overlap_after or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=mask_without_appendages), title="Overlapping mask after appendage removal")

        ## NEW: second appendage removal step
        if masks.overlap(source.background_mask, mask_without_appendages) and config.remove_appendages:

            # Show a plot for debugging
            if config.debug.overlap_before or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=mask_without_appendages), title="Overlapping mask before second appendage removal")

            if special: log.debug("second appendage removal step")

            # Do a second appendage removal
            mask_without_appendages = mask_without_appendages.remove_appendages(super=True)

            # Show a plot for debugging
            if config.debug.overlap_after or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=mask_without_appendages), title="Overlapping mask after second appendage removal")

    # Check if the mask hits the boundary of the cutout or overlaps with the background mask (depending on the configuration settings)
    if segmentation_expand_condition(mask_without_appendages, source.background_mask, config, special):

        # If expanding is enabled
        if config.expand:

            # Add a snapshot of the source to the track record for debugging
            if track_record is not None: track_record.append(copy.deepcopy(source))

            # Show a plot for debugging
            if config.debug.expand or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=masks.union(mask, source.background_mask)), title="Masked segment hits boundary [expansion level = " + str(expansion_level) + "]")

            # If the maximum expansion level has been reached, no source could be found
            if expansion_level >= config.max_expansion_level:

                if special:

                    log.debug("maximum expansion level reached (", expansion_level, ")")

                    # To visualize the case where maximum expansion has been reached
                    plotting.plot_box(np.ma.masked_array(source.cutout, mask=mask))

                # No source can be found
                return None

            else:

                # Calculate the expanded parameters
                ellipse *= config.expansion_factor
                expansion_level += 1

                if special: log.debug("expanding to level", expansion_level + 1, " (maximum level =", config.max_expansion_level)

                # Repeat the procedure for the expanded ellipse
                return find_source_segmentation(frame, ellipse, config, track_record=track_record, expansion_level=expansion_level, special=special)

        # If expanding is disabled, no source can be found
        else: return None

    else:

        if special: log.debug("center segment does not overlap with background mask")

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.success or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=masks.union(mask, source.background_mask)), title="Masked segment doesn't hit boundary")

        # -- Fill holes --

        if special: log.debug("fixing holes in segment mask")
        source.mask = mask.fill_holes()

        # Show a plot for debugging
        if config.debug.holes or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=source.mask), title="Removed holes")

        # -- Dilation --

        # Dilate the mask if requested
        #if config.dilate:
        if False:

            if special: log.debug("dilating the mask")

            source = source.zoom_out(1.5, frame, keep_original_mask=True)

            if special: source.plot(title="zoomed-out source before mask dilation")

            # Dilate the mask
            source.mask = source.mask.disk_dilation(radius=10, iterations=expansion_level)
            #mask = mask.dilated(connectivity=config.connectivity, iterations=config.iterations)

        if config.dilate:

            if special: log.debug("dilating the mask")

            source = source.zoom_out(config.dilation_factor, frame, keep_original_mask=True)

            mask_area = np.sum(source.mask)
            area_dilation_factor = config.dilation_factor ** 2.
            new_area = mask_area * area_dilation_factor

            ## Circular mask approximation

            #ellipse = find_contour(source.mask.astype(float), source.mask)
            #radius = ellipse.radius.norm

            mask_radius = math.sqrt(mask_area / math.pi)
            new_radius = math.sqrt(new_area / math.pi)

            kernel_radius = new_radius - mask_radius

            if special: log.debug("dilation disk radius:" + str(kernel_radius))

            source.mask = source.mask.disk_dilation(radius=kernel_radius)

        # Show a plot for debugging
        if config.debug.dilated or special: plotting.plot_box(np.ma.masked_array(source.cutout, mask=source.mask), title="Dilated mask")

        # -- Final source --

        # Inform the user
        if special:
            log.debug("source was found")
            log.debug("Final expansion level: " + str(expansion_level))

        # Return the source
        return source

# -----------------------------------------------------------------

def find_source_peaks(frame, ellipse, config, track_record=None, level=0, special=False):

    """
    This function ...
    :param frame:
    :param ellipse:
    :param config:
    :param track_record:
    :param level:
    :param special:
    :return:
    """

    # If the maximum or minimum level is reached, return without source
    if level < config.min_level or level > config.max_level: return None

    # Create a source object
    source = Detection.from_ellipse(frame, ellipse, config.background_outer_factor)

    # If the frame is zero in this box, continue to the next object
    if not np.any(source.cutout): return None

    # If the box is too small, skip this object
    if source.cutout.xsize < config.minimum_pixels or source.cutout.ysize < config.minimum_pixels: return None

    # If always subtract background is enabled
    if config.always_subtract_background: source.estimate_background(config.background_est_method, config.sigma_clip_background)

    # Check if a FWHM is defined for convolving the source cutout before looking for peaks
    if config.convolution_fwhm is not None:

        # Create a Gaussian convolution kernel and return it
        sigma = config.convolution_fwhm * statistics.fwhm_to_sigma
        kernel = Gaussian2DKernel(sigma)
        kernel.normalize() # to suppress warning

    # Else, set the kernel to None
    else: kernel = None

    # Find the location of peaks in the box (do not remove gradient yet for performance reasons)
    peaks = source.locate_peaks(config.sigma_level, kernel=kernel)

    # If no peaks could be detected, remove a potential background gradient from the box before detection
    if len(peaks) == 0 and not config.always_subtract_background:

        # Add a snapshot of the source to the track record for debugging
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.zero_peaks_before or special: source.plot(title="0 peaks, gradient background will be removed")

        # Estimate and subtract the background (remove the background gradient)
        source.estimate_background(config.background_est_method, config.sigma_clip_background)

        # Find the location of peaks in the box
        peaks = source.locate_peaks(config.sigma_level, kernel=kernel)

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
        ellipse *= config.scale_factor

        if special: log.debug("zooming in to find peak")

        # Find a source in the zoomed-out region
        return find_source_peaks(frame, ellipse, config, track_record=track_record, level=level+1, special=special)

    # If more than one source was detected
    elif len(peaks) > 1:

        # Add a snapshot of the source to the track record
        if track_record is not None: track_record.append(copy.deepcopy(source))

        # Show a plot for debugging
        if config.debug.more_peaks or special: source.plot(title="More peaks", peaks=peaks)

        # If the level was positive, no source can be found
        if level > 0: return None

        # Scale the ellipse in which to look for a source
        ellipse /= config.scale_factor

        # Find a source in the zoomed-in region
        return find_source_peaks(frame, ellipse, config, track_record=track_record, level=level-1, special=special)

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
        if not (np.isclose(x_peak, ellipse.center.x, atol=config.peak_offset_tolerance) and np.isclose(y_peak, ellipse.center.y, atol=config.peak_offset_tolerance)):

            # Show a plot for debugging
            if config.debug.off_center or special: source.plot(title="Peak and center position do not match")

            # No source was found
            return None

        # Else, return the source
        else: return source

# -----------------------------------------------------------------

def segmentation_expand_condition(mask, background_mask, config, special=False):

    """
    This function ...
    :return:
    """

    # If:
    #  - if overlapping with the source's background mask is allowed, check whether the mask does not hit the boundary.
    #     -> if it hits the boundary, (config.allow_overlap and not hits_boundary) evaluates to False, so that not (..) evaluates to True --> enter the if if it also overlaps (will be True because it also hits the boundary of the box)
    #     -> if it doesn't hit the boundary, this first part evaluates to False --> do not enter the if
    #  - if overlapping with the source's background mask is not allowed, config.allow_overlap = False -> (config.allow_overlap and ... ) = False --> not ( ... ) = True --> check right part to enter if ()
    if config.allow_overlap:

        if special: log.debug("Overlapping allowed; checking whether center segment hits boundary")

        if mask.hits_boundary(min_pixels=2):

            if special: log.debug("Center segment hits boundary of box: expand")
            return True

        else:

            if special: log.debug("Mask from center segment does not hit boundary: keep this mask for the saturation source")
            return False

    # Overlapping not allowed
    else:

        if special: log.debug("Overlapping no allowed; checking wether center segment overlaps source's background mask")

        if masks.overlap(background_mask, mask):

            if special: log.debug("Center segment overlaps background mask: expand")
            return True

        else:

            if special: log.debug("Mask from center segment does not overlap background mask: keep this mask for the saturation source")
            return False

# -----------------------------------------------------------------
