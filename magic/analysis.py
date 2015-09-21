#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

# Import image modules
import fitting
import plotting
import regions
import statistics
import masks
from tools import coordinates, cropping, interpolation

# Import astronomical modules
from astropy.table import Table
from photutils import find_peaks
from find_galaxy import find_galaxy
from astropy import log
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
from photutils import daofind
from astropy.stats import sigma_clipped_stats
from photutils import detect_threshold

# *****************************************************************

def find_segments(data, kernel_fwhm=2.0, kernel_size=3.0, threshold=None, signal_to_noise=2.0):

    """
    This function ...
    :param data:
    :param shape:
    :return:
    """

    if threshold is None: threshold = detect_threshold(data, snr=signal_to_noise)

    sigma = kernel_fwhm * gaussian_fwhm_to_sigma    # FWHM = 2.
    kernel = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size)

    segments = detect_sources(data, threshold, npixels=5, filter_kernel=kernel)

    return segments

# *****************************************************************

def find_sources_in_region(data, region, model_name, detection_method, plot=False, plot_custom=[False, False, False, False]):

    """
    This function searches for sources within a given frame and a region of ellipses
    :param data:
    :param region:
    :param method: "peaks", "segmentation", "dao", "iraf"
    :param plot:
    :return:
    """

    # Initialize an empty list of sources
    sources = []

    # Inform the user
    log.info("Looking for sources...")

    # Loop over all objects
    for shape in region:

        # Find a source
        source = find_source_in_shape(data, shape, model_name, detection_method, level=0, plot=plot, plot_custom=plot_custom)

        # If a source was found, add it to the list of sources
        if source is not None: sources.append(source)

    # Return the list of sources
    return sources

# *****************************************************************

def find_source_in_shape(data, shape, model_name, detection_method, level, min_level=-2, max_level=2, aid_detection_by_removing_gradient=True, peak_offset_tolerance=3.0, plot=False, plot_custom=[False, False, False, False]):

    """
    This function ...
    :param data:
    :param shape:
    :param detection_method:
    :param plot:
    :return:
    """

    if level < min_level or level > max_level: return None

    # Find the parameters of this ellipse (or circle)
    x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

    # Cut out a box of the primary image around the star
    box, x_min, x_max, y_min, y_max = cropping.crop(data, x_center, y_center, x_radius, y_radius)

    # If the frame is zero in this box, continue to the next object
    if not np.any(box): return None

    # If the box is too small, skip this object
    if box.shape[0] < 5 or box.shape[1] < 5: return None

    # Find the location of sources (peaks)
    # TODO: find peaks only within the shape, instead of the entire box cut around the shape
    positions = find_source_positions(box, detection_method, x_shift=x_min, y_shift=y_min, plot=False)

    # If no sources could be detected, remove a potential background gradient from the box before detection
    if len(positions) == 0 and aid_detection_by_removing_gradient:

        if plot_custom[0]: plotting.plot_box(box, title="0 peaks")

        # Create a central mask that covers the expected source
        central_radius = 5.0
        central_mask = np.zeros_like(box, dtype=bool)
        central_mask[int(round(box.shape[0]/2.0-central_radius)):int(round(box.shape[0]/2.0+central_radius)),int(round(box.shape[1]/2.0-central_radius)):int(round(box.shape[1]/2.0+central_radius))] = True

        if np.all(central_mask): return None

        # Fit a polynomial to the background
        poly = fitting.fit_polynomial(box, 3, mask=central_mask)
        polynomial_box = fitting.evaluate_model(poly, 0, box.shape[1], 0, box.shape[0])

        if plot_custom[0]: plotting.plot_difference_model(box, poly)

        # Detect source positions with the background gradient removed
        positions = find_source_positions(box-polynomial_box, detection_method, x_shift=x_min, y_shift=y_min, plot=False)

        if plot_custom[0]: plotting.plot_peaks(box, positions['x_peak']-x_min, positions['y_peak']-y_min, title=str(len(positions))+" peak(s) found after removing gradient background")

    # If no sources were detected
    if len(positions) == 0:

        if plot_custom[1]: plotting.plot_box(box, title="0 peaks")

        if level < 0: return None

        new_shape = regions.scale_circle(shape, 2.0)
        return find_source_in_shape(data, new_shape, model_name, detection_method, level=level+1, min_level=min_level, max_level=max_level, plot=plot, plot_custom=[plot_custom[1],plot_custom[1],plot_custom[1],plot_custom[1]])

    # If one source was detected
    elif len(positions) == 1:

        if plot_custom[2]: plotting.plot_peaks(box, positions['x_peak']-x_min, positions['y_peak']-y_min, title="1 peak")

        # Create the mask to estimate the background in the box
        background_inner_factor = 1.0
        background_outer_factor = 1.5
        fit_factor = 1.0
        sigma_clip_background = True

        # Get the x and y coordinate of the peak
        x_peak = positions['x_peak'][0]
        y_peak = positions['y_peak'][0]

        # TODO: check that peak position corresponds with center of shape
        if not (np.isclose(x_peak, x_center, atol=peak_offset_tolerance) and np.isclose(y_peak, y_center, atol=peak_offset_tolerance)):

            log.warning("Peak position and shape center do not match")

            #plotting.plot_peaks(box, positions['x_peak']-x_min, positions['y_peak']-y_min, title="Peak and shape do not match")
            return None

        # Try to fit a Gaussian model to the data within the shape
        source = make_gaussian_star_model(shape, data, x_peak, y_peak, background_inner_factor, background_outer_factor,
                                          fit_factor, "polynomial", peak_offset_tolerance,
                                          sigma_clip_background=sigma_clip_background, upsample_factor=1.0,
                                          use_shape_or_peak="peak", model_name=model_name, plot=plot)

        return source

    # If more than one source was detected
    elif len(positions) > 1:

        if plot_custom[3]: plotting.plot_peaks(box, positions['x_peak']-x_min, positions['y_peak']-y_min, title="more peaks")

        if level > 0: return None

        new_shape = regions.scale_circle(shape, 0.5)
        return find_source_in_shape(data, new_shape, model_name, detection_method, level=level-1, min_level=min_level, max_level=max_level, plot=plot, plot_custom=[plot_custom[3],plot_custom[3],plot_custom[3],plot_custom[3]])

# *****************************************************************

def find_source_positions(data, detection_method, x_shift=0.0, y_shift=0.0, plot=False):

    """
    This function ...
    :param data:
    :param detection_method:
    :param x_shift:
    :param y_shift:
    :param plot:
    :return:
    """

    if detection_method == "peaks": peaks, median = locate_sources_peaks(data, plot=plot)
    elif detection_method == "segmentation": peaks, median = locate_sources_segmentation(data)
    elif detection_method == "dao": peaks, median = locate_sources_daofind(data, plot=plot)
    elif detection_method == "iraf": peaks, median = locate_sources_iraf(data)
    elif detection_method == "all":

        peaks, median = locate_sources_peaks(data)
        peaks1, median1 = locate_sources_daofind(data)

        peaks += peaks1
        median = np.mean(median, median1)

    #elif detection_method is None: pass
    else: raise ValueError("Unknown source extraction method")

    peaks['x_peak'] += x_shift
    peaks['y_peak'] += y_shift

    return peaks

# *****************************************************************

def locate_sources_peaks(data, threshold_sigmas=7.0, plot=False):

    """
    This function looks for peaks in the given frame
    :param data:
    :param sigma:
    :return:
    """

    # Calculate the sigma-clipped statistics of the frame and find the peaks
    mean, median, stddev = statistics.sigma_clipped_statistics(data, sigma=3.0)
    threshold = median + (threshold_sigmas * stddev)
    peaks = find_peaks(data, threshold, box_size=5)

    # For some reason, once in a while, an ordinary list comes out of the find_peaks routine instead of an Astropy Table instance. We assume we need an empty table in this case
    if type(peaks) is list: peaks = Table([[], []], names=('x_peak', 'y_peak'))

    # If requested, make a plot with the source(s) indicated
    if plot and len(peaks) > 0: plotting.plot_peaks(data, peaks['x_peak'], peaks['y_peak'])

    # Return the list of peaks
    return peaks, median

# *****************************************************************

def locate_sources_daofind(data, threshold_sigmas=5.0, plot=False):

    """
    This function ...
    :param data:
    :return:
    """

    # Calculate the sigma-clipped statistics of the data
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    result_table = daofind(data - median, fwhm=3.0, threshold=threshold_sigmas*std)

    result_table.rename_column('xcentroid', 'x_peak')
    result_table.rename_column('ycentroid', 'y_peak')

    # If requested, make a plot with the source(s) indicated
    if plot: plotting.plot_peaks(data, result_table['x_peak'], result_table['y_peak'], radius=4.0)

    # Return the list of source positions
    return result_table, median

# *****************************************************************

def locate_sources_iraf(data):

    """
    This function ...
    :param data:
    :return:
    """

    pass

# *****************************************************************

def locate_sources_segmentation(data):

    """
    This function ...
    :param data:
    :return:
    """

    sigma = 4.0 * gaussian_fwhm_to_sigma

    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)

    frame_name = self.frames.get_selected(require_single=True)

    from photutils import detect_threshold

    total_mask = self.combine_masks(return_mask=True)

    threshold = detect_threshold(self.frames[frame_name].data, snr=3., mask=total_mask)

    segm = detect_sources(self.frames[frame_name].data, threshold, npixels=5, filter_kernel=kernel)

# *****************************************************************

def find_sources_peaks(data, peaks, plot=None, x_shift=0.0, y_shift=0.0):

    """
    This function searches for sources in the given frame, based on a list of peaks detected in this frame
    :param data: the frame
    :param peaks: the list of peaks
    :param plot: optional, specific plotting flags
    :param x_shift: the shift in the x direction that is required to ...
    :param y_shift: the shift in the y direction that is required to ...
    :return:
    """

    # Initialize an empty list of sources
    sources = []

    if plot is None: plot=[False, False, False, False]

    # No peaks are found above the threshold
    if len(peaks) == 0:

        # Plot the data
        if plot[0]: plotting.plot_box(data)

        # Look for a star
        sources += find_sources_nopeak(data, x_shift, y_shift, plot=plot[0])

    # Exactly one peak is found
    elif len(peaks) == 1:

        x_peak = peaks['x_peak'][0]
        y_peak = peaks['y_peak'][0]

        # Look for a star corresponding to this peak
        sources += find_sources_1peak(data, x_peak, y_peak, x_shift, y_shift, plot=plot[1])

    # Two peaks are found in the box
    elif len(peaks) == 2:

        x_peak1 = peaks['x_peak'][0]
        y_peak1 = peaks['y_peak'][0]

        x_peak2 = peaks['x_peak'][1]
        y_peak2 = peaks['y_peak'][1]

        # Look for stars corresponding to the two found peaks
        sources += find_sources_2peaks(data, x_peak1, y_peak1, x_peak2, y_peak2, x_shift, y_shift, plot=plot[2])

    # More than two peaks are found
    elif len(peaks) > 2:

        # Look for multiple stars
        sources += find_sources_multiplepeaks(data, peaks['x_peak'], peaks['y_peak'], x_shift, y_shift, plot=plot[3])

    # Return the list of sources
    return sources

# *****************************************************************

def find_sources_nopeak(box, x_min, y_min, plot=False):

    """
    This function ...
    :param box:
    :param x_min:
    :param y_min:
    :param plot:
    :return:
    """

    # Initiate an empty list of stars
    sources = []

    # If the box is too small, don't bother looking for stars (anymore)
    if box.shape[0] < 5 or box.shape[1] < 5: return []

    # Fit a 2D Gaussian to the data
    model = fitting.fit_2D_Gaussian(box)

    # Get the parameters of the Gaussian function
    amplitude = model.amplitude
    x_mean = model.x_mean.value
    y_mean = model.y_mean.value

    # Skip non-stars (negative amplitudes)
    if amplitude < 0:

        log.warning("Source profile has negative amplitude")
        return []

    # If the center of the Gaussian falls out of the square, skip this star
    if round(x_mean) < 0 or round(x_mean) >= box.shape[0] - 1 or round(y_mean) < 0 or round(y_mean) >= box.shape[1] - 1:

        log.warning("Center of source profile lies outside of the region where it was expected")
        return []

    # Determine the coordinates of the center pixel of the box
    x_center = 0.5*box.shape[1]
    y_center = 0.5*box.shape[0]

    # Plot if requested
    if plot: plotting.plot_peak_model(box, x_center, y_center, model)

    # Change the coordinates of the model to absolute coordinates
    model.x_mean.value = model.x_mean.value + x_min
    model.y_mean.value = model.y_mean.value + y_min

    # Add the model for this source to the list of sources
    sources.append(model)

    # Return the list of sources (containing only one source in this case)
    return sources

# *****************************************************************

def find_sources_1peak(box, x_peak, y_peak, x_min, y_min, model_name="Gaussian", initial_radius=None, center_offset_tolerance=None, plot=False):

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

    # Initiate an empty list of stars
    sources = []

    # If the box is too small, don't bother looking for stars (anymore)
    if box.shape[0] < 5 or box.shape[1] < 5: return []

    #if initial_radius is None: initial_radius = 0.2*box.shape[0]

    # Fit a 2D Gaussian to the data
    if model_name == "Gaussian":

        model = fitting.fit_2D_Gaussian(box, center=(x_peak, y_peak), radius=initial_radius)
        x_mean = model.x_mean.value
        y_mean = model.y_mean.value

    elif model_name == "Airy":

        model = fitting.fit_2D_Airy(box, center=(x_peak, y_peak))
        x_mean = model.x_0.value
        y_mean = model.y_0.value

    else: raise ValueError("Other models are not supported yet")

    amplitude = model.amplitude

    # Skip non-stars (negative amplitudes)
    if amplitude < 0:

        log.warning("Source profile has negative amplitude")
        return []

    # Check if the center of the Gaussian corresponds to the position of the detected peak
    diff_x = x_mean - x_peak
    diff_y = y_mean - y_peak

    distance = np.sqrt(diff_x**2+diff_y**2)

    if distance > center_offset_tolerance:

        log.warning("Center of source profile lies outside of the region where it was expected (distance = "+str(distance))

        if plot: plotting.plot_peak_model(box, x_peak, y_peak, model)

        # Calculate the size of the smaller box
        smaller_box_ysize = int(round(box.shape[0] / 4.0))
        smaller_box_xsize = int(round(box.shape[1] / 4.0))

        # Create a smaller box
        smaller_box, smaller_xmin, smaller_xmax, smaller_ymin, smaller_ymax = cropping.crop(box, x_peak, y_peak, smaller_box_xsize, smaller_box_ysize)

        # Calculate x and y delta
        x_delta = x_min + smaller_xmin
        y_delta = y_min + smaller_ymin

        # Find relative coordinate of peak
        rel_xpeak, rel_ypeak = coordinates.relative_coordinate(x_peak, y_peak, smaller_xmin, smaller_ymin)

        # Try again (iterative procedure of zooming in, stops if the size of the box becomes too small)
        sources += find_sources_1peak(smaller_box, rel_xpeak, rel_ypeak, x_delta, y_delta, initial_radius=initial_radius, center_offset_tolerance=center_offset_tolerance, plot=plot)

    # Potential star detected!
    else:

        # Plot if requested
        if plot: plotting.plot_peak_model(box, x_peak, y_peak, model)

        # Change the coordinates of the model to absolute coordinates
        if model_name == "Gaussian":

            model.x_mean.value = model.x_mean.value + x_min
            model.y_mean.value = model.y_mean.value + y_min

        elif model_name == "Airy":

            model.x_0.value = model.x_0.value + x_min
            model.y_0.value = model.y_0.value + y_min

        # Add the model for this source to the list of sources
        sources.append(model)

    # Return the list of sources (containing only one source in this case)
    return sources

# *****************************************************************

def find_sources_2peaks(box, x_peak1, y_peak1, x_peak2, y_peak2, x_min, y_min, plot=False):

    """
    This function searches for sources ...
    :param box:
    :param x_peak1:
    :param y_peak1:
    :param x_peak2:
    :param y_peak2:
    :param x_min:
    :param y_min:
    :param plot:
    :return:
    """

    # Initiate an empty list of stars
    sources = []

    # Calculate the distance between the two peaks
    diff_x = x_peak1 - x_peak2
    diff_y = y_peak1 - y_peak2
    distance = np.sqrt(diff_x**2 + diff_y**2)

    # Calculate the midpoint between the two peaks
    mid_x = 0.5*(x_peak1 + x_peak2)
    mid_y = 0.5*(y_peak1 + y_peak2)

    # The peaks are probably part of the same object
    # So, perform a fit with a Gaussian
    if distance < 4:

        # Do the fit
        model = fitting.fit_2D_Gaussian(box, center=(mid_x, mid_y))

        # Plot if requested
        if plot: plotting.plot_peaks_models(box, [x_peak1, x_peak2], [y_peak1, y_peak2], [model])

        # Change the coordinates of the model to absolute coordinates
        model.x_mean.value = model.x_mean.value + x_min
        model.y_mean.value = model.y_mean.value + y_min

        # Add the model for this source to the list of sources
        sources.append(model)

    # Split the box
    else:

        # Calculate the absolute distance between the two peaks
        delta_x = np.abs(diff_x)
        delta_y = np.abs(diff_y)

        # Cut along the axis where the distance is largest
        if delta_x > delta_y:

            # Create the boxes
            smaller_box_1 = box[:][0:mid_x]
            smaller_box_2 = box[:][mid_x:]

            # Determine the relative positions of the found peaks
            rel_peak1_x, rel_peak1_y = coordinates.relative_coordinate(x_peak1, y_peak1, 0, 0)
            rel_peak2_x, rel_peak2_y = coordinates.relative_coordinate(x_peak2, y_peak2, mid_x, 0)

            # Determine the translation needed for the second cut
            rel_xmin2 = x_min + mid_x
            rel_ymin2 = y_min

        else:

            # Create the boxes
            smaller_box_1 = box[0:mid_y][:]
            smaller_box_2 = box[mid_y:][:]

            # Determine the relative positions of the found peaks
            rel_peak1_x, rel_peak1_y = coordinates.relative_coordinate(x_peak1, y_peak1, 0, 0)
            rel_peak2_x, rel_peak2_y = coordinates.relative_coordinate(x_peak2, y_peak2, 0, mid_y)

            # Determine the translation needed for the second cut
            rel_xmin2 = x_min
            rel_ymin2 = y_min + mid_y

        # Find stars in both of the cuts
        sources += find_sources_1peak(smaller_box_1, rel_peak1_x, rel_peak1_y, x_min, y_min, plot=plot)
        sources += find_sources_1peak(smaller_box_2, rel_peak2_x, rel_peak2_y, rel_xmin2, rel_ymin2, plot=plot)

    # Return the list of sources
    return sources

# *****************************************************************

def find_sources_multiplepeaks(box, x_peaks, y_peaks, x_min, y_min, plot=False):

    """
    This function searches for sources ...
    :param box:
    :param x_peaks:
    :param y_peaks:
    :param x_min:
    :param y_min:
    :param plot:
    :return:
    """

    # Initialize an empty list of stars
    sources = []

    # Calculate the new x and y size of the smaller boxes
    smaller_ysize = int(round(box.shape[0] / 2.0))
    smaller_xsize = int(round(box.shape[1] / 2.0))

    # Define the boundaries of 4 new boxes
    boundaries_list = []
    boundaries_list.append((0,smaller_ysize, 0,smaller_xsize))
    boundaries_list.append((0,smaller_ysize, smaller_xsize,box.shape[1]))
    boundaries_list.append((smaller_ysize,box.shape[0], 0,smaller_xsize))
    boundaries_list.append((smaller_ysize,box.shape[0], smaller_xsize,box.shape[1]))

    # For each new box
    for boundaries in boundaries_list:

        # Create the box
        smaller_box = box[boundaries[0]:boundaries[1],boundaries[2]:boundaries[3]]

        # Initialize empty lists to contain the coordinates of the peaks that fall within the current box
        x_peaks_inside = []
        y_peaks_inside = []

        # Calculate the translation needed for this smaller box
        abs_x_min = x_min + boundaries[2]
        abs_y_min = y_min + boundaries[0]

        # Check how many peaks fall within this box
        for x_peak, y_peak in zip(x_peaks, y_peaks):

            if x_peak >= boundaries[2] and x_peak < boundaries[3] and y_peak >= boundaries[0] and y_peak < boundaries[1]:

                x_peaks_inside.append(x_peak)
                y_peaks_inside.append(y_peak)

        # No peaks in this box
        if len(x_peaks_inside) == 0: continue

        # One peak in this box
        if len(x_peaks_inside) == 1:

            # Determine the relative coordinate of the peak
            x_peak_rel, y_peak_rel = coordinates.relative_coordinate(x_peaks_inside[0], y_peaks_inside[0], boundaries[2], boundaries[0])

            # Find the source
            sources += find_sources_1peak(smaller_box, x_peak_rel, y_peak_rel, abs_x_min, abs_y_min, plot=plot)

        # Two peaks in this box
        elif len(x_peaks_inside) == 2:

            # Determine the relative coordinates of the peaks
            x_peak1_rel, y_peak1_rel = coordinates.relative_coordinate(x_peaks_inside[0], y_peaks_inside[0], boundaries[2], boundaries[0])
            x_peak2_rel, y_peak2_rel = coordinates.relative_coordinate(x_peaks_inside[1], y_peaks_inside[1], boundaries[2], boundaries[0])

            # Find sources
            sources += find_sources_2peaks(smaller_box, x_peak1_rel, y_peak1_rel, x_peak2_rel, y_peak2_rel, abs_x_min, abs_y_min, plot=plot)

        # More than 2 peaks
        else:

            # Initialize empty lists for the relative coordinates of the peaks within this smaller box
            x_peaks_rel = []
            y_peaks_rel = []

            # Calculate these relative coordinates
            for x_peak_inside, y_peak_inside in zip(x_peaks_inside, y_peaks_inside):

                x_peak_rel, y_peak_rel = coordinates.relative_coordinate(x_peak_inside, y_peak_inside, boundaries[2], boundaries[0])
                x_peaks_rel.append(x_peak_rel)
                y_peaks_rel.append(y_peaks_rel)

            # Find sources
            sources += find_sources_multiplepeaks(smaller_box, x_peaks_rel, y_peaks_rel, abs_x_min, abs_y_min, plot=plot)

    # Return the parameters of the found stars
    return sources

# *****************************************************************

def remove_duplicate_sources(sources):

    """
    This function ...
    :param sources:
    :return:
    """

    # Inform the user
    log.info("Removing duplicates...")

    # Initialize a list that only contains unique sources
    unique_sources = []

    # Loop over all sources, adding them to the list of unique sources if not yet present
    for source in sources:

        # Check if another source on the same location is already present in the unique_sources list
        for unique_source in unique_sources:

            distance = coordinates.distance_models(source, unique_source)
            if distance < 1.0: break

        # No break statement was encountered; add the source to the unique_sources list
        else: unique_sources.append(source)

    # Return the list of unique sources
    return unique_sources

# *****************************************************************

def estimate_background(data, mask, interpolate=True, sigma_clip=True):

    """
    This function ...
    :param data:
    :param mask:
    :return:
    """

    # Sigma clipping
    if sigma_clip: mask = statistics.sigma_clip_mask(data, sigma=3.0, mask=mask)

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

def make_gaussian_star_model(shape, data, x_peak, y_peak, background_inner_factor, background_outer_factor, fit_factor,
                             background_est_method, center_offset_tolerance, sigma_clip_background=True,
                             upsample_factor=1.0, use_shape_or_peak="peak", model_name="Gaussian", plot=False):

    """
    This function ...
    :param shape:
    :param data:
    :param background_outer_factor:
    :param fit_factor:
    :param background_est_method: 'median', 'interpolation' or 'polynomial'
    :param sigma_clip_background:
    :param plot:
    :return:
    """

    # Find the parameters of this ellipse (or circle)
    x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

    # Create a mask that
    annulus_mask = masks.annulus_around(shape, background_inner_factor, background_outer_factor, data.shape[1], data.shape[0])

    # Set the radii for cutting out the background box
    radius = 0.5*(x_radius + y_radius)
    x_radius_outer = background_outer_factor*x_radius
    y_radius_outer = background_outer_factor*y_radius

    # Cut out the background and the background mask
    background, x_min_back, x_max_back, y_min_back, y_max_back = cropping.crop(data, x_center, y_center, x_radius_outer, y_radius_outer)
    background_mask = cropping.crop_check(annulus_mask, x_min_back, x_max_back, y_min_back, y_max_back)

    # Estimate the background
    if background_est_method == "median":

        pass

    elif background_est_method == "mean":

        pass

    if background_est_method == "polynomial":

        background_mask_beforeclipping = np.copy(background_mask)
        poly, background_mask = fitting.fit_polynomial(background, 3, mask=background_mask, sigma_clip_background=sigma_clip_background)
        est_background = fitting.evaluate_model(poly, 0, background.shape[1], 0, background.shape[0])

    elif background_est_method == "interpolation":

        background_mask_beforeclipping = np.copy(background_mask)
        est_background, background_mask = estimate_background(background, background_mask, interpolate=True, sigma_clip=sigma_clip_background)

    else: raise ValueError("Unknown background estimation method")

    # Cut out a box of the primary image around the star
    box, x_min, x_max, y_min, y_max = cropping.crop(data, x_center, y_center, fit_factor*x_radius, fit_factor*y_radius)

    # Crop the interpolated background to the frame of the box
    box_background = cropping.crop_check(est_background, x_min-x_min_back, x_max-x_min_back, y_min-y_min_back, y_max-y_min_back)

    # Subtract background from the box
    box = box - box_background

    # Find source
    if use_shape_or_peak == "shape":
        x_center_rel, y_center_rel, = coordinates.relative_coordinate(x_center, y_center, x_min, y_min)
    elif use_shape_or_peak == "peak":
        x_center_rel, y_center_rel = coordinates.relative_coordinate(x_peak, y_peak, x_min, y_min)
    else:
        raise ValueError("Invalid option (should be 'shape' or 'peak')")

    sources = find_sources_1peak(box, x_center_rel, y_center_rel, x_min, y_min, model_name=model_name, center_offset_tolerance=center_offset_tolerance, plot=plot)

    # Search again, with an Airy model
    if not sources: sources = find_sources_1peak(box, x_center_rel, y_center_rel, x_min, y_min, model_name="Airy", center_offset_tolerance=center_offset_tolerance, plot=plot)

    # If a source was found, return it
    if sources:

        source = sources[0]

        if plot:

            evaluated_model = fitting.evaluate_model(source, x_min, x_max, y_min, y_max)
            plotting.plot_star_model(background=np.ma.masked_array(background,mask=background_mask_beforeclipping),
                                     background_clipped=np.ma.masked_array(background,mask=background_mask),
                                     est_background=est_background,
                                     star=np.ma.masked_array(box, mask=None),
                                     est_background_star= box_background,
                                     fitted_star=evaluated_model)

        return source

    else:

        if plot:

            plotting.plot_background_subtraction(background=np.ma.masked_array(background,mask=background_mask_beforeclipping),
                                                 background_clipped=np.ma.masked_array(background,mask=background_mask),
                                                 est_background=est_background,
                                                 star=np.ma.masked_array(box, mask=None), est_background_star= box_background)

        return None

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

def crop_and_mask_for_background(data, shape, inner_factor, outer_factor):

    """
    This function ...
    :param data:
    :param shape:
    :param inner_factor:
    :param outer_factor:
    :return:
    """

    # Find the parameters of this ellipse (or circle)
    x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

    # Set the radii for cutting out the background box
    x_radius_outer = outer_factor*x_radius
    y_radius_outer = outer_factor*y_radius

    # Cut out the background and the background mask
    background, x_min, x_max, y_min, y_max = cropping.crop(data, x_center, y_center, x_radius_outer, y_radius_outer)

    x_center_rel = x_center - x_min
    y_center_rel = y_center - y_min

    # Create the mask
    mask = masks.create_disk_mask(background.shape[1], background.shape[0], x_center_rel, y_center_rel, x_radius)

    try:

        background[1,2]

    except IndexError:

        print outer_factor
        print x_radius_outer
        print y_radius_outer
        print background.shape
        print x_min
        print x_max
        print y_min
        print y_max

    return np.ma.masked_array(background, mask=mask), x_min, x_max, y_min, y_max

# *****************************************************************

def find_center_segment_in_shape(data, shape, kernel_fwhm, kernel_size, threshold_sigmas, expand=True, expansion_factor=1.5, expansion_level=1, max_expansion_level=4, plot=False):

    """
    This function ...
    :return:
    """

    # Get the parameters of this shape
    x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

    # Crop
    box, x_min, x_max, y_min, y_max = cropping.crop(data, x_center, y_center, x_radius, y_radius)

    background, x_min_back, x_max_back, y_min_back, y_max_back = crop_and_mask_for_background(data, shape, 1.0, 1.2)

    #plotting.plot_box(box_background, title="Box for background estimation")

    # Remove gradient
    poly = fitting.fit_polynomial(background.data, 3, mask=background.mask)
    polynomial = fitting.evaluate_model(poly, 0, background.shape[1], 0, background.shape[0])

    #plotting.plot_difference(background, polynomial)

    background = background - polynomial

    mean, median, stddev = statistics.sigma_clipped_statistics(background.data, mask=background.mask)

    #print mean, median, stddev

    threshold = mean + threshold_sigmas*stddev

    #print y_min-y_min_back >= 0, y_max-y_min_back >= 0, x_min-x_min_back >= 0, x_max-x_min_back >= 0

    evaluated_poly = polynomial[y_min-y_min_back:y_max-y_min_back, x_min-x_min_back:x_max-x_min_back]

    #print box.shape, evaluated_poly.shape

    plotting.plot_difference(box, evaluated_poly)

    #oldbox = box
    box = box - evaluated_poly

    #plotting.plot_difference(oldbox, box)

    # Find segments
    segments = find_segments(box, kernel_fwhm=kernel_fwhm, kernel_size=kernel_size, threshold=threshold)

    #label_im, nb_labels = ndimage.label(mask)

    label = segments[y_center-y_min, x_center-x_min]

    box_mask = (segments == label)

    # Test whether the mask reaches the boundaries
    hits_boundary = False
    for x in range(box_mask.shape[1]):

        if box_mask[0, x] or box_mask[box_mask.shape[0]-1, x]:

            # If this already happened with another pixel, break the loop
            if hits_boundary:

                hits_boundary = True
                break

            # If this is the first pixel for which this occurs, continue (one masked pixel on the edge is tolerated)
            else: hits_boundary = True

    for y in range(1, box_mask.shape[0]-1):

        if box_mask[y, 0] or box_mask[y, box_mask.shape[1]-1]:

            # If this already happened with another pixel, break the loop
            if hits_boundary:

                hits_boundary = True
                break

            # If this is the first pixel for which this occurs, continue (one masked pixel on the edge is tolerated)
            else: hits_boundary = True

    if hits_boundary and expand:

        if plot: plotting.plot_box(np.ma.masked_array(box, mask=box_mask), title="Masked segment hits boundary")

        if expansion_level == max_expansion_level: box_mask.fill(False)
        else:

            shape = regions.scale_circle(shape, expansion_factor)
            box_mask, x_min, x_max, y_min, y_max = find_center_segment_in_shape(data, shape, kernel_fwhm, kernel_size,
                                                                                threshold_sigmas, expand=expand,
                                                                                expansion_factor=expansion_factor,
                                                                                expansion_level=expansion_level+1,
                                                                                max_expansion_level=max_expansion_level,
                                                                                plot=plot)
    else:

        if plot: plotting.plot_box(np.ma.masked_array(box, mask=box_mask), title="Masked segment doesn't hit boundary")

    # Return the mask
    return box_mask, x_min, x_max, y_min, y_max

# *****************************************************************