#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np
import copy

# Import astronomical modules
from astropy import units as u
from astropy.table import Table
from photutils import find_peaks
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
#from photutils import detect_threshold

# Import Astromagic modules
from .masks import Mask
from .box import Box
from ..tools import plotting
from ..core import masks
from ..tools import statistics
from .vector import Position

# *****************************************************************

class Source(object):

    """
    This class...
    """

    def __init__(self, frame, center, radius, angle, outer_factor):

        """
        The constructor ...
        """

        # Set attributes
        self.center = center
        self.radius = radius
        self.angle = angle

        # Create the cutout and background
        outer_radius = radius * outer_factor
        self.cutout = Box.from_ellipse(frame, center, radius, angle=angle)
        self.background = Box.from_ellipse(frame, center, outer_radius, angle=angle)

        # Calculate the relative coordinate of the center for the cutout and background boxes
        rel_center = self.cutout.rel_position(center)
        rel_center_background = self.background.rel_position(center)

        # Create masks for the background box that cover the given ellipse and
        #self.background_mask = masks.create_ellipse_mask(self.background.xsize, self.background.ysize, rel_center_background, radius, angle)
        self.background_mask = masks.create_annulus_mask(self.background.xsize, self.background.ysize, rel_center_background, radius, outer_radius, angle)

        # Create a mask that covers the given ellipse in the cutout box (can be overwritten by a mask for the center segment)
        self.mask = masks.create_ellipse_mask(self.cutout.xsize, self.cutout.ysize, rel_center, radius, angle)

        # Set subtracted box to None
        self.estimated_background = None
        self.estimated_background_cutout = None
        self.subtracted = None
        self.removed = None

        # Set peak position to None initially
        self.peak = None

    # *****************************************************************

    @property
    def is_subtracted(self):

        """
        This function ...
        :return:
        """

        return self.subtracted is not None

    # *****************************************************************

    @property
    def is_estimated(self):

        """
        This function ...
        :return:
        """

        return self.estimated_background is not None

    # *****************************************************************

    @property
    def flux(self):

        """
        This function ...
        :return:
        """

        # Calculate the relative coordinate of the center in the background box
        #rel_center = self.background.rel_position(self.center)

        # Method with photutils:
        #aperture = CircularAperture((x_center, y_center), r=x_radius)
        #phot_table = aperture_photometry(self.frames[frame_name], aperture, mask=np.isnan(self.frames[frame_name]))
        #flux = phot_table[0]["aperture_sum"]

        # Sum the pixel values

        #plotting.plot_box(np.ma.masked_array(np.asarray(self.subtracted), mask=self.cutout_mask.inverse()))

        return np.ma.sum(np.ma.masked_array(np.asarray(self.subtracted), mask=self.mask.inverse()))

    # *****************************************************************

    def estimate_background(self, method, sigma_clip=True, sigma_level=3.0):

        """
        This function ...
        :return:
        """

        # Perform sigma-clipping on the background if requested
        if sigma_clip: mask = statistics.sigma_clip_mask(self.background, sigma_level=sigma_level, mask=self.background_mask)
        else: mask = self.background_mask

        if method == "polynomial":

            # Estimate the background
            self.estimated_background = self.background.fit_polynomial(3, mask=mask)

        elif method == "interpolation":

            # Interpolate over the mask
            self.estimated_background = self.background.interpolate(mask)

        elif method == "mean":

            mean = np.ma.mean(np.ma.masked_array(self.background, mask=mask))
            self.estimated_background = self.background.fill(mean)

        elif method == "median":

            median = np.ma.median(np.ma.masked_array(self.background, mask=mask))
            self.estimated_background = self.background.fill(median)

        else: raise ValueError("Unknown background estimation method")

        # Create an estimated background box for the cutout
        self.estimated_background_cutout = self.estimated_background.box_like(self.cutout)

    # *****************************************************************

    def subtract_background(self):

        """
        This function ...
        :return:
        """

        # Subtract the background from the box
        self.subtracted = self.cutout - self.estimated_background_cutout

    # *****************************************************************

    def find_center_segment(self, threshold_sigmas, kernel_fwhm, kernel_size, min_pixels=5):

        """
        This function ...
        :return:
        """

        # Calculate threshold for segmentation
        mean, median, stddev = statistics.sigma_clipped_statistics(self.background, mask=self.background_mask)
        threshold = mean + stddev * threshold_sigmas

        #if threshold is None: threshold = detect_threshold(data, snr=signal_to_noise) #snr=2.0

        # Create a kernel
        sigma = kernel_fwhm * gaussian_fwhm_to_sigma
        kernel = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size)

        # Perform the segmentation
        segments = detect_sources(self.cutout, threshold, npixels=min_pixels, filter_kernel=kernel)

        # Get the label of the center segment
        rel_center = self.cutout.rel_position(self.center)
        label = segments[rel_center.y, rel_center.x]

        # If the center pixel is identified as being part of the background, create an empty mask (the center does not
        # correspond to a segment)
        if label == 0: self.mask = Mask(np.zeros_like(self.cutout, dtype=bool))
        else:

            # Create a mask of the center segment
            self.mask = Mask((segments == label))

    # *****************************************************************

    def locate_peaks(self, threshold_sigmas):

        """
        This function ...
        :return:
        """

        # If a subtracted box is present, use it to locate the peaks
        box = self.subtracted if self.subtracted is not None else self.cutout

        # Calculate the sigma-clipped statistics of the frame and find the peaks
        mean, median, stddev = statistics.sigma_clipped_statistics(box, sigma=3.0)
        threshold = median + (threshold_sigmas * stddev)

        # Find peaks
        peaks = find_peaks(box, threshold, box_size=5)

        # For some reason, once in a while, an ordinary list comes out of the find_peaks routine instead of an
        # Astropy Table instance. We assume we need an empty table in this case
        if type(peaks) is list: peaks = Table([[], []], names=('x_peak', 'y_peak'))

        # Initialize a list to contain the peak positions
        positions = []

        # Loop over the peaks
        for peak in peaks:

            # Calculate the absolute x and y coordinate of the peak
            x = peak['x_peak'] + self.cutout.x_min
            y = peak['y_peak'] + self.cutout.y_min

            # Add the coordinates to the positions list
            positions.append(Position(x=x,y=y))

        # If exactly one peak was found, set the self.peak attribute accordingly
        if len(positions) == 1: self.peak = positions[0]

        # Return the list of peak positions
        return positions

    # *****************************************************************

    def zoom(self, factor):

        """
        This function ...
        :return:
        """

        # Create a copy of this object
        source = copy.deepcopy(self)

        # Zoom in on the cutout
        source.cutout = self.cutout.zoom(self.center, factor)

        # Zoom in on the background
        source.background = self.background.zoom(self.center, factor)

        # Set derived properties to None
        source.estimated_background = None
        source.estimated_background_cutout = None
        source.subtracted = None
        source.removed = None

        # Return the new source
        return source

    # *****************************************************************

    def plot(self, title=None, peaks=None):

        """
        This function ...
        :return:
        """

        if peaks is not None:

            x_positions = []
            y_positions = []

            # Loop over all peaks
            for peak in peaks:

                rel_peak = self.cutout.rel_position(peak)
                x_positions.append(rel_peak.x)
                y_positions.append(rel_peak.y)

            peak_coordinates = [x_positions, y_positions]

        elif self.peak is not None:

            rel_peak = self.cutout.rel_position(self.peak)
            peak_coordinates = [[rel_peak.x], [rel_peak.y]]

        else: peak_coordinates = None

        # If the background has been estimated for this source
        if self.estimated_background_cutout is not None:

            # Do the plotting
            plotting.plot_source(self.background, self.background_mask, self.estimated_background, self.cutout,
                                 self.estimated_background_cutout, self.mask, peaks=peak_coordinates, title=title)

        # Else, we just have a background and cutout box
        else:

            # Do the plotting
            plotting.plot_background_center(self.background, self.background_mask, self.cutout, peaks=peak_coordinates, title=title)

# *****************************************************************