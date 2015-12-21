#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.source Contains the Source class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import copy

# Import astronomical modules
from astropy.table import Table
from photutils import find_peaks
from photutils import detect_sources
from photutils import detect_threshold
from astropy.convolution import convolve, convolve_fft

# Import the relevant AstroMagic classes and modules
from .box import Box
from ..basics import Position, Mask
from ..tools import plotting, statistics

# -----------------------------------------------------------------

class Source(object):

    """
    This class...
    """

    def __init__(self, frame, center, radius, angle, factor):

        """
        The constructor ...
        """

        # Set attributes
        self.center = center
        self.radius = radius
        self.angle = angle
        self.factor = factor

        # Create the cutout box
        self.cutout = Box.from_ellipse(frame, center, self.radius*self.factor, self.angle)

        # Calculate the relative coordinate of the center for the cutout box
        rel_center = self.cutout.rel_position(center)

        # Create masks
        self.mask = Mask.from_ellipse(self.cutout.xsize, self.cutout.ysize, rel_center, radius, self.angle)

        # Set (estimated) background and removed to None
        self.background = None
        self.removed = None

        # Set peak position to None initially
        self.peak = None

    # -----------------------------------------------------------------

    @property
    def background_mask(self):

        """
        This function ...
        :return:
        """

        return self.mask.inverse()

    # -----------------------------------------------------------------

    @property
    def subtracted(self):

        """
        This function ...
        :return:
        """

        # Return the box with the background subtracted
        return self.cutout - self.background

    # -----------------------------------------------------------------

    @property
    def has_background(self):

        """
        This function ...
        :return:
        """

        return self.background is not None

    # -----------------------------------------------------------------

    @property
    def has_peak(self):

        """
        This function ...
        :return:
        """

        return self.peak is not None

    # -----------------------------------------------------------------

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

        return np.ma.sum(np.ma.masked_array(np.asarray(self.subtracted), mask=self.background_mask))

    # -----------------------------------------------------------------

    def estimate_background(self, method, sigma_clip=True, sigma_level=3.0):

        """
        This function ...
        :return:
        """

        # Perform sigma-clipping on the background if requested
        if sigma_clip: mask = statistics.sigma_clip_mask(self.cutout, sigma_level=sigma_level, mask=self.mask)
        else: mask = self.mask

        # Perform the interpolation
        self.background = self.cutout.interpolated(mask, method)

    # -----------------------------------------------------------------

    def find_center_segment(self, threshold_sigmas, kernel=None, min_pixels=5):

        """
        This function ...
        :return:
        """

        # If a subtracted box is present, use it to find the center segment
        box = self.subtracted if self.has_background else self.cutout

        if not np.all(self.mask):

            #print(box)
            #print(self.mask)

            # Calculate threshold for segmentation
            try:
                mean, median, stddev = statistics.sigma_clipped_statistics(box, mask=self.mask)
                threshold = mean + stddev * threshold_sigmas
            except TypeError:

                print(box)
                print(self.mask)

                print("not_nan=", np.sum(np.logical_not(np.isnan(box))))

                exit()

        else:

            threshold = detect_threshold(box, snr=2.0) #snr=2.0

        # Perform the segmentation
        segments = detect_sources(box, threshold, npixels=min_pixels, filter_kernel=kernel)

        # To plot the multiple segments that are detected
        #if segments.max() > 1: plotting.plot_box(np.ma.masked_array(box, mask=segments.astype(bool)))

        # Get the label of the center segment
        rel_center = self.cutout.rel_position(self.center)
        label = segments[rel_center.y, rel_center.x]
        #except IndexError:
        #plotting.plot_box(self.cutout)
        #plotting.plot_peak(self.cutout, rel_center.x, rel_center.y)

        # If the center pixel is identified as being part of the background, create an empty mask (the center does not
        # correspond to a segment)
        if label == 0: return Mask(np.zeros_like(self.cutout, dtype=bool))

        # Create a mask of the center segment
        else: return Mask((segments == label))

    # -----------------------------------------------------------------

    def locate_peaks(self, threshold_sigmas, kernel=None):

        """
        This function ...
        :return:
        """

        # If a subtracted box is present, use it to locate the peaks
        box = self.subtracted if self.has_background else self.cutout

        # Calculate the sigma-clipped statistics of the box
        mean, median, stddev = statistics.sigma_clipped_statistics(box, sigma=3.0)
        threshold = median + (threshold_sigmas * stddev)

        # Convolve the box with the given kernel, if any
        if kernel is not None: box = convolve_fft(box, kernel, normalize_kernel=True)

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
            positions.append(Position(x, y))

        # If exactly one peak was found, set the self.peak attribute accordingly
        if len(positions) == 1: self.peak = positions[0]

        # Return the list of peak positions
        return positions

    # -----------------------------------------------------------------

    def zoom(self, factor):

        """
        This function ...
        :return:
        """

        # Create a copy of this object
        source = copy.deepcopy(self)

        # Zoom in on the cutout
        source.cutout = self.cutout.zoom(self.center, factor)

        # Calculate the relative coordinate of the center for the cutout and background boxes
        rel_center = source.cutout.rel_position(self.center)

        # Decrease the radius
        source.radius = self.radius / factor

        # Create smaller mask
        source.mask = Mask.from_ellipse(source.cutout.xsize, source.cutout.ysize, rel_center, source.radius, source.angle)

        # Set other properties to None
        source.background = None
        source.removed = None

        # Return the new source
        return source

    # -----------------------------------------------------------------

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
        if self.has_background:

            # Do the plotting
            plotting.plot_source(self.cutout, self.mask, self.background, peaks=peak_coordinates, title=title)

        # Else, we just have a background and cutout box
        else:

            # Do the plotting
            plotting.plot_background_center(self.cutout, self.mask, peaks=peak_coordinates, title=title)

# -----------------------------------------------------------------
