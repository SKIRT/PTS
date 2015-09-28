#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy import units as u
from astropy.table import Table
from photutils import find_peaks

# Import Astromagic modules
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

    def __init__(self, frame, center, radius, angle, inner_factor, outer_factor):

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
        #rel_center = self.cutout.rel_position(center)
        rel_center_background = self.background.rel_position(center)

        # Create the masks that cover the source
        inner_radius = radius * inner_factor
        self.background_mask = masks.create_ellipse_mask(self.background.xsize, self.background.ysize, rel_center_background, inner_radius, angle)

        # Set mask for the source (e.g. segmentation) to None
        self.mask = None

        # Set subtracted box to None
        self.estimated_background = None
        self.estimated_background_cutout = None
        self.subtracted = None
        self.removed = None

        # Set peak position to None initially
        self.peak = None

    # *****************************************************************

    def estimate_background(self, sigma_clip=True, sigma=3.0):

        """
        This function ...
        :return:
        """

        # Perform sigma-clipping on the background if requested
        if sigma_clip: mask = statistics.sigma_clip_mask(self.background, sigma=sigma, mask=self.background_mask)
        else: mask = self.background_mask

        # Estimate the background
        self.estimated_background = self.background.fit_polynomial(3, mask=mask)

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