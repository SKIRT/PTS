#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np

# Import Astromagic modules
from ..tools import cropping
from ..tools import fitting
from ..tools import statistics
from .vector import Position

# Import astronomical modules
from astropy.table import Table
from photutils import find_peaks

# *****************************************************************

class Box(np.ndarray):

    """
    This class ...
    """

    def __new__(cls, frame, center, radius, angle):

        """
        This function ...
        :param cls:
        :param input_array:
        :param info:
        :return:
        """

        x_extent = np.abs(np.sin(angle)*radius.x + np.cos(angle)*radius.y)
        y_extent = np.abs(np.cos(angle)*radius.x - np.sin(angle)*radius.y)

        # Crop the frame
        cropped, x_min, x_max, y_min, y_max = cropping.crop(frame, center.x, center.y, x_extent, y_extent)

        # Check that the center position lies within the box
        assert (x_min <= center.x < x_max and y_min <= center.y < y_max)

        # Create an object of the Box class
        obj = np.asarray(cropped).view(cls)

        # Set attributes of the object
        obj.x_min = x_min
        obj.x_max = x_max
        obj.y_min = y_min
        obj.y_max = y_max

        # Return the object
        return obj

    # *****************************************************************

    def __array_finalize__(self, obj):

        # see InfoArray.__array_finalize__ for comments
        if obj is None: return

        self.x_min = getattr(obj, 'x_min', None)
        self.x_max = getattr(obj, 'x_max', None)
        self.y_min = getattr(obj, 'y_min', None)
        self.y_max = getattr(obj, 'y_max', None)

    # *****************************************************************

    def plot(self):
        
        pass

    # *****************************************************************

    @property
    def xsize(self):

        return self.shape[1]

    # *****************************************************************

    @property
    def ysize(self):

        return self.shape[0]

    # *****************************************************************

    def rel_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the relative position
        return Position(x=position.x-self.x_min, y=position.y-self.y_min)

    # *****************************************************************

    def abs_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the absolute position
        return Position(x=position.x+self.x_min, y=position.y+self.y_min)

    # *****************************************************************

    def locate_peaks(self, threshold_sigmas, remove_gradient=False, central_radius=5.0):

        """
        This function ...
        :return:
        """

        # If requested, remove the gradient background before peak detection
        if remove_gradient:

            # Create a central mask
            central_mask = np.zeros_like(self, dtype=bool)
            central_mask[int(round(self.ysize/2.0-central_radius)):int(round(self.ysize/2.0+central_radius)), int(round(self.xsize/2.0-central_radius)):int(round(self.xsize/2.0+central_radius))] = True

            if np.all(central_mask): return []

            # Fit a polynomial to the background
            poly = fitting.fit_polynomial(self, 3, mask=central_mask)
            polynomial_box = fitting.evaluate_model(poly, 0, self.xsize, 0, self.ysize)

            #if debug: plotting.plot_difference_model(box, poly)

            # Subtract the gradient from the box
            self -= polynomial_box

        # Calculate the sigma-clipped statistics of the frame and find the peaks
        mean, median, stddev = statistics.sigma_clipped_statistics(self, sigma=3.0)
        threshold = median + (threshold_sigmas * stddev)
        peaks = find_peaks(self, threshold, box_size=5)

        # For some reason, once in a while, an ordinary list comes out of the find_peaks routine instead of an
        # Astropy Table instance. We assume we need an empty table in this case
        if type(peaks) is list: peaks = Table([[], []], names=('x_peak', 'y_peak'))

        # Initialize a list to contain the peak positions
        positions = []

        # Loop over the peaks
        for peak in peaks:

            # Calculate the absolute x and y coordinate of the peak
            x = peak['x_peak'] + self.x_min
            y = peak['y_peak'] + self.y_min

            # Add the coordinates to the positions list
            positions.append(Position(x=x,y=y))

        # Return the list of peak positions
        return positions

# *****************************************************************
