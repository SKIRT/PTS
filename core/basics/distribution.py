#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.distribution Contains the Distribution class.

# -----------------------------------------------------------------

# Import standard modules
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from scipy.signal import argrelextrema
from scipy.signal import find_peaks_cwt

# -----------------------------------------------------------------

class Distribution(object):
    
    """
    This class ...
    """
    
    def __init__(self, values, bins):
        
        """
        The constructor ...
        """

        self.values = np.array(values)
        self.bins = bins
        self.counts, self.edges = np.histogram(values, bins=bins, density=True)

    # -----------------------------------------------------------------

    @property
    def bin_width(self):

        """
        This function ...
        :return:
        """

        return (self.max_value - self.min_value) / self.bins

    # -----------------------------------------------------------------

    @property
    def mean(self):

        """
        This function ...
        :return:
        """

        return np.mean(self.values)

    # -----------------------------------------------------------------

    @property
    def median(self):

        """
        This function ...
        :return:
        """

        return np.median(self.values)

    # -----------------------------------------------------------------

    @property
    def most_frequent(self):

        """
        This function ...
        :return:
        """

        index = np.argmax(self.counts)
        return self.centers[index]

    # -----------------------------------------------------------------

    @property
    def min_value(self):

        """
        This function ...
        :return:
        """

        return self.edges[0]

    # -----------------------------------------------------------------

    @property
    def max_value(self):

        """
        This function ...
        :return:
        """

        return self.edges[-1]

    # -----------------------------------------------------------------

    @property
    def max_count(self):

        """
        This function ...
        :return:
        """

        return np.max(self.counts)

    # -----------------------------------------------------------------

    @property
    def centers(self):

        """
        This function ...
        :return:
        """

        centers = []

        for i in range(len(self.edges)-1):

            centers.append(0.5*(self.edges[i]+self.edges[i+1]))

        return centers

    # -----------------------------------------------------------------

    @property
    def local_maxima(self):

        """
        This function ...
        :return:
        """

        m = argrelextrema(self.counts, np.greater)

        centers = np.array(self.centers)

        x = [centers[i] for i in m]
        y = [self.counts[i] for i in m]

        return x, y

    # -----------------------------------------------------------------

    @property
    def local_minima(self):

        """
        This function ...
        :return:
        """

        m = argrelextrema(self.counts, np.less)

        centers = np.array(self.centers)

        x = [centers[i] for i in m]
        y = [self.counts[i] for i in m]

        return x, y

    # -----------------------------------------------------------------

    @property
    def smooth(self):

        """
        This function ...
        :return:
        """

        x_smooth = np.linspace(self.min_value, self.max_value, 200)
        y_smooth = spline(self.centers, self.counts, x_smooth)

        return x_smooth, y_smooth

    # -----------------------------------------------------------------

    @property
    def local_maxima_smooth(self):

        """
        This function ...
        :return:
        """

        x_smooth, y_smooth = self.smooth

        m = argrelextrema(y_smooth, np.greater)

        x = [x_smooth[i] for i in m]
        y = [y_smooth[i] for i in m]

        return x, y

    # -----------------------------------------------------------------

    @property
    def local_minima_smooth(self):

        """
        This function ...
        :return:
        """

        x_smooth, y_smooth = self.smooth

        m = argrelextrema(y_smooth, np.less)

        x = [x_smooth[i] for i in m]
        y = [y_smooth[i] for i in m]

        return x, y

    # -----------------------------------------------------------------

    def plot(self, title=None):

        """
        This function ...
        :return:
        """

        # Create a canvas to place the subgraphs
        canvas = plt.figure()
        rect = canvas.patch
        rect.set_facecolor('white')

        sp1 = canvas.add_subplot(1, 1, 1, axisbg='w')

        sp1.bar(self.edges[:-1], self.counts, linewidth=0, width=self.bin_width, alpha=0.75)

        sp1.set_xlim(0.8 * self.min_value, 1.2 * self.max_value)
        sp1.set_ylim(0, 1.1 * self.max_count)

        x_smooth, y_smooth = self.smooth
        sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)

        x, y = self.local_maxima
        sp1.plot(x, y, 'rs')

        sp1.axvline(self.mean, color="green", linestyle="dashed")
        sp1.axvline(self.median, color="purple", linestyle="dashed")
        sp1.axvline(self.most_frequent, color="orange", linestyle="dashed")

        # Colorcode the tick tabs
        sp1.tick_params(axis='x', colors='red')
        sp1.tick_params(axis='y', colors='red')

        # Colorcode the spine of the graph
        sp1.spines['bottom'].set_color('r')
        sp1.spines['top'].set_color('r')
        sp1.spines['left'].set_color('r')
        sp1.spines['right'].set_color('r')

        # Put the title and labels
        if title is not None: sp1.set_title(title, color='red')
        sp1.set_xlabel('Values', color='red')
        sp1.set_ylabel('Frequency', color='red')

        plt.tight_layout()
        plt.grid(alpha=0.8)

        plt.show()

# -----------------------------------------------------------------
