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
from matplotlib.figure import Figure
from scipy.interpolate import spline
from scipy.signal import argrelextrema
from scipy.signal import find_peaks_cwt

# -----------------------------------------------------------------

class Distribution(object):
    
    """
    This class ...
    """
    
    def __init__(self, counts, edges, centers, mean, median, percentile_16=None, percentile_84=None):
        
        """
        The constructor ...
        :param counts:
        :param edges:
        :param centers:
        :param mean:
        :param median:
        """

        self.counts = counts
        self.edges = edges
        self.centers = centers
        self.mean = mean
        self.median = median
        self.percentile_16 = percentile_16
        self.percentile_84 = percentile_84

    # -----------------------------------------------------------------

    @classmethod
    def from_values(cls, values, bins):

        """
        This function ...
        :param values:
        :param bins:
        :return:
        """

        counts, edges = np.histogram(values, bins=bins, density=True)

        centers = []
        for i in range(len(edges) - 1):
            centers.append(0.5 * (edges[i] + edges[i + 1]))

        mean = np.mean(values)
        median = np.median(values)

        percentile_16 = np.percentile(values, 15.86)
        percentile_84 = np.percentile(values, 84.14)

        return cls(counts, edges, centers, mean, median, percentile_16, percentile_84)

    # -----------------------------------------------------------------

    @property
    def bins(self):

        """
        This function ...
        :return:
        """

        return len(self.centers)

    # -----------------------------------------------------------------

    @property
    def bin_widths(self):

        """
        This function ...
        :return:
        """

        widths = []
        for i in range(len(self.edges) - 1):
            widths.append(self.edges[i+1] - self.edges[i])
        return widths

    # -----------------------------------------------------------------

    @property
    def bin_width(self):

        """
        This function ...
        :return:
        """

        widths = self.bin_widths

        if not all_close(widths):

            print(widths)
            raise RuntimeError("Bin widths not equal")

        width = (self.max_value - self.min_value) / self.bins

        assert np.isclose(width, widths[0])

        return width

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
    def least_frequent(self):

        """
        This function ...
        :return:
        """

        index = np.argmin(self.counts)
        return self.centers[index]

    # -----------------------------------------------------------------

    @property
    def least_frequent_non_zero(self):

        """
        This function ...
        :return:
        """

        ma = np.ma.masked_equal(self.counts, 0.0, copy=False)
        index = np.argmin(ma)
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
    def min_count(self):

        """
        This function ...
        :return:
        """

        return np.min(self.counts)

    # -----------------------------------------------------------------

    @property
    def min_count_nonzero(self):

        """
        This function ...
        :return:
        """

        ma = np.ma.masked_equal(self.counts, 0.0, copy=False)
        return np.min(ma)

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

    def plot(self, title=None, path=None, logscale=False):

        """
        This function ...
        :param title:
        :param path:
        :param logscale:
        :return:
        """

        # Create a canvas to place the subgraphs
        canvas = plt.figure()
        #canvas = Figure()
        rect = canvas.patch
        rect.set_facecolor('white')

        sp1 = canvas.add_subplot(1, 1, 1, axisbg='w')

        sp1.bar(self.edges[:-1], self.counts, linewidth=0, width=self.bin_width, alpha=0.75)

        sp1.set_xlim(0.8 * self.min_value, 1.2 * self.max_value)

        if logscale: sp1.set_ylim(0.5*self.min_count_nonzero, 2.0*self.max_count)
        else: sp1.set_ylim(0, 1.1 * self.max_count)

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
        sp1.set_ylabel('Probability', color='red')

        if logscale: sp1.set_yscale("log", nonposx='clip')

        plt.tight_layout()
        plt.grid(alpha=0.8)

        if path is None: plt.show()
        else: canvas.savefig(path)

# -----------------------------------------------------------------

def all_equal(array):

    """
    This function ...
    :param array:
    :return:
    """

    first = array[0]

    for i in range(1, len(array)):
        if array[i] != first: return False

    return True

# -----------------------------------------------------------------

def all_close(array):

    """
    This function ...
    :param array:
    :return:
    """

    first = array[0]

    for i in range(1, len(array)):
        if not np.isclose(array[i], first): return False

    return True

# -----------------------------------------------------------------
