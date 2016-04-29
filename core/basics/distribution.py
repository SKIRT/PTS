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
from scipy.stats import rv_continuous
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy.interpolate import spline
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from scipy.signal import find_peaks_cwt
from scipy.interpolate import InterpolatedUnivariateSpline

# Import astronomical modules
from astropy.table import Table

# -----------------------------------------------------------------

class Distribution(object):
    
    """
    This class ...
    """
    
    def __init__(self, counts, edges, centers, mean, median, percentile_16, percentile_84, name=None):
        
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
        self.mean = float(mean)
        self.median = float(median)
        self.percentile_16 = float(percentile_16) if percentile_16 is not None else None
        self.percentile_84 = float(percentile_84) if percentile_84 is not None else None

        self.name = name

        self._rv = None

    # -----------------------------------------------------------------

    @classmethod
    def from_probabilities(cls, probabilities, values, name=None):

        """
        This function ...
        :param probabilities:
        :param values:
        :param name:
        :return:
        """

        # Calculate the mean value (weigh each value by its probability and calculate the normalized sum)
        mean = np.sum(values * probabilities) / np.sum(probabilities)

        # Calculate the edges
        edges = [0]
        for i in range(len(values) - 1):
            edges.append(0.5 * (values[i] + values[i + 1]))
        edges[0] = values[0] - (edges[1] - values[0])
        edges.append(values[-1] + (values[-1] - edges[-1]))

        # Calculate the 16, 50 and 84 percentiles
        percentile_16, median, percentile_84 = find_percentiles(values, probabilities)

        # Return a new Distribution instance
        centers = values
        return cls(probabilities, edges, centers, mean, median, percentile_16, percentile_84, name)

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
        for i in range(len(edges) - 1): centers.append(0.5 * (edges[i] + edges[i + 1]))

        mean = np.mean(values)
        median = np.median(values)

        percentile_16 = np.percentile(values, 15.86)
        percentile_84 = np.percentile(values, 84.14)

        return cls(counts, edges, centers, mean, median, percentile_16, percentile_84)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Read the table from file
        fill_values = [('--', '0')]
        table = Table.read(path, fill_values=fill_values, format="ascii.ecsv")

        parameter_name = table.colnames[0]

        values = table[parameter_name]
        probabilities = np.array(table["Probability"])

        mean = np.sum(values * probabilities) / np.sum(probabilities)

        edges = [0]

        for i in range(len(values) - 1):
            edges.append(0.5 * (values[i] + values[i + 1]))

        edges[0] = values[0] - (edges[1] - values[0])
        edges.append(values[-1] + (values[-1] - edges[-1]))

        centers = values

        median = table.meta["median"]
        percentile_16 = table.meta["percentile16"]
        percentile_84 = table.meta["percentile84"]

        return cls(probabilities, edges, centers, mean, median, percentile_16, percentile_84)

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        data = [self.centers, self.counts]
        names = [self.name, "Probability"]
        meta = {"mean": self.mean, "median": self.median, "percentile16": self.percentile_16, "percentile84": self.percentile_84}

        table = Table(data, names=names, meta=meta, masked=True)

        # Save the table
        table.write(path, format="ascii.ecsv")

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

        m = argrelextrema(self.counts, np.greater)[0].tolist()

        centers = np.array(self.centers)

        # Find the index of the absolute maximum (should also be included, is not for example when it is at the edge)
        index = np.argmax(self.counts)
        if index not in m: m.append(index)

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

        m = argrelextrema(self.counts, np.less)[0].tolist()

        centers = np.array(self.centers)

        # Find the indx of the absolute minimum (should also be included, is not for example when it is at the edge)
        index = np.argmin(self.counts)
        if index not in m: m.append(index)

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

        order = 2
        s = InterpolatedUnivariateSpline(self.centers, self.counts, k=order)

        # Return the spline curve
        return s

    # -----------------------------------------------------------------

    @property
    def smooth_log(self):

        """
        This function ...
        :return:
        """

        centers = np.array(self.centers)
        counts = np.array(self.counts)

        not_zero = counts != 0

        centers = centers[not_zero]
        counts = counts[not_zero]

        order = 2
        s = InterpolatedUnivariateSpline(centers, np.log10(counts), k=order)

        # Return the spline curve
        return s

    # -----------------------------------------------------------------

    @property
    def smooth_values(self):

        """
        This function ...
        :return:
        """

        x_smooth = np.linspace(self.min_value, self.max_value, 200)

        s = self.smooth
        y_smooth = s(x_smooth)

        return x_smooth, y_smooth

    # -----------------------------------------------------------------

    @property
    def smooth_values_log(self):

        """
        This function ...
        :return:
        """

        x_smooth = np.linspace(self.min_value, self.max_value, 200)

        s = self.smooth_log
        y_smooth_log = s(x_smooth)

        return x_smooth, 10.**y_smooth_log

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

    def get_rv(self):

        """
        This function ...
        :return:
        """

        if self._rv is None:

            class RVDistribution(rv_continuous):
                def _pdf(self, x): return np.exp(-x ** 2 / 2.) / np.sqrt(2.0 * np.pi) # TODO

            self._rv = RVDistribution()

        return self._rv

    # -----------------------------------------------------------------

    def random(self):

        """
        This function ...
        :return:
        """

        # Get the distribution function (in Scipy's API)
        rv_function = self.get_rv()

        # Get a random variate
        return rv_function.rvs()

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

        if logscale:

            x_smooth, y_smooth = self.smooth_values_log
            sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)

        else:

            x_smooth, y_smooth = self.smooth_values
            sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)

        x, y = self.local_maxima
        sp1.plot(x, y, 'g^')

        x, y = self.local_minima
        sp1.plot(x, y, 'rv')

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

def find_percentiles(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    if len(values) > 1: return find_percentile_16(values, probabilities), find_percentile_50(values, probabilities), find_percentile_84(values, probabilities)
    else: return None, None, None

# -----------------------------------------------------------------

def find_percentile_16(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    find_percentile(values, probabilities, 15.86)

# -----------------------------------------------------------------

def find_percentile_50(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    return find_percentile(values, probabilities, 50.)

# -----------------------------------------------------------------

def find_percentile_84(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities:
    :return:
    """

    return find_percentile(values, probabilities, 84.14)

# -----------------------------------------------------------------

def find_percentile(values, probabilities, percentile):

    """
    This function ...
    :param values:
    :param probabilities:
    :param percentile:
    :return:
    """

    npoints = 10000
    interpfunc = interp1d(values, probabilities, kind='linear')

    parRange = np.linspace(min(values), max(values), npoints)
    interProb = interpfunc(parRange)

    cumInteg = np.zeros(npoints-1)

    for i in range(1,npoints-1):
        cumInteg[i] = cumInteg[i-1] + (0.5*(interProb[i+1] + interProb[i]) * (parRange[i+1] - parRange[i]))

    cumInteg = cumInteg / cumInteg[-1]
    idx = (np.abs(cumInteg-percentile/100.)).argmin()

    return parRange[idx]

# -----------------------------------------------------------------
