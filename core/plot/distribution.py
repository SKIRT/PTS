#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.distribution Contains the DistributionPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import matplotlib.pyplot as plt
from collections import OrderedDict
from textwrap import wrap
#import seaborn as sns
import matplotlib.gridspec as gridspec

# Import the relevant PTS classes and modules
from ..basics.log import log
from ...magic.tools.plotting import line_styles, filled_markers, pretty_colours
from ..tools import strings

# -----------------------------------------------------------------

def plot_distribution(distribution, path=None, logscale=False, logfrequency=False, title=None, x_limits=None,
                      y_limits=None, color="xkcd:sky blue", x_label="Values", y_label="Frequency"):

    """
    This function ...
    :param distribution:
    :param path:
    :param logscale:
    :param logfrequency:
    :param title:
    :param x_limits:
    :param y_limits:
    :param color:
    :param x_label:
    :param y_label:
    :return:
    """

    #def plot(self, title=None, path=None, logscale=False, xlogscale=False, x_limits=None, y_limits=None,
    #         add_smooth=False, format=None, add_extrema=False, model=None):

    # Create a canvas to place the subgraphs
    figure = plt.figure()
    rect = figure.patch
    rect.set_facecolor('white')

    #sp1 = canvas.add_subplot(1, 1, 1, axisbg='w')
    #sp1 = canvas.add_subplot(111)
    sp1 = figure.gca()

    #sp1.bar(self.edges[:-1], self.counts, linewidth=0, width=self.bin_width, alpha=0.5)

    #print(distribution.values)
    #print(distribution.edges_log)
    #print(distribution.edges)
    #print(distribution.frequencies)

    alpha = None
    edgecolor = "black"
    linewidth = 1

    edgecolor_list = [edgecolor for _ in range(distribution.nvalues)]
    linewidth_list = [linewidth for _ in range(distribution.nvalues)]

    if logscale: sp1.bar(distribution.edges_log[:-1], distribution.frequencies, width=distribution.bin_widths_log, linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)
    else: sp1.bar(distribution.edges[:-1], distribution.frequencies, width=distribution.bin_widths, linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)

    #print("min", distribution.min_value)
    #print("max", distribution.max_value)

    # Determine the x limits
    if x_limits is None:
        #x_min = distribution.min_value
        #x_max = distribution.max_value
        if logscale:
            x_min = distribution.min_edge_log
            x_max = distribution.max_edge_log
        else:
            x_min = distribution.min_edge
            x_max = distribution.max_edge
    else:
        x_min = x_limits[0]
        x_max = x_limits[1]

    # Determine the y limits
    if y_limits is None:
        y_min = 0. if not logfrequency else 0.5 * distribution.min_frequency_nonzero
        y_max = 1.1 * distribution.max_frequency if not logfrequency else 2. * distribution.max_frequency
    else:
        y_min = y_limits[0]
        y_max = y_limits[1]

    # Set the axis limits
    sp1.set_xlim(x_min, x_max)
    sp1.set_ylim(y_min, y_max)

    # Add smooth
    #if add_smooth:
    #    if logscale:
    #        x_smooth, y_smooth = self.smooth_values_log(x_min=x_min, x_max=x_max)
    #        sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)
    #    else:
    #        x_smooth, y_smooth = self.smooth_values(x_min=x_min, x_max=x_max)
    #        sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)

    #if add_extrema:
    #    x, y = self.local_maxima
    #    sp1.plot(x, y, 'g^')
    #    x, y = self.local_minima
    #    sp1.plot(x, y, 'rv')

    #if model is not None: sp1.plot(self.centers, model(self.centers), label='Model')

    #if logscale: print("mean", distribution.geometric_mean)
    #else: print("mean", distribution.mean)
    #print("median", distribution.median)
    #print("max", distribution.most_frequent)

    #print(self.mean, self.median, self.most_frequent)
    if logscale: mean_line = sp1.axvline(distribution.geometric_mean, color="green", linestyle="dashed", label="Mean")
    else: mean_line = sp1.axvline(distribution.mean, color="green", linestyle="dashed", label="Mean")
    median_line = sp1.axvline(distribution.median, color="purple", linestyle="dashed", label="Median")
    max_line = sp1.axvline(distribution.most_frequent, color="orange", linestyle="dashed", label="Most frequent")
    plt.legend()

    # Colorcode the tick tabs
    #sp1.tick_params(axis='x', colors='red')
    #sp1.tick_params(axis='y', colors='red')

    # Colorcode the spine of the graph
    #sp1.spines['bottom'].set_color('r')
    #sp1.spines['top'].set_color('r')
    #sp1.spines['left'].set_color('r')
    #sp1.spines['right'].set_color('r')

    if distribution.unit is not None: xlabel = x_label + " [" + str(distribution.unit) + "]"
    else: xlabel = x_label

    # Put the title and labels
    if title is not None: sp1.set_title(strings.split_in_lines(title))
    sp1.set_xlabel(xlabel)
    sp1.set_ylabel(y_label)

    if logfrequency: sp1.set_yscale("log", nonposx='clip')
    if logscale: sp1.set_xscale("log")

    plt.tight_layout()
    #plt.grid(alpha=0.8)

    if path is None: plt.show()
    else: figure.savefig(path, format=format)

    # Close the figure
    plt.close()

# -----------------------------------------------------------------

class DistributionPlotter(object):
    
    """
    This class ...
    """

    def __init__(self, title=None):

        """
        This function ...
        :return:
        """

        # Set the title
        self.title = title

        # The name of the variable
        self.name = None

        # The different distributions
        self.distributions = OrderedDict()

        # Keep track of the minimal and maximal value and count encountered during the plotting
        self._min_value = None
        self._max_value = None
        self._min_count = None
        self._max_count = None

        # Input
        self.min_value = None
        self.max_value = None
        self.min_count = None
        self.max_count = None
        self.output_path = None
        self.format = None
        self.logscale = False
        self.legend = True
        self.add_smooth = False
        self.add_extrema = False
        self.add_statistics = False

        # Store the figure and its axes as references
        self._figure = None

        # Properties
        self.transparent = False

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

    # -----------------------------------------------------------------

    def set_variable_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        self.name = name

    # -----------------------------------------------------------------

    def add_distribution(self, distribution, label):

        """
        This function ...
        :param distribution:
        :param label:
        :return:
        """

        self.distributions[label] = distribution

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Make the plot
        self.plot()

    # -----------------------------------------------------------------

    def clear(self, clear_distributions=True):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the distribution plotter ...")

        # Set default values for all attributes
        self.title = None
        if clear_distributions: self.distributions = OrderedDict()

        self._min_value = None
        self._max_value = None
        self._min_count = None
        self._max_count = None

        # Input
        self.min_value = None
        self.max_value = None
        self.min_count = None
        self.max_count = None
        self.output_path = None
        self.format = None
        self.logscale = False
        self.legend = True
        self.add_smooth = False
        self.add_extrema = False
        self.add_statistics = False

        self._figure = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Set the axis limits
        self.min_value = kwargs.pop("min_value", None)
        self.max_value = kwargs.pop("max_value", None)
        self.min_count = kwargs.pop("min_count", None)
        self.max_count = kwargs.pop("max_count", None)

        # Get ...
        self.output_path = kwargs.pop("output_path", None)
        self.format = kwargs.pop("format", None)
        self.logscale = kwargs.pop("logscale", False)
        self.legend = kwargs.pop("legend", False)

        self.add_smooth = kwargs.pop("add_smooth", False)
        self.add_extrema = kwargs.pop("add_extrema", False)
        self.add_statistics = kwargs.pop("add_statistics", False)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        path = self.output_path
        logscale = self.logscale
        format = self.format
        legend = self.legend
        add_smooth = self.add_smooth
        add_extrema = self.add_extrema
        add_statistics = self.add_statistics

        # Inform the user
        log.info("Making the distribution plot ...")

        # Plot the attenuation curves
        plt.figure(figsize=(10, 10))

        # Create the figure
        self._figure = plt.figure()

        colors = iter(pretty_colours)

        # Plot the distributions
        for label in self.distributions:

            distribution = self.distributions[label]

            # Plot the distribution as a histogram
            plt.bar(distribution.edges[:-1], distribution.counts, linewidth=0, width=distribution.bin_width, alpha=0.5, color=next(colors), align="edge")

            min_value = distribution.min_value
            max_value = distribution.max_value
            min_count = 0. if not logscale else 0.5 * distribution.min_count_nonzero
            max_count = 1.1 * distribution.max_count if not logscale else 2. * distribution.max_count

            # Keep track of minimum and maximum value
            if self._min_value is None or min_value < self._min_value: self._min_value = min_value
            if self._max_value is None or max_value > self._max_value: self._max_value = max_value

            # Keep track of the minimal and maximal count
            if self._min_count is None or min_count < self._min_count: self._min_count = min_count
            if self._max_count is None or max_count > self._max_count: self._max_count = max_count

            # Add smooth
            if add_smooth:

                if logscale:
                    x_smooth, y_smooth = distribution.smooth_values_log(x_min=self._min_value, x_max=self._max_value)
                    plt.plot(x_smooth, y_smooth, 'red', linewidth=1)
                else:
                    x_smooth, y_smooth = distribution.smooth_values(x_min=self._min_value, x_max=self._max_value)
                    plt.plot(x_smooth, y_smooth, 'red', linewidth=1)

            if add_extrema:

                x, y = distribution.local_maxima
                plt.plot(x, y, 'g^')

                x, y = distribution.local_minima
                plt.plot(x, y, 'rv')

            if add_statistics:

                plt.axvline(distribution.mean, color="green", linestyle="dashed")
                plt.axvline(distribution.median, color="purple", linestyle="dashed")
                plt.axvline(distribution.most_frequent, color="orange", linestyle="dashed")

        # Axis limits are now definite
        if self.min_value is None: self.min_value = self._min_value
        if self.max_value is None: self.max_value = self._max_value
        if self.min_count is None: self.min_count = self._min_count
        if self.max_count is None: self.max_count = self._max_count

        axes = plt.gca()

        # Set the axis limits
        axes.set_xlim(self.min_value, self.max_value)
        axes.set_ylim(self.min_count, self.max_count)

        # Set the axis labels
        if self.name is not None: axes.set_xlabel(self.name)
        axes.set_ylabel("Normalized count")

        # Set the title
        if self.title is not None: axes.set_title("\n".join(wrap(self.title, 60)))

        if logscale: axes.set_yscale("log", nonposx='clip')

        # Add legend
        if len(self.distributions) > 1 and legend: plt.legend()

        plt.tight_layout()

        # Debugging
        if path is None: log.debug("Showing the plot interactively ...")
        elif type(path).__name__ == "BytesIO": log.debug("Saving the distribution plot to a buffer ...")
        else: log.debug("Saving the distribution plot to " + str(path) + " ...")

        # Finish
        if path is None: plt.show()
        else: plt.savefig(path, bbox_inches='tight', pad_inches=0.25, format=format, transparent=self.transparent) # Save the figure
        plt.close()

# -----------------------------------------------------------------

class DistributionGridPlotter(object):

    """
    This class ...
    """

    def __init__(self, title=None):

        """
        The constructor ...
        """

        # Set the title
        self.title = title

        # The different distributions
        self.distributions = OrderedDict()
        self.extra_distributions = dict()

        # Properties
        self.format = None
        self.transparent = False
        self.ncols = 7
        self.width = 16

    # -----------------------------------------------------------------

    def add_distribution(self, distribution, label):

        """
        This function ...
        :return:
        """

        if label in self.distributions: self.extra_distributions[label] = distribution
        else: self.distributions[label] = distribution

    # -----------------------------------------------------------------

    def run(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Make the plot
        self.plot(path)

    # -----------------------------------------------------------------

    @property
    def npanels(self):

        """
        This function ...
        :return:
        """

        return len(self.distributions)

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :parm path:
        :return:
        """

        # Determine the necessary number of rows
        nrows = int(math.ceil(self.npanels / self.ncols))

        ratio = float(nrows) / float(self.ncols)
        height = ratio * self.width

        # Create the figure
        self._figure = plt.figure(figsize=(self.width, height))

        self._figure.subplots_adjust(hspace=0.0, wspace=0.0)


        gs = gridspec.GridSpec(nrows, self.ncols, wspace=0.0, hspace=0.0)

        # Loop over the distributions
        counter = 0
        ax = None
        for label in self.distributions:

            row = int(counter / self.ncols)
            col = counter % self.ncols

            subplotspec = gs[row, col]

            ax = plt.subplot(subplotspec)

            distribution = self.distributions[label]

            # Plot the distribution as a histogram
            ax.bar(distribution.edges[:-1], distribution.counts, linewidth=0, width=distribution.bin_width, alpha=0.5, color=pretty_colours[0], align="edge")

            if label in self.extra_distributions:

                extra_distribution = self.extra_distributions[label]

                # Plot the distribution as a histogram
                ax.bar(extra_distribution.edges[:-1], extra_distribution.counts, linewidth=0, width=extra_distribution.bin_width, alpha=0.5, color=pretty_colours[1], align="edge")

            counter += 1

        # Finish
        self.finish_plot(path)

    # -----------------------------------------------------------------

    def finish_plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Set the title
        if self.title is not None: self._figure.suptitle("\n".join(wrap(self.title, 60)))

        # plt.tight_layout()

        # Debugging
        if type(path).__name__ == "BytesIO": log.debug("Saving the distribution plot to a buffer ...")
        elif path is None: log.debug("Showing the distribution plot ...")
        else: log.debug("Saving the distribution plot to " + str(path) + " ...")

        if path is not None:
            # Save the figure
            plt.savefig(path, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent, format=self.format)
        else: plt.show()
        plt.close()

# -----------------------------------------------------------------
