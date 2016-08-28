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
from ..tools.logging import log
from ...magic.tools.plotting import line_styles, filled_markers, pretty_colours

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

        # The definite axes limits
        self.min_value = None
        self.max_value = None
        self.min_count = None
        self.max_count = None

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

    def run(self, output_path=None, min_value=None, max_value=None, min_count=None, max_count=None, format=None, logscale=False, legend=True):

        """
        This function ...
        :param output_path:
        :param min_value:
        :param max_value:
        :param min_count:
        :param max_count:
        :param format:
        :param logscale:
        :return:
        """

        # Set the axis limits
        self.min_value = min_value
        self.max_value = max_value
        self.min_count = min_count
        self.max_count = max_count

        # Make the plot
        self.plot(output_path, format=format, logscale=logscale, legend=legend)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the distribution plotter ...")

        # Set default values for all attributes
        self.title = None
        self.distributions = OrderedDict()
        self._min_value = None
        self._max_value = None
        self._min_count = None
        self._max_count = None
        self.min_value = None
        self.max_value = None
        self.min_count = None
        self.max_count = None
        self._figure = None

    # -----------------------------------------------------------------

    def plot(self, path=None, logscale=False, add_smooth=False, add_extrema=False, format=None, add_statistics=True, legend=True):

        """
        This function ...
        :param path:
        :param logscale:
        :param add_smooth:
        :param add_extrema:
        :param format:
        :param add_statistics:
        :param legend:
        :return:
        """

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
            plt.bar(distribution.edges[:-1], distribution.counts, linewidth=0, width=distribution.bin_width, alpha=0.5, color=next(colors))

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
            ax.bar(distribution.edges[:-1], distribution.counts, linewidth=0, width=distribution.bin_width, alpha=0.5, color=pretty_colours[0])

            if label in self.extra_distributions:

                extra_distribution = self.extra_distributions[label]

                # Plot the distribution as a histogram
                ax.bar(extra_distribution.edges[:-1], extra_distribution.counts, linewidth=0, width=extra_distribution.bin_width, alpha=0.5, color=pretty_colours[1])

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
