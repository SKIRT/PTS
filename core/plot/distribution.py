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
from ...magic.tools.plotting import pretty_colours
from ..tools import strings
from ..basics.configurable import Configurable
from ..basics.containers import DefaultOrderedDict
from ..tools import filesystem as fs
from ..basics.distribution import Distribution
from ..tools import types
from ..basics.plot import MPLFigure, BokehFigure, BokehPlot, mpl, bokeh
from ..tools.stringify import tostr

# -----------------------------------------------------------------

def plot_distribution(distribution, path=None, logscale=False, logfrequency=False, title=None, x_limits=None,
                      y_limits=None, color="xkcd:sky blue", x_label=None, y_label=None):

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

    if x_label is None: x_label = distribution.value_name
    if y_label is None: y_label = distribution.y_name
    if distribution.unit is not None: x_label += " [" + str(distribution.unit) + "]"

    # Put the title and labels
    if title is not None: sp1.set_title(strings.split_in_lines(title))
    sp1.set_xlabel(x_label)
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

def plot_cumulative(distribution, title=None, path=None, logscale=False, x_limits=None, y_limits=None, npoints=200):

    """
    This function ...
    :param distribution:
    :param title:
    :param path:
    :param logscale:
    :param x_limits:
    :param y_limits:
    :param npoints:
    :return:
    """

    pass

# -----------------------------------------------------------------

def plot_cumulative_smooth(distribution, title=None, path=None, logscale=False, x_limits=None, y_limits=None, npoints=200):

    """
    This function ...
    :param distribution:
    :param title:
    :param path:
    :param logscale:
    :param x_limits:
    :param y_limits:
    :param npoints:
    :return:
    """

    # Create a canvas to place the subgraphs
    figure = plt.figure()
    rect = figure.patch
    rect.set_facecolor('white')

    #sp1 = canvas.add_subplot(1, 1, 1, axisbg='w')
    #sp1 = figure.add_subplot(111)
    sp1 = figure.gca()

    # Determine the x limits
    if x_limits is None:
        x_min = 0.8 * self.min_value
        x_max = 1.2 * self.max_value
    else:
        x_min = x_limits[0]
        x_max = x_limits[1]

    # Determine the y limits
    if y_limits is None:
        y_min = 0. if not logscale else 0.5 * self.min_count_nonzero
        y_max = 1.1 * self.max_count if not logscale else 2. * self.max_count
    else:
        y_min = y_limits[0]
        y_max = y_limits[1]

    # Set the axis limits
    sp1.set_xlim(x_min, x_max)
    sp1.set_ylim(y_min, y_max)

    if logscale:
        x_smooth, y_smooth = self.cumulative_log_smooth(x_min=x_min, x_max=x_max, npoints=npoints)
        sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)
    else:
        x_smooth, y_smooth = self.cumulative_smooth(x_min=x_min, x_max=x_max, npoints=npoints)
        sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)

    #print(self.mean, self.median, self.most_frequent)
    mean_line = sp1.axvline(self.mean, color="green", linestyle="dashed", label="Mean")
    median_line = sp1.axvline(self.median, color="purple", linestyle="dashed", label="Median")
    max_line = sp1.axvline(self.most_frequent, color="orange", linestyle="dashed", label="Max")
    plt.legend()

    # Colorcode the tick tabs
    sp1.tick_params(axis='x', colors='red')
    sp1.tick_params(axis='y', colors='red')

    # Colorcode the spine of the graph
    sp1.spines['bottom'].set_color('r')
    sp1.spines['top'].set_color('r')
    sp1.spines['left'].set_color('r')
    sp1.spines['right'].set_color('r')

    if self.unit is not None: xlabel = "Values [" + str(self.unit) + "]"
    else: xlabel = "Values"

    # Put the title and labels
    if title is not None: sp1.set_title(title, color='red')
    sp1.set_xlabel(xlabel, color='red')
    sp1.set_ylabel('Cumulative probability', color='red')

    if logscale: sp1.set_yscale("log", nonposx='clip')

    plt.tight_layout()
    plt.grid(alpha=0.8)

    if path is None: plt.show()
    else: figure.savefig(path)

# -----------------------------------------------------------------

def plot_smooth(distribution, title=None, path=None, logscale=False, xlogscale=False, x_limits=None, y_limits=None, npoints=200):

    """
    This function ...
    :param distribution:
    :param title:
    :param path:
    :param logscale:
    :param xlogscale:
    :param x_limits:
    :param y_limits:
    :param npoints:
    :return:
    """

    # Create a canvas to place the subgraphs
    figure = plt.figure()
    rect = figure.patch
    rect.set_facecolor('white')

    #sp1 = canvas.add_subplot(1, 1, 1, axisbg='w')
    #sp1 = canvas.add_subplot(111)
    sp1 = figure.gca()

    # Determine the x limits
    if x_limits is None:
        x_min = 0.8 * self.min_value
        x_max = 1.2 * self.max_value
    else:
        x_min = x_limits[0]
        x_max = x_limits[1]

    # Determine the y limits
    if y_limits is None:
        y_min = 0. if not logscale else 0.5 * self.min_count_nonzero
        y_max = 1.1 * self.max_count if not logscale else 2. * self.max_count
    else:
        y_min = y_limits[0]
        y_max = y_limits[1]

    # Set the axis limits
    sp1.set_xlim(x_min, x_max)
    sp1.set_ylim(y_min, y_max)

    if logscale:
        x_smooth, y_smooth = self.smooth_values_log(x_min=x_min, x_max=x_max, npoints=npoints)
        sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)
    else:
        x_smooth, y_smooth = self.smooth_values(x_min=x_min, x_max=x_max, npoints=npoints)
        sp1.plot(x_smooth, y_smooth, 'red', linewidth=1)

    x, y = get_local_maxima(x_smooth, y_smooth)
    sp1.plot(x, y, 'g^')

    x, y = get_local_minima(x_smooth, y_smooth)
    sp1.plot(x, y, 'rv')

    # print(self.mean, self.median, self.most_frequent)
    mean_line = sp1.axvline(self.mean, color="green", linestyle="dashed", label="Mean")
    median_line = sp1.axvline(self.median, color="purple", linestyle="dashed", label="Median")
    max_line = sp1.axvline(self.most_frequent, color="orange", linestyle="dashed", label="Max")
    plt.legend()

    # Colorcode the tick tabs
    sp1.tick_params(axis='x', colors='red')
    sp1.tick_params(axis='y', colors='red')

    # Colorcode the spine of the graph
    sp1.spines['bottom'].set_color('r')
    sp1.spines['top'].set_color('r')
    sp1.spines['left'].set_color('r')
    sp1.spines['right'].set_color('r')

    if self.unit is not None: xlabel = "Values [" + str(self.unit) + "]"
    else: xlabel = "Values"

    # Put the title and labels
    if title is not None: sp1.set_title(title, color='red')
    sp1.set_xlabel(xlabel, color='red')
    sp1.set_ylabel('Probability', color='red')

    if logscale: sp1.set_yscale("log", nonposx='clip')
    if xlogscale: sp1.set_xscale("log")

    plt.tight_layout()
    plt.grid(alpha=0.8)

    if path is None: plt.show()
    else: figure.savefig(path)

# -----------------------------------------------------------------

class DistributionPlotter(Configurable):
    
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(DistributionPlotter, self).__init__(*args, **kwargs)

        # Set the title
        self.title = None

        # The different distributions
        self.distributions = DefaultOrderedDict(OrderedDict)

        # Keep track of the minimal and maximal value and count encountered during the plotting
        self._min_value = None
        self._max_value = None
        self._min_frequency = None
        self._max_frequency = None

        # Min and max values per panel
        self.min_values = dict()
        self.max_values = dict()

        # Input
        self.min_value = None
        self.max_value = None
        self.min_frequency = None
        self.max_frequency = None
        self.format = None
        self.legend = True

        # Add features
        self.add_smooth = False
        self.add_extrema = False
        self.add_statistics = False

        # Logscale
        self.logscale = False
        self.logfrequency = False

        # The figure
        self.figure = None

        # The output path
        self.out_path = None

        # Properties
        self.transparent = False

        # Plots?
        self.main_plot = None
        self.panel_plots = OrderedDict()

    # -----------------------------------------------------------------

    def add_distribution(self, distribution, label, panel="standard"):

        """
        This function ...
        :param distribution:
        :param label:
        :param panel
        :return:
        """

        self.distributions[panel][label] = distribution

    # -----------------------------------------------------------------

    @property
    def npanels(self):

        """
        This function ...
        :return:
        """

        return len(self.distributions)

    # -----------------------------------------------------------------

    @property
    def panels(self):

        """
        This function ...
        :return:
        """

        return self.distributions.keys()

    # -----------------------------------------------------------------

    @property
    def has_single_panel(self):

        """
        This function ...
        :return:
        """

        return self.npanels == 1

    # -----------------------------------------------------------------

    @property
    def single_panel(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_panel: raise ValueError("Not a single panel")
        return self.panels[0]

    # -----------------------------------------------------------------

    def ndistributions_for_panel(self, panel):

        """
        This function ...
        :param panel:
        :return:
        """

        return len(self.distributions[panel])

    # -----------------------------------------------------------------

    @property
    def ndistributions(self):

        """
        This function ...
        :return:
        """

        total = 0
        for panel in self.panels: total += self.ndistributions_for_panel(panel)
        return total

    # -----------------------------------------------------------------

    @property
    def has_single_distribution(self):

        """
        This function ...
        :return:
        """

        return self.ndistributions == 1

    # -----------------------------------------------------------------

    @property
    def single_distribution_label(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_distribution: raise ValueError("Not a single distribution")
        return self.distributions[self.single_panel].keys()[0]

    # -----------------------------------------------------------------

    @property
    def single_distribution(self):

        """
        This function ...
        :return:
        """

        return self.distributions[self.single_panel][self.single_distribution_label]

    # -----------------------------------------------------------------

    @property
    def no_distributions(self):

        """
        This function ...
        :return:
        """

        return self.ndistributions == 0

    # -----------------------------------------------------------------

    def get_distributions_panel(self, panel):

        """
        Thisf unction ...
        :param panel:
        :return:
        """

        return self.distributions[panel].values()

    # -----------------------------------------------------------------

    def get_ndistributions_panel(self, panel):

        """
        This function ...
        :param panel:
        :return:
        """

        return len(self.distributions[panel])

    # -----------------------------------------------------------------

    def has_single_distribution_panel(self, panel):

        """
        Thisf unction ...
        :param panel:
        :return:
        """

        return self.get_ndistributions_panel(panel) == 1

    # -----------------------------------------------------------------

    def get_single_distribution_label_panel(self, panel):

        """
        This function ...
        :param panel:
        :return:
        """

        if not self.has_single_distribution_panel(panel): raise ValueError("Not a single distribution")
        return self.distributions[panel].keys()[0]

    # -----------------------------------------------------------------

    def get_single_distribution_panel(self, panel):

        """
        This function ...
        :param panel:
        :return:
        """

        return self.distributions[panel][self.get_single_distribution_label_panel(panel)]

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Make the plot
        self.plot()

    # -----------------------------------------------------------------

    def clear(self, clear_distributions=True):

        """
        This function ...
        :param clear_distributions:
        :return:
        """

        # Inform the user
        log.info("Clearing the distribution plotter ...")

        # Set default values for all attributes
        self.title = None
        if clear_distributions: self.distributions = DefaultOrderedDict(OrderedDict)

        self._min_value = None
        self._max_value = None
        self._min_frequency = None
        self._max_frequency = None

        # Input
        self.min_value = None
        self.max_value = None
        self.min_frequency = None
        self.max_frequency = None
        self.format = None

        self.logscale = False
        self.logfrequency = False

        self.legend = True
        self.add_smooth = False
        self.add_extrema = False
        self.add_statistics = False

        self.figure = None
        self.main_plot = None
        self.panel_plots = OrderedDict()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DistributionPlotter, self).setup(**kwargs)

        # Set the title
        self.title = kwargs.pop("title", None)

        # Set the axis limits
        self.min_value = kwargs.pop("min_value", None)
        self.max_value = kwargs.pop("max_value", None)
        self.min_frequency = kwargs.pop("min_frequency", None)
        self.max_frequency = kwargs.pop("max_frequency", None)

        # Get ...
        self.format = kwargs.pop("format", None)
        self.legend = kwargs.pop("legend", False)

        # Logscale
        self.logscale = kwargs.pop("logscale", self.config.logscale)
        self.logfrequency = kwargs.pop("logfrequency", self.config.logfrequency)

        # Features
        self.add_smooth = kwargs.pop("add_smooth", self.config.smooth)
        self.add_extrema = kwargs.pop("add_extrema", self.config.extrema)
        self.add_statistics = kwargs.pop("add_statistics", self.config.statistics)

        # Set the output path
        self.out_path = kwargs.pop("output", None)
        if self.out_path is None and "output" in self.config and self.config.output is not None:
            full_output_path = fs.absolute_or_in(self.config.output, self.config.path)
            if fs.has_extension(full_output_path):
                directory_path = fs.directory_of(full_output_path)
                if not fs.is_directory(directory_path): fs.create_directory(directory_path)
            elif not fs.is_directory(full_output_path):
                fs.create_directory(full_output_path)
            self.out_path = full_output_path

        # Add distribution files present in the current working directory (if nothing is added manually)
        if self.no_distributions: self.load_distributions()

        # Initialize min and max values for panels
        for panel in self.panels: self.min_values[panel] = None
        for panel in self.panels: self.max_values[panel] = None

        # Initialize the plot
        if self.has_single_panel: figsize = self.config.plot.figsize
        else:
            xsize = self.config.plot.figsize[0] * self.npanels
            ysize = self.config.plot.figsize[1]
            figsize = (xsize, ysize)
            #print("figsize:", figsize)

        # Create the plot
        if self.config.library == mpl: self.figure = MPLFigure(size=figsize)
        elif self.config.library == bokeh: self.figure = BokehFigure()
        else: raise ValueError("Invalid libary: " + self.config.library)

        # Set the 'show' flag
        if self.config.show is None:
            if self.out_path is not None: self.config.show = False
            else: self.config.show = True

    # -----------------------------------------------------------------

    def load_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading distributions from file ...")

        # From file paths
        if self.config.distributions is not None: self.load_distribution_files()

        # From files in working directory
        else: self.load_distributions_cwd()

    # -----------------------------------------------------------------

    def load_distribution_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading distribution files ...")

        # Loop over the distribution files
        for path in self.config.distributions:

            # Get the SED name
            name = fs.strip_extension(fs.name(path))

            # Debugging
            log.debug("Loading '" + name + "' distribution ...")

            # Load the distribution
            distribution = Distribution.from_file(path)

            # Add the distribution
            if self.config.panels: self.add_distribution(distribution, label=name, panel=name)
            else: self.add_distribution(distribution, label=name)

    # -----------------------------------------------------------------

    def load_distributions_cwd(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading distribution files in the current working directory ...")

        # Loop over the files
        for path, name in fs.files_in_path(self.config.path, extension="dat", returns=["path", "name"]):

            # Debugging
            log.debug("Loading '" + name + "' distribution ...")

            # Load the distribution
            distribution = Distribution.from_file(path)

            # Add the distribution
            if self.config.panels: self.add_distribution(distribution, label=name, panel=name)
            else: self.add_distribution(distribution, label=name)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot
        if self.has_single_panel: self.plot_single_panel()
        else: self.plot_more_panels()

    # -----------------------------------------------------------------

    def plot_single_panel(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.debug("Making a single-panel distribution plot ...")

        # Setup the figure
        self.main_plot = self.figure.create_one_plot()

        colors = iter(pretty_colours)

        # Plot the distributions
        for label in self.distributions:

            distribution = self.distributions[label]

            # Plot the distribution as a histogram
            alpha = None
            edgecolor = "black"
            linewidth = 1

            edgecolor_list = [edgecolor for _ in range(distribution.nvalues)]
            linewidth_list = [linewidth for _ in range(distribution.nvalues)]

            # Pick color
            color = next(colors)

            if self.logscale:
                self.main_plot.bar(distribution.edges_log[:-1], distribution.frequencies, width=distribution.bin_widths_log,
                        linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)
            else:
                self.main_plot.bar(distribution.edges[:-1], distribution.frequencies, width=distribution.bin_widths,
                        linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)

            # Add frequencies
            if self.config.frequencies:
                for value, frequency in zip(distribution.values, distribution.frequencies):
                    y_value = frequency + 0.02
                    self.main_plot.text(value, y_value, tostr(frequency, round=True, ndigits=2), color=color, fontweight='bold', va='center')

            min_value = distribution.min_edge
            if self.config.frequencies: max_value = distribution.max_edge + 0.02
            else: max_value = distribution.max_edge
            min_frequency = 0. if not self.logfrequency else 0.5 * distribution.min_frequency_nonzero
            max_frequency = 1.1 * distribution.max_frequency if not self.logfrequency else 2. * distribution.max_frequency

            # Keep track of minimum and maximum value
            if self._min_value is None or min_value < self._min_value: self._min_value = min_value
            if self._max_value is None or max_value > self._max_value: self._max_value = max_value

            # Keep track of the minimal and maximal count
            if self._min_frequency is None or min_frequency < self._min_frequency: self._min_frequency = min_frequency
            if self._max_frequency is None or max_frequency > self._max_frequency: self._max_frequency = max_frequency

            # Add smooth
            if self.add_smooth:

                if self.logfrequency:
                    x_smooth, y_smooth = distribution.smooth_values_log(x_min=self._min_value, x_max=self._max_value)
                    self.main_plot.plot(x_smooth, y_smooth, color='red', linewidth=1)
                else:
                    x_smooth, y_smooth = distribution.smooth_values(x_min=self._min_value, x_max=self._max_value)
                    self.main_plot.plot(x_smooth, y_smooth, color='red', linewidth=1)

            # Add extrema
            if self.add_extrema:

                # Local maxima
                x, y = distribution.local_maxima
                self.main_plot.plot(x, y, color='green', marker='^')

                # Local minima
                x, y = distribution.local_minima
                self.main_plot.plot(x, y, color='red', marker='v')

            # Add statistics
            if self.add_statistics:

                if self.logscale: mean_line = plt.axvline(distribution.geometric_mean, color="green", linestyle="dashed", label="Mean")
                else: mean_line = plt.axvline(distribution.mean, color="green", linestyle="dashed", label="Mean")
                median_line = plt.axvline(distribution.median, color="purple", linestyle="dashed", label="Median")
                max_line = plt.axvline(distribution.most_frequent, color="orange", linestyle="dashed", label="Most frequent")

        # Finish plot
        self.finish_main_plot()

    # -----------------------------------------------------------------

    @property
    def x_scale(self):

        """
        This function ...
        :return:
        """

        if self.logscale: return "log"
        else: return "linear"

    # -----------------------------------------------------------------

    @property
    def x_scales(self):

        """
        This function ...
        :return:
        """

        return [self.x_scale] * self.npanels

    # -----------------------------------------------------------------

    @property
    def y_scale(self):

        """
        This function ...
        :return:
        """

        if self.logfrequency: return "log"
        else: return "linear"

    # -----------------------------------------------------------------

    @property
    def y_label(self):

        """
        This function ...
        :return:
        """

        return "Frequency"

    # -----------------------------------------------------------------

    @property
    def x_label(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_panel: raise RuntimeError("Function shouldn't be called: multiple panels")
        return "Value"

    # -----------------------------------------------------------------

    @property
    def x_labels(self):

        """
        This function ...
        :return:
        """

        from ..tools import sequences

        labels = []
        for panel in self.panels:
            if self.has_single_distribution_panel(panel):
                distribution = self.get_single_distribution_panel(panel)
                x_label = distribution.value_name
                if distribution.unit is not None: x_label += " [" + str(distribution.unit) + "]"
                label = x_label
            else:
                distributions = self.get_distributions_panel(panel)
                value_names = [distribution.value_name for distribution in distributions]
                if not sequences.all_equal(value_names): label = value_names[0]
                else: label = "Value"

            # Add label
            labels.append(label)

        # Return labels
        return labels

    # -----------------------------------------------------------------

    def plot_more_panels(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Making a multi-panel distribution plot ...")

        width_ratios = [1.] * self.npanels

        # Create row
        plots = self.figure.create_row(self.npanels, share_axis=True, width_ratios=width_ratios, x_scales=self.x_scales, y_scale=self.y_scale, y_label=self.y_label, x_labels=self.x_labels)

        # Set panel plots
        for panel, plot in zip(self.panels, plots): self.panel_plots[panel] = plot

        colors = iter(pretty_colours)

        # Loop over the panels
        for panel in self.panels:

            # Get the plot
            plot = self.panel_plots[panel]

            # Loop over the distributions
            for label in self.distributions[panel]:

                # Get the distribution
                distribution = self.distributions[panel][label]

                # Plot the distribution as a histogram
                alpha = None
                edgecolor = "black"
                linewidth = 1

                edgecolor_list = [edgecolor for _ in range(distribution.nvalues)]
                linewidth_list = [linewidth for _ in range(distribution.nvalues)]

                # Pick color
                color = next(colors)

                if self.logscale:
                    plot.bar(distribution.edges_log[:-1], distribution.frequencies, width=distribution.bin_widths_log,
                            linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)
                else:
                    plot.bar(distribution.edges[:-1], distribution.frequencies, width=distribution.bin_widths,
                            linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)

                # Add frequencies
                if self.config.frequencies:
                    for value, frequency in zip(distribution.values, distribution.frequencies):
                        y_value = frequency + 0.02
                        plot.text(value, y_value, tostr(frequency, round=True, ndigits=2), color=color, horizontalalignment='center', size="small") # fontweight='bold'?

                # Get min and max values
                min_value = distribution.min_edge
                if self.config.frequencies: max_value = distribution.max_edge + 0.02
                else: max_value = distribution.max_edge
                min_frequency = 0. if not self.logfrequency else 0.5 * distribution.min_frequency_nonzero
                max_frequency = 1.1 * distribution.max_frequency if not self.logfrequency else 2. * distribution.max_frequency

                # Keep track of minimum and maximum value
                if self.min_values[panel] is None or min_value < self.min_values[panel]: self.min_values[panel] = min_value
                if self.max_values[panel] is None or max_value > self.max_values[panel]: self.max_values[panel] = max_value

                # Keep track of the minimal and maximal frequency
                if self._min_frequency is None or min_frequency < self._min_frequency: self._min_frequency = min_frequency
                if self._max_frequency is None or max_frequency > self._max_frequency: self._max_frequency = max_frequency

                # Add smooth
                if self.add_smooth:

                    if self.logfrequency:
                        x_smooth, y_smooth = distribution.smooth_values_log(x_min=self.min_value, x_max=self.max_value)
                        plot.plot(x_smooth, y_smooth, color='red', linewidth=1)
                    else:
                        x_smooth, y_smooth = distribution.smooth_values(x_min=self.min_value, x_max=self.max_value)
                        plot.plot(x_smooth, y_smooth, color='red', linewidth=1)

                # Add extrema
                if self.add_extrema:

                    # Local maxima
                    x, y = distribution.local_maxima
                    plot.plot(x, y, color="green", marker='^')

                    # Local minima
                    x, y = distribution.local_minima
                    plot.plot(x, y, color="red", marker='v')

                # Add statistics
                if self.add_statistics:

                    if self.logscale: mean_line = plot.axvline(distribution.geometric_mean, color="green", linestyle="dashed", label="Mean")
                    else: mean_line = plot.axvline(distribution.mean, color="green", linestyle="dashed", label="Mean")
                    median_line = plot.axvline(distribution.median, color="purple", linestyle="dashed", label="Median")
                    max_line = plot.axvline(distribution.most_frequent, color="orange", linestyle="dashed", label="Most frequent")

        # Finish plot
        self.finish_panels()

    # -----------------------------------------------------------------

    def finish_main_plot(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Finishing plot ...")

        # Axis limits are now definite
        if self.min_value is None: self.min_value = self._min_value
        if self.max_value is None: self.max_value = self._max_value
        if self.min_frequency is None: self.min_frequency = self._min_frequency
        if self.max_frequency is None: self.max_frequency = self._max_frequency

        # Set axes limits
        self.main_plot.set_xlim((self.min_value, self.max_value))
        self.main_plot.set_ylim((self.min_frequency, self.max_frequency))

        # Add legend
        if len(self.distributions) > 1 and self.legend: plt.legend()

        # Set log scales
        if self.logscale: self.main_plot.set_xscale("log")
        if self.logfrequency: self.main_plot.set_yscale("log")

        # Ticks
        self.main_plot.set_xticks(fontsize=self.config.plot.ticks_fontsize)
        self.main_plot.set_yticks(fontsize=self.config.plot.ticks_fontsize)

        # Labels
        self.main_plot.set_xlabel(self.x_label, fontsize='small')
        self.main_plot.set_ylabel(self.y_label, fontsize='small')

        # Set grid
        self.figure.set_grid(self.config.plot, which="both")

        # Set borders
        self.figure.set_borders(self.config.plot)

        # Add title if requested
        if self.title is not None: self.figure.set_title(self.title)

        # Save or show the plot
        if self.config.show: self.figure.show()

        # Save the figure
        if self.out_path is not None: self.save_figure()

    # -----------------------------------------------------------------

    @property
    def first_panel(self):

        """
        This function ...
        :return:
        """

        return self.panels[0]

    # -----------------------------------------------------------------

    @property
    def first_panel_plot(self):

        """
        This function ....
        :return:
        """

        return self.panel_plots[self.first_panel]

    # -----------------------------------------------------------------

    def finish_panels(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Finishing plot with panels ...")

        # Set y limits
        self.first_panel_plot.set_ylim(self.min_frequency, self.max_frequency)

        # Set x limits
        for panel in self.panels:
            plot = self.panel_plots[panel]
            plot.set_xlim((self.min_values[panel], self.max_values[panel]))

        # Set grid
        #self.figure.set_grid(self.config.plot, which="both")

        # Set borders
        #self.figure.set_borders(self.config.plot)

        # Add title if requested
        if self.title is not None: self.figure.set_title(self.title)

        # Save or show the plot
        if self.config.show: self.figure.show()

        # Save the figure
        if self.out_path is not None: self.save_figure()

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Close the figure
        if self.figure is not None: self.figure.close()

    # -----------------------------------------------------------------

    def save_figure(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Saving figure ...")

        # Determine path
        if types.is_string_type(self.out_path):
            if fs.is_directory(self.out_path): path = fs.join(self.out_path, "seds" + self.config.format)
            else: path = self.out_path
        else: path = self.out_path

        # Save
        self.figure.saveto(path)

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
