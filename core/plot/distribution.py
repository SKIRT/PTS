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
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from collections import OrderedDict
from textwrap import wrap
from matplotlib import rc
import matplotlib.gridspec as gridspec
from matplotlib.ticker import LinearLocator, LogLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import get_cmap

# Import the relevant PTS classes and modules
from ..basics.log import log
from ...magic.tools.plotting import pretty_colours
from ..tools import strings
from ..basics.configurable import Configurable
from ..basics.containers import DefaultOrderedDict
from ..tools import filesystem as fs
from ..basics.distribution import Distribution, Distribution2D
from ..tools import types
from ..basics.plot import MPLFigure, BokehFigure, BokehPlot, mpl, bokeh
from ..tools.stringify import tostr

# -----------------------------------------------------------------

# Use LaTeX rendering
rc('text', usetex=True)

# -----------------------------------------------------------------

#patterns = ('-', '+', 'x', '\\', '*', 'o', 'O', '.')
patterns = ("/", "\\", '-', '+', 'x', "//", '*', 'o', 'O', '.')
line_styles = ['-', '--', '-.', ':']

# -----------------------------------------------------------------

def plot_distribution_from_file(path, **kwargs):

    """
    This function ...
    :param path:
    :param kwargs:
    :return:
    """

    # Load
    distr = Distribution.from_file(path)

    # Plot
    plot_distribution(distr, **kwargs)

# -----------------------------------------------------------------

def plot_distribution(distribution, path=None, logscale=False, logfrequency=False, title=None, x_limits=None,
                      y_limits=None, color="xkcd:sky blue", x_label=None, y_label=None, format="pdf", axes=None,
                      xsize=5, ysize=5, colors=None, statistics=True, soft_xmin=False,
                      soft_xmax=False, soft_ymin=False, soft_ymax=False, mean=None, median=None, most_frequent=None,
                      stddev=None, fwhm=None, show_mean=None, show_median=None, show_most_frequent=None,
                      show_stddev=False, show_fwhm=False, alpha=None, plot=None, cmap=None, cmap_interval=None,
                      mean_color="green", median_color="purple", most_frequent_color="orange"):

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
    :param format:
    :param axes:
    :param xsize:
    :param ysize:
    :param colors:
    :param statistics:
    :param soft_xmin:
    :param soft_xmax:
    :param soft_ymin:
    :param soft_ymax:
    :param mean:
    :param median:
    :param most_frequent:
    :param stddev:
    :param fwhm:
    :param show_mean:
    :param show_median:
    :param show_most_frequent:
    :param show_stddev:
    :param show_fwhm:
    :param alpha:
    :param plot:
    :param cmap:
    :param cmap_interval:
    :param mean_color:
    :param median_color:
    :param most_frequent_color:
    :return:
    """

    # Output
    from ...core.basics.map import Map
    output = Map()

    # Get format from path
    if path is not None: format = fs.get_extension(path)

    # Set color (or sequences of colors)
    if colors is not None: color = colors

    # Create figure if necessary, get the axes
    if plot is not None: axes = plot.axes
    only_axes = False
    if axes is None:
        figure = plt.figure(figsize=(xsize, ysize))
        axes = plt.gca()
        rect = figure.patch
        rect.set_facecolor('white')
    else: only_axes = True

    edgecolor = "black"
    linewidth = 1

    edgecolor_list = [edgecolor for _ in range(distribution.nvalues)]
    linewidth_list = [linewidth for _ in range(distribution.nvalues)]

    # Plot the bar graph
    if logscale:

        if cmap is not None:

            if cmap_interval is None: norm = LogNorm(vmin=distribution.edges_log[0], vmax=distribution.edges_log[-1])
            else: norm = LogNorm(vmin=cmap_interval[0], vmax=cmap_interval[1])

            # Get colours
            cmap = get_cmap(cmap)
            color = [cmap(norm(value)) for value in distribution.values]

        # Plot
        patch = axes.bar(distribution.edges_log[:-1], distribution.frequencies, width=distribution.bin_widths_log, linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)
        output.patch = patch

    else:

        if cmap is not None:

            if cmap_interval is None: norm = Normalize(vmin=distribution.edges[0], vmax=distribution.edges[-1])
            else: norm = Normalize(vmin=cmap_interval[0], vmax=cmap_interval[1])

            # Get colours
            cmap = get_cmap(cmap)
            color = [cmap(norm(value)) for value in distribution.values]

        # Plot
        patch = axes.bar(distribution.edges[:-1], distribution.frequencies, width=distribution.bin_widths, linewidth=linewidth_list, alpha=alpha, align="edge", color=color, edgecolor=edgecolor_list)
        output.patch = patch

    # Determine automatic x_min and x_max
    x_min_auto = distribution.min_edge_log if logscale else distribution.min_edge
    x_max_auto = distribution.max_edge_log if logscale else distribution.max_edge

    # Use automatic xmin and xmax
    if x_limits is None: x_min, x_max = x_min_auto, x_max_auto

    # From arguments
    else:

        if x_limits[0] is None: x_min = x_min_auto
        elif soft_xmin: x_min = max(x_limits[0], x_min_auto)
        else: x_min = x_limits[0]

        if x_limits[1] is None: x_max = x_max_auto
        elif soft_xmax: x_max = min(x_limits[1], x_max_auto)
        else: x_max = x_limits[1]

    # Calculate automatic y_min and y_max
    y_min_auto = 0. if not logfrequency else 0.5 * distribution.min_frequency_nonzero
    y_max_auto = 1.1 * distribution.max_frequency if not logfrequency else 2. * distribution.max_frequency

    # Use automatic ymin and ymax
    if y_limits is None: y_min, y_max = y_min_auto, y_max_auto

    # From arguments
    else:

        if y_limits[0] is None: y_min = y_min_auto
        elif soft_ymin: y_min = max(y_limits[0], y_min_auto)
        else: y_min = y_limits[0]

        if y_limits[1] is None: y_max = x_max_auto
        elif soft_ymax: y_max = min(y_limits[1], y_max_auto)
        else: y_max = y_limits[1]

    #print(x_min, x_max)
    #print(y_min, y_max)

    # Set the axis limits
    axes.set_xlim(x_min, x_max)
    axes.set_ylim(y_min, y_max)

    # Set flags
    if show_mean is None: show_mean = statistics
    if show_median is None: show_median = statistics
    if show_most_frequent is None: show_most_frequent = statistics
    if show_stddev is None: show_stddev = statistics
    if show_fwhm is None: show_fwhm = statistics

    # Show mean
    if show_mean:

        # Get mean
        if mean is None:
            if logscale: mean = distribution.geometric_mean_value
            else: mean = distribution.mean_value

        # Show
        mean_line = axes.axvline(mean, color=mean_color, linestyle="dashed", label="Mean")

    # Don't show mean
    else: mean_line = None

    # Show median
    if show_median:

        # Get median
        if median is None: median = distribution.median_value

        # Show
        median_line = axes.axvline(median, color=median_color, linestyle="dashed", label="Median")

    # Don't show median
    else: median_line = None

    # Show most frequent
    if show_most_frequent:

        # Get most frequent
        if most_frequent is None: most_frequent = distribution.most_frequent_value

        # Show
        most_frequent_line = axes.axvline(most_frequent, color=most_frequent_color, linestyle="dashed", label="Most frequent")

    # Don't show most frequent
    else: most_frequent_line = None

    # Show stddev
    if show_stddev:

        # Get the stddev
        if stddev is None: stddev = distribution.stddev_value

    # Show FWHM
    if show_fwhm:

        # Get the FWHM
        if fwhm is None: fwhm = distribution.fwhm_value

    # Add legend
    axes.legend()

    # Axes were not provided, but figure was created here
    if not only_axes:

        # Create legend
        #plt.legend()

        if x_label is None: x_label = distribution.value_name
        if y_label is None: y_label = distribution.y_name
        #print("UNIT:", distribution.unit)
        if distribution.unit is not None: x_label += " [" + str(distribution.unit) + "]"

        # Put the title and labels
        if title is not None: axes.set_title(strings.split_in_lines(title))
        axes.set_xlabel(x_label)
        axes.set_ylabel(y_label)

        if logfrequency: axes.set_yscale("log", nonposx='clip')
        if logscale: axes.set_xscale("log")

        # Set locator
        if logscale:

            # Set locator
            log_subs = (1., 2., 5.)
            major_locator = LogLocator(subs=log_subs)
            # minor_locator = None

            xaxis = axes.get_xaxis()

            # Set the locators
            xaxis.set_major_locator(major_locator)

            # Set formatter
            major_formatter = FormatStrFormatter('%g')
            xaxis.set_major_formatter(major_formatter)

        plt.tight_layout()
        #plt.grid(alpha=0.8)

        if path is None: plt.show()
        else: figure.savefig(path, format=format)

        # Close the figure
        plt.close()

    # Return
    return output

# -----------------------------------------------------------------

def plot_distributions(distributions, panels=False, smooth=False, statistics=False, extrema=False,
                       maxima=None, minima=None, edges=False, frequencies=False, path=None, logscale=False,
                       logfrequency=False, alpha=None):

    """
    This function ...
    :param distributions:
    :param panels:
    :param smooth:
    :param statistics:
    :param extrema:
    :param maxima:
    :param minima:
    :param edges:
    :param frequencies:
    :param path:
    :param logscale:
    :param logfrequency:
    :param alpha:
    :return:
    """

    # Check input
    ndistributions = len(distributions)
    if ndistributions == 0: raise ValueError("No distributions are passed")

    # Initialize plotter
    plotter = DistributionPlotter()

    # Set options
    plotter.config.smooth = smooth
    plotter.config.statistics = statistics
    plotter.config.extrema = extrema
    plotter.config.maxima = maxima
    plotter.config.minima = minima
    plotter.config.edges = edges
    plotter.config.frequencies = frequencies
    plotter.config.alpha = alpha

    # Logscales
    plotter.config.logscale = logscale
    plotter.config.logfrequency = logfrequency

    # Add distributions
    for label in distributions:

        # Get the distribution
        distribution = distributions[label]

        # Add the distribution
        if panels: plotter.add_distribution(distribution, label, panel=label)
        else: plotter.add_distribution(distribution, label)

    # Run the plotter
    plotter.run(output=path)

# -----------------------------------------------------------------

def plot_2d_distribution_from_file(path, **kwargs):

    """
    This function ...
    :param path:
    :param kwargs:
    :return:
    """

    # Load
    distr = Distribution2D.from_file(path)

    # Plot
    plot_2d_distribution(distr, **kwargs)

# -----------------------------------------------------------------

def plot_2d_distribution(distribution, title=None, path=None, x_lines=None, y_lines=None, xlim=None, ylim=None,
                         mode="colormesh", add_average=False, format=None, cmap="Blues"):

    """
    This function ...
    :param distribution:
    :param title:
    :param path:
    :param x_lines:
    :param y_lines:
    :param xlim:
    :param ylim:
    :param mode:
    :param add_average:
    :return:
    """

    # Get format from path
    if path is not None: format = fs.get_extension(path)

    # Create a figure
    fig = plt.figure()

    #ax.set_ylabel('$\mathcal{F}_\mathrm{unev.}^\mathrm{abs}$', fontsize=18)
    #ax.set_xlabel('R (kpc)', fontsize=18)
    # ax.hexbin(r/1000.,F_abs_yng,gridsize=150,bins='log',cmap=plt.cm.autumn, mincnt=1,linewidths=0)
    #print("x edges:", self.x_edges, self.x_edges.shape)
    #print("y edges:", self.y_edges, self.y_edges.shape)
    #print("counts:", self.counts, self.counts.shape)
    #ax.pcolor(self.x_edges, self.y_edges, self.counts)
    #ax.imshow(self.counts)

    if mode == "colormesh":

        ax = fig.gca()

        # Plot the colored 2D distribution
        # Plot 2D histogram using pcolor
        ax.pcolormesh(distribution.x, distribution.y, distribution.counts, cmap=cmap)

        # Plot the running average
        # if self.rBins_F is not None and self.FBins_r is not None:
        #    #ax.plot(self.rBins_F, self.FBins_r, 'k-', linewidth=2)
        #    ax.plot(self.rBins_F, self.FBins_r, 'w-', linewidth=1)

        # ax.errorbar(1.7, 0.88, xerr=1.4, color='k')
        # ax.text(1.8, 0.90, 'Bulge', ha='center')
        # ax.errorbar(11., 0.88, xerr=2.75, color='k')
        # ax.text(11., 0.90, 'main SF ring', ha='center')
        # ax.errorbar(16., 0.88, xerr=1, color='k')
        # ax.text(15., 0.90, r'$2^\mathrm{nd}$ SF ring', ha='left')

        # Plot vertical lines
        if x_lines is not None:
            for x in x_lines: plt.axvline(x=x)

        # Plot horizontal lines
        if y_lines is not None:
            for y in y_lines: plt.axhline(y=y)

        # Set axes limits
        if xlim is not None: ax.set_xlim(xlim[0], xlim[1])
        if ylim is not None: ax.set_ylim(ylim[0], ylim[1])

        # Set labels
        ax.set_xlabel(distribution.x_name)
        ax.set_ylabel(distribution.y_name)

    elif mode == "3d":

        #fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Construct arrays for the anchor positions of the bars.
        # Note: np.meshgrid gives arrays in (ny, nx) so we use 'F' to flatten xpos,
        # ypos in column-major order. For numpy >= 1.7, we could instead call meshgrid
        # with indexing='ij'.
        xpos, ypos = np.meshgrid(distribution.x[:-1] + 0.25, distribution.y[:-1] + 0.25)
        xpos = xpos.flatten('F')
        ypos = ypos.flatten('F')
        zpos = np.zeros_like(xpos)

        # Construct arrays with the dimensions for the 16 bars.
        dx = 0.5 * np.ones_like(zpos)
        dy = dx.copy()
        dz = distribution.counts.flatten()
        #dz = distribution.counts

        # Plot 3D bars
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

        #plt.show()

    # Set plot title
    if title is not None: ax.set_title(title)

    # Save the figure
    if path is not None: plt.savefig(path, format=format)
    else: plt.show()

    # Close
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

standard_panel = "standard"

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

        # The plotting properties for the distributions
        self.properties = DefaultOrderedDict(OrderedDict)

        # Panel properties
        self.panel_properties = OrderedDict()

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
        self.add_minima = False
        self.add_maxima = False
        self.add_statistics = False
        self.add_hatches = False
        self._ylabel = None

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

    @property
    def add_lower(self):
        return self.config.lower_than is not None

    # -----------------------------------------------------------------

    @property
    def lower_than(self):
        return self.config.lower_than

    # -----------------------------------------------------------------

    @property
    def add_higher(self):
        return self.config.higher_than is not None

    # -----------------------------------------------------------------

    @property
    def higher_than(self):
        return self.config.higher_than

    # -----------------------------------------------------------------

    def add_distribution(self, distribution, label, panel=standard_panel, properties=None):

        """
        This function ...
        :param distribution:
        :param label:
        :param panel:
        :param properties:
        :return:
        """

        # Set distribution
        self.distributions[panel][label] = distribution

        # Set properties
        if properties is not None: self.properties[panel][label] = properties

    # -----------------------------------------------------------------

    def set_panel_properties(self, label, properties):

        """
        This function ...
        :param label:
        :param properties:
        :return:
        """

        if label in self.panel_properties: self.panel_properties[label].update(properties)
        else: self.panel_properties[label] = properties

    # -----------------------------------------------------------------

    @property
    def npanels(self):
        return len(self.distributions)

    # -----------------------------------------------------------------

    @property
    def panels(self):
        return self.distributions.keys()

    # -----------------------------------------------------------------

    @property
    def has_single_panel(self):
        return self.npanels == 1

    # -----------------------------------------------------------------

    @property
    def single_panel(self):
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
        total = 0
        for panel in self.panels: total += self.ndistributions_for_panel(panel)
        return total

    # -----------------------------------------------------------------

    @property
    def has_single_distribution(self):
        return self.ndistributions == 1

    # -----------------------------------------------------------------

    @property
    def single_distribution_label(self):
        if not self.has_single_distribution: raise ValueError("Not a single distribution")
        return self.distributions[self.single_panel].keys()[0]

    # -----------------------------------------------------------------

    @property
    def single_distribution(self):
        return self.distributions[self.single_panel][self.single_distribution_label]

    # -----------------------------------------------------------------

    @property
    def no_distributions(self):
        return self.ndistributions == 0

    # -----------------------------------------------------------------

    def get_distribution_labels_panel(self, panel):

        """
        This funtion ...
        :param panel:
        :return:
        """

        return self.distributions[panel].keys()

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

    @property
    def do_normalize(self):
        return self.config.normalize is not None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Normalize?
        if self.do_normalize: self.normalize()

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
        self.add_minima = False
        self.add_maxima = False
        self.add_statistics = False
        self.add_hatches = False
        self._ylabel = None

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
        if "legend" in kwargs: self.legend = kwargs.pop("legend")

        # Logscale
        self.logscale = kwargs.pop("logscale", self.config.logscale)
        self.logfrequency = kwargs.pop("logfrequency", self.config.logfrequency)

        # Features
        self.add_smooth = kwargs.pop("add_smooth", self.config.smooth)
        self.add_minima = kwargs.pop("add_minima", self.config.minima)
        self.add_maxima = kwargs.pop("add_maxima", self.config.maxima)
        if kwargs.get("add_extrema", None) is not None: self.add_minima = self.add_maxima = kwargs.pop("add_extrema")
        self.add_statistics = kwargs.pop("add_statistics", self.config.statistics)
        self.add_hatches = kwargs.pop("add_hatches", self.config.hatches)
        self._ylabel = kwargs.pop("y_label", self.config.y_label)

        # Set the output path
        self.out_path = kwargs.pop("output", None)
        if self.out_path is None and "output" in self.config and self.config.output is not None:
            full_output_path = fs.absolute_or_in(self.config.output, self.config.path)
            if fs.has_extension(full_output_path):
                directory_path = fs.directory_of(full_output_path)
                if not fs.is_directory(directory_path): fs.create_directory(directory_path)
            elif not fs.is_directory(full_output_path): fs.create_directory(full_output_path)
            self.out_path = full_output_path

        # Add distribution files present in the current working directory (if nothing is added manually)
        if self.no_distributions: self.load_distributions()

        # Initialize min and max values for panels
        for panel in self.panels: self.min_values[panel] = None
        for panel in self.panels: self.max_values[panel] = None

        # Create the plot
        if self.config.library == mpl: self.figure = MPLFigure(size=self.figsize)
        elif self.config.library == bokeh: self.figure = BokehFigure()
        else: raise ValueError("Invalid libary: " + self.config.library)

        # Set the 'show' flag
        if self.config.show is None:
            if self.out_path is not None: self.config.show = False
            else: self.config.show = True

    # -----------------------------------------------------------------

    @property
    def xsize(self):
        if self.has_single_panel: return self.config.plot.xsize
        else: return self.config.plot.xsize * self.npanels

    # -----------------------------------------------------------------

    @property
    def ysize(self):
        return self.config.plot.ysize

    # -----------------------------------------------------------------

    @property
    def figsize(self):
        return (self.xsize, self.ysize)

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
        for path, name in fs.files_in_path(self.config.path, extension="dat", returns=["path", "name"], recursive=self.config.recursive):

            # Check whether actually a distribution
            if fs.get_ncolumns(path) != 2:
                log.warning("Assuming the file '" + name + "' is not a distribution, skipping ...")
                continue

            # Debugging
            log.debug("Loading '" + name + "' distribution ...")

            # Load the distribution
            distribution = Distribution.from_file(path)

            # Add the distribution
            if self.config.panels:

                current_labels = self.get_distribution_labels_panel(name)
                if name in current_labels:
                    log.warning("A distribution with the name '" + name + "' is already added for panel '" + name + "', adjusting distribution label ...")
                    label = name + str(len(current_labels))
                else: label = name

                # Add
                self.add_distribution(distribution, label=label, panel=name)

            else:
                current_labels = self.get_distribution_labels_panel(standard_panel)
                if name in current_labels:
                    log.warning("A distribution with the name '" + name + "' is already added, adjusting distribution label ...")
                    label = name + str(len(current_labels))
                else: label = name

                # Ad
                self.add_distribution(distribution, label=label)

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the distributions ...")

        # Loop over the panels
        for panel in self.panels:

            # Loop over the distributions
            for label in self.distributions[panel]:

                # Get the distribution
                distribution = self.distributions[panel][label]

                # Normalize and replace
                self.distributions[panel][label] = distribution.normalized(value=self.config.normalization_value, method=self.config.normalize)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Single panel
        if self.has_single_panel: self.plot_single_panel()

        # Multiple panels
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

        all_bars = []

        # Plot the distributions
        for label in self.distributions[self.single_panel]:

            # Get the distribution
            distribution = self.distributions[self.single_panel][label]

            # Get the properties
            properties = self.properties[self.single_panel][label] if self.single_panel in self.properties and label in self.properties[self.single_panel] else None

            # Plot the distribution as a histogram
            #alpha = None
            edgecolor = "black"
            linewidth = 1

            edgecolor_list = [edgecolor for _ in range(distribution.nvalues)]
            linewidth_list = [linewidth for _ in range(distribution.nvalues)]

            # Pick color
            color = next(colors)

            #print(distribution.edges_log)
            #print(distribution.frequencies)
            #print(distribution.bin_widths_log)

            if self.logscale:
                #bars = self.main_plot.bar(distribution.edges_log[:-1], distribution.frequencies, width=distribution.bin_widths_log,
                #        linewidth=linewidth_list, alpha=self.config.alpha, align="edge", color=color, edgecolor=edgecolor_list, label=label)
                widths = np.array(distribution.bin_widths_log)
                bars = self.main_plot.bar(distribution.edges_log[:-1], distribution.frequencies,
                                          width=widths,
                                          linewidth=linewidth_list, alpha=self.config.alpha, align="center", color=color,
                                          edgecolor=edgecolor_list, label=label)
            else:
                #bars = self.main_plot.bar(distribution.edges[:-1], distribution.frequencies, width=distribution.bin_widths,
                #        linewidth=linewidth_list, alpha=self.config.alpha, align="edge", color=color, edgecolor=edgecolor_list, label=label)
                widths = np.array(distribution.bin_widths) * self.config.bar_width
                bars = self.main_plot.bar(distribution.edges[:-1], distribution.frequencies,
                                          width=widths,
                                          linewidth=linewidth_list, alpha=self.config.alpha, align="center", color=color,
                                          edgecolor=edgecolor_list, label=label)
            all_bars.append(bars)

            # Add frequencies
            if self.config.frequencies:
                for value, frequency in zip(distribution.values, distribution.frequencies):
                    y_value = frequency + 0.02
                    self.main_plot.text(value, y_value, tostr(frequency, round=True, ndigits=2), color=color, fontweight='bold', va='center')

            # Get min and max value
            min_value = distribution.min_edge
            max_value = distribution.max_edge

            # Set min and max frequency
            min_frequency = 0. if not self.logfrequency else 0.5 * distribution.min_frequency_nonzero
            max_frequency = 1.1 * distribution.max_frequency if not self.logfrequency else 2. * distribution.max_frequency
            #if self.config.frequencies: max_frequency += 0.05

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

            # Add maxima
            if self.add_maxima:

                # Local maxima
                x, y = distribution.local_maxima
                self.main_plot.scatter(x, y, color='green', marker='^', s=100)

            # Add minima
            if self.add_minima:

                # Local minima
                x, y = distribution.local_minima
                self.main_plot.scatter(x, y, color='red', marker='v', s=100)

            # Add lower than
            if self.add_lower:

                from ...magic.tools.plotting import align_marker

                # Get
                x, y = distribution.less_frequent_than(self.lower_than)
                self.main_plot.scatter(x, y, color=color, marker=align_marker('v', valign="bottom"), s=100)
                if self.config.add_lower_frequencies:
                    for xi, yi in zip(x, y):
                        y_value = yi + 0.05
                        if self.config.lower_frequency_percentages: text = tostr(yi * 100, round=True, decimal_places=1) + "%"
                        else: text = tostr(yi, round=True, ndigits=2)
                        self.main_plot.text(xi, y_value, text, color=color, horizontalalignment='center', size="small")

            # Add higher than
            if self.add_higher:

                from ...magic.tools.plotting import align_marker

                # Get
                x, y = distribution.more_frequent_than(self.higher_than)
                self.main_plot.scatter(x, y, color=color, marker=align_marker('^', valign="bottom"), s=100)
                if self.config.add_higher_frequencies:
                    for xi, yi in zip(x, y):
                        y_value = yi + 0.05
                        if self.config.higher_frequency_percentages: text = tostr(yi * 100, round=True, decimal_places=1) + "%"
                        else: text = tostr(yi, round=True, ndigits=2)
                        self.main_plot.text(xi, y_value, text, color=color, horizontalalignment='center', size="small")

            # Add statistics
            if self.add_statistics:

                if self.logscale: mean_line = plt.axvline(distribution.geometric_mean, color="green", linestyle="dashed", label="Mean")
                else: mean_line = plt.axvline(distribution.mean, color="green", linestyle="dashed", label="Mean")
                median_line = plt.axvline(distribution.median, color="purple", linestyle="dashed", label="Median")
                max_line = plt.axvline(distribution.most_frequent, color="orange", linestyle="dashed", label="Most frequent")

            # Add vlines?
            if properties is not None and properties.vlines is not None:

                # Sequence
                if types.is_sequence(properties.vlines):
                    styles = iter(line_styles)
                    for value in properties.vlines:
                        linestyle = styles.next()
                        plt.axvline(value, color="black", linestyle=linestyle)

                # Dictionary
                elif types.is_dictionary(properties.vlines):
                    styles = iter(line_styles)
                    for line_label in properties.vlines:
                        linestyle = styles.next()
                        plt.axvline(properties.vlines[line_label], color="black", linestyle=linestyle, label=line_label)

        # Set hatches
        if self.add_hatches:
            for bars, pattern in zip(all_bars, patterns):
                for bar in bars: bar.set_hatch(pattern)

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

        if self._ylabel is not None: return self._ylabel

        from ..tools import sequences

        labels = []
        for panel in self.panels:
            distributions = self.get_distributions_panel(panel)
            panel_labels = [distribution.y_name for distribution in distributions]
            if not sequences.all_equal(panel_labels): raise ValueError("Y labels of panel distributions are not equal")
            panel_label = panel_labels[0]
            labels.append(panel_label)
        if not sequences.all_equal(labels): raise ValueError("Y labels of panels are not equal")
        return labels[0].replace("_", "\_")

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

        #if self._xlabels is not None: return self._xlabels

        from ..tools import sequences

        labels = []
        for panel in self.panels:

            # Panel label is defined in panel properties
            if panel in self.panel_properties and "x_label" in self.panel_properties[panel]: label = self.panel_properties[panel].x_label

            # Single distribution for this panel
            elif self.has_single_distribution_panel(panel):

                distribution = self.get_single_distribution_panel(panel)

                # Use the name of the distribution?
                if self.config.use_name_xlabel: x_label = self.get_single_distribution_label_panel(panel)

                # Use the label of the x axis of the distribution object
                else: x_label = distribution.value_name

                # Add unit
                if distribution.unit is not None: x_label += " [" + str(distribution.unit) + "]"
                label = x_label

            # Multiple distributions for this panel
            else:

                distributions = self.get_distributions_panel(panel)
                value_names = [distribution.value_name for distribution in distributions]
                value_units = [distribution.unit for distribution in distributions]

                # Determine name
                if self.config.use_name_xlabel:
                    distribution_labels = self.get_distribution_labels_panel(panel)
                    label = strings.common_part(*distribution_labels)
                    if label is None: label = tostr(distribution_labels)
                else:
                    if sequences.all_equal(value_names): label = value_names[0]
                    else: label = "Value"

                # Determine unit
                if not sequences.all_equal(value_units): raise ValueError("Value units are different")
                value_unit = sequences.get_first_not_none_value(value_units)
                if value_unit is not None: label += " [" + str(value_unit) + "]"

            # Specify magnitude?
            if panel in self.panel_properties and "magnitude" in self.panel_properties[panel]:
                magnitude = self.panel_properties[panel].magnitude
                label += " ($\\times 10^{" + str(magnitude) + "}$)"

            # Add label
            labels.append(label.replace("_", "\_"))

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

        # Lines
        lines = OrderedDict()

        # Unique values for the x axes, for the different panels
        unique_values = defaultdict(set)

        # Loop over the panels
        for panel in self.panels:

            # Re-iterate the colours for each panel
            if self.config.colours_per_panel: colors = iter(pretty_colours)

            lines_panel = OrderedDict()

            all_bars = []

            # Get the plot
            plot = self.panel_plots[panel]

            # Loop over the distributions
            for label in self.distributions[panel]:

                # Get the distribution
                distribution = self.distributions[panel][label]

                # Get the properties
                properties = self.properties[panel][label] if panel in self.properties and label in self.properties[panel] else None

                # Plot the distribution as a histogram
                edgecolor = "black"
                linewidth = 1

                edgecolor_list = [edgecolor for _ in range(distribution.nvalues)]
                linewidth_list = [linewidth for _ in range(distribution.nvalues)]

                # Pick color
                color = next(colors)

                #print("Alpha", self.config.alpha)

                # Plot
                if self.logscale:
                    #bars = plot.bar(distribution.edges_log[:-1], distribution.frequencies, width=distribution.bin_widths_log,
                    #        linewidth=linewidth_list, alpha=self.config.alpha, align="edge", color=color, edgecolor=edgecolor_list, label=label)
                    widths = np.array(distribution.bin_widths_log) * self.config.bar_width
                    bars = plot.bar(distribution.values, distribution.frequencies,
                                    width=widths,
                                    linewidth=linewidth_list, alpha=self.config.alpha, align="center", color=color,
                                    edgecolor=edgecolor_list, label=label)
                else:
                    #bars = plot.bar(distribution.edges[:-1], distribution.frequencies, width=distribution.bin_widths,
                    #        linewidth=linewidth_list, alpha=self.config.alpha, align="edge", color=color, edgecolor=edgecolor_list, label=label)
                    widths = np.array(distribution.bin_widths)
                    bars = plot.bar(distribution.values, distribution.frequencies, width=widths,
                                    linewidth=linewidth_list, alpha=self.config.alpha, align="center", color=color,
                                    edgecolor=edgecolor_list, label=label)
                all_bars.append(bars)
                for value in distribution.values: unique_values[panel].add(value)

                # Add frequencies
                if self.config.frequencies:
                    for value, frequency in zip(distribution.values, distribution.frequencies):
                        #print(value, frequency)
                        y_value = frequency + 0.02
                        plot.text(value, y_value, tostr(frequency, round=True, ndigits=2), color=color, horizontalalignment='center', size="small") # fontweight='bold'?

                # Get min and max values
                min_value = distribution.min_edge
                max_value = distribution.max_edge

                # Get min and max frequency
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

                # Add maxima
                if self.add_maxima:

                    # Local maxima
                    x, y = distribution.local_maxima
                    plot.scatter(x, y, color="green", marker='^', s=100)

                # Add minima
                if self.add_minima:

                    # Local minima
                    x, y = distribution.local_minima
                    plot.scatter(x, y, color="red", marker='v', s=100)

                # Add lower than
                if self.add_lower:

                    from ...magic.tools.plotting import align_marker

                    # Get
                    x, y = distribution.less_frequent_than(self.lower_than)
                    plot.scatter(x, y, color=color, marker=align_marker('v', valign='bottom'), s=100)
                    if self.config.add_lower_frequencies:
                        for xi, yi in zip(x, y):
                            y_value = yi + 0.05
                            if self.config.lower_frequency_percentages: text = tostr(yi * 100, round=True, decimal_places=1) + "%"
                            else: text = tostr(yi, round=True, ndigits=2)
                            plot.text(xi, y_value, text, color=color, horizontalalignment='center', size="small")

                # Add higher than
                if self.add_higher:

                    from ...magic.tools.plotting import align_marker

                    # Get
                    x, y = distribution.more_frequent_than(self.higher_than)
                    plot.scatter(x, y, color=color, marker=align_marker('^', valign='bottom'), s=100)
                    if self.config.add_higher_frequencies:
                        for xi, yi in zip(x, y):
                            y_value = yi + 0.05
                            if self.config.higher_frequency_percentages: text = tostr(yi * 100, round=True, decimal_places=1) + "%"
                            else: text = tostr(yi, round=True, ndigits=2)
                            plot.text(xi, y_value, text, color=color, horizontalalignment='center', size="small")

                # Add statistics
                if self.add_statistics:

                    if self.logscale: mean_line = plot.axvline(distribution.geometric_mean, color="green", linestyle="dashed", label="Mean")
                    else: mean_line = plot.axvline(distribution.mean, color="green", linestyle="dashed", label="Mean")
                    median_line = plot.axvline(distribution.median, color="purple", linestyle="dashed", label="Median")
                    max_line = plot.axvline(distribution.most_frequent, color="orange", linestyle="dashed", label="Most frequent")

                    lines_panel["mean"] = mean_line
                    lines_panel["median"] = median_line
                    lines_panel["max"] = max_line

                # Add vlines?
                if properties is not None and properties.vlines is not None:

                    # Sequence
                    if types.is_sequence(properties.vlines):
                        styles = iter(line_styles)
                        for value in properties.vlines:
                            linestyle = styles.next()
                            plot.axvline(value, color="black", linestyle=linestyle)

                    # Dictionary
                    elif types.is_dictionary(properties.vlines):
                        styles = iter(line_styles)
                        for line_label in properties.vlines:
                            linestyle = styles.next()
                            line = plot.axvline(properties.vlines[line_label], color="black", linestyle=linestyle, label=line_label)
                            lines_panel[line_label] = line
                    #plot.legend()

            # Add lines
            lines[panel] = lines_panel

            # Set hatches
            if self.add_hatches:
                for bars, pattern in zip(all_bars, patterns):
                    for bar in bars: bar.set_hatch(pattern)

        # Set x axes values
        if self.config.distribution_ticks:
            x_ticks = dict()
            for panel in unique_values: x_ticks[panel] = list(sorted(unique_values[panel]))
        else: x_ticks = None

        # Finish plot
        self.finish_panels(legend_patches=lines, x_ticks=x_ticks)

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
        #print(self.legend)
        #print(len(self.distributions))
        #if len(self.distributions) > 1 and self.legend:
        if self.ndistributions > 1 and self.legend:
            log.debug("Adding legend ...")
            plt.legend()

        # Set log scales
        if self.logscale: self.main_plot.set_xscale("log")
        if self.logfrequency: self.main_plot.set_yscale("log")

        # Ticks
        self.main_plot.set_xticks(fontsize=self.config.plot.ticks_fontsize)
        if self.config.y_ticks: self.main_plot.set_yticks(fontsize=self.config.plot.ticks_fontsize)
        else: self.main_plot.hide_yticks()

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

    def finish_panels(self, legend_patches=None, x_ticks=None):

        """
        This function ...
        :param legend_patches:
        :param x_ticks:
        :return:

        """

        # Debugging
        log.debug("Finishing plot with panels ...")

        # Set y limits
        if self.min_frequency is None: self.min_frequency = self._min_frequency
        if self.max_frequency is None: self.max_frequency = self._max_frequency
        #print(self.min_frequency, self.max_frequency)
        self.first_panel_plot.set_ylim(self.min_frequency, self.max_frequency)

        # Set x limits
        for panel in self.panels:
            plot = self.panel_plots[panel]
            plot.set_xlim((self.min_values[panel], self.max_values[panel]))

        # Set x ticks
        for panel in self.panels:

            plot = self.panel_plots[panel]

            # Get specific ticks
            ticks = x_ticks[panel] if x_ticks is not None else None

            # Get magnitude
            if panel in self.panel_properties and "magnitude" in self.panel_properties[panel]: magnitude = self.panel_properties[panel].magnitude
            else: magnitude = None

            # Set ticks
            # plot.set_xticks will set things nice also for logscales
            plot.set_xticks(ticks, magnitude=magnitude)

        # Set y ticks
        if not self.config.y_ticks:
            for panel in self.panels:
                plot = self.panel_plots[panel]
                plot.hide_yticks()

        # Add legend?
        if legend_patches is not None:
            for panel in legend_patches:
                if len(legend_patches[panel]) == 0: continue
                plot = self.panel_plots[panel]
                #print(legend_patches[panel].values())
                #legend = plot.legend(legend_patches[panel].values()) #for_legend_parameters, **legend_properties)
                legend = plot.legend(legend_patches[panel])
                #legends.append(legend)

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
            if fs.is_directory(self.out_path): path = fs.join(self.out_path, "distributions" + self.config.format)
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
        :param distribution:
        :param label:
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
