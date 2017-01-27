#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.plot Contains the Plot class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
from io import BytesIO
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as patches
from collections import OrderedDict
from textwrap import wrap
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
from matplotlib import rc
from scipy.interpolate import interp1d

# Import the relevant PTS classes and modules
from ..tools.logging import log

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd'] * 6

distinguishable_colormaps = ["spring", "winter", "copper", "cool", "PRGn", "coolwarm"]

other_colormaps = ["YlGn", "YlGnBu", "YlOrRd", "Purples"]

# http://matplotlib.org/examples/color/named_colors.html
color_hex = colors.cnames

# Add the single letter colors
for name, rgb in colors.ColorConverter.colors.items():
    hex_ = colors.rgb2hex(rgb)
    color_hex[name] = hex_

color_plt_identifiers = color_hex.keys()

pretty_colors = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class Plot(object):
        
    """
    This class ...
    """

    def __init__(self, size):

        """
        This function ...
        :param size:
        """

        # Setup the figure
        self.figure = plt.figure(figsize=size)
        self.ax = self.figure.gca()

        # The plot path
        self.path = None

        # Properties
        self.add_borders = False
        self.transparent = True
        self.format = None

    # -----------------------------------------------------------------

    def set_borders(self, config):

        """
        This function ...
        :return:
        """

        # Set border width
        if config.add_borders: [i.set_linewidth(config.borderwidth) for i in self.ax.spines.itervalues()]

        # Remove borders
        else:
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['right'].set_visible(False)
            self.ax.spines['bottom'].set_visible(False)
            self.ax.spines['left'].set_visible(False)
            self.ax.tick_params(axis=u'both', which=u'both', length=0)

    # -----------------------------------------------------------------

    def add_legend(self, ax, config, legend_title=None):

        """
        This function ...
        :param ax:
        :param config:
        :param legend_title:
        :return:
        """

        # if nlegends > 1: percentage = 25.
        # else: percentage = 10.
        percentage = 20.

        # Shrink current axis's height by a certain percentage on the bottom
        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0 + box.height * percentage / 100., box.width, box.height * (100 - percentage) / 100.])

        # Plot legend
        legend_title = r"\underline{" + legend_title + "}"
        legend = plt.legend(loc="upper center", title=legend_title, bbox_to_anchor=(0.5, -0.25), fancybox=False,
                            shadow=False, ncol=4)

        frame = legend.get_frame()

        # Set legend frame color and line width
        if config.add_legend_border:
            frame.set_linewidth(config.legend_borderwidth)
        else:
            frame.set_linewidth(0)
        frame.set_edgecolor(config.legend_bordercolor)

        # Set background color
        frame.set_facecolor('0.85')
        legend.legendPatch.set_alpha(0.75)

        index = 0

        # Move to foreground
        legend.set_zorder(100 + index)

        # Set fontsize
        plt.setp(legend.get_title(), fontsize=str(config.legend_title_fontsize))

    # -----------------------------------------------------------------

    def set_ticks(self, config, x_range, nxticks):

        """
        This function ...
        :param config:
        :param x_range:
        :param nxticks:
        :return:
        """

        # Format the axis ticks and create a grid
        ticks = x_range.log(nxticks, fancy=True)
        self.ax.set_xticks(ticks)
        self.ax.set_xticklabels(ticks)

        # Tick formatter
        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))

        # Set ticks fontsize
        plt.setp(self.ax.get_xticklabels(), rotation='horizontal', fontsize=config.ticks_fontsize)
        plt.setp(self.ax.get_yticklabels(), rotation='horizontal', fontsize=config.ticks_fontsize)

    # -----------------------------------------------------------------

    def set_grid(self, config, which="major"):

        """
        This function ...
        :param config:
        :param which:
        :return:
        """

        if config.add_grid: plt.grid(linewidth=config.grid_linewidth, linestyle=config.grid_linestyle, color=config.grid_color, which=which)

    # -----------------------------------------------------------------

    def set_background(self, config):

        """
        This function ...
        :param config:
        :return:
        """

        if config.transparent_background:

            # Set transparent background
            for item in [self.figure, self.ax]:
                item.patch.set_visible(False)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the plot ...")

        # Show
        plt.show()

    # -----------------------------------------------------------------

    def saveto(self, out):

        """
        This function ...
        :param out:
        :return:
        """

        # Inform the user
        log.info("Saving the plot ...")

        # Differentiate
        if isinstance(out, BytesIO): self.to_buffer(out)
        else: self.to_file(out)

    # -----------------------------------------------------------------

    def to_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the plot to " + str(path) + " ...")

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent, format=self.format)

    # -----------------------------------------------------------------

    def to_buffer(self, buf):

        """
        This function ...
        :param buf:
        :return:
        """

        # Inform the user
        log.info("Saving the SED plot to a buffer ...")

        # Save to buffer
        plt.savefig(buf, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent, format=self.format)

    # -----------------------------------------------------------------

    def close(self):

        """
        This function ...
        :return:
        """

        plt.close()

# -----------------------------------------------------------------

def get_plot_wavelength_limits(min_wavelength, max_wavelength):

    """
    This function ...
    :param min_wavelength:
    :param max_wavelength:
    :return:
    """

    log_min_wavelength = np.log10(min_wavelength)
    log_max_wavelength = np.log10(max_wavelength)

    plot_min_log_wavelength = math.floor(log_min_wavelength)
    plot_max_log_wavelength = math.ceil(log_max_wavelength)

    plot_min_wavelength = 10.**plot_min_log_wavelength
    plot_max_wavelength = 10.**plot_max_log_wavelength

    return plot_min_wavelength, plot_max_wavelength

# -----------------------------------------------------------------
