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
from abc import ABCMeta, abstractmethod
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
from ..basics.log import log

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

#pretty_colors = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# More:
pretty_colors = ["dodgerblue", "r", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen",
                 "lightcoral", "crimson", "saddlebrown", "mediumslateblue", "lightslategrey", "firebrick",
                 "orange", "darkcyan", "hotpink", "indianred", "aqua"]

# -----------------------------------------------------------------

mpl = "matplotlib"
bokeh = "bokeh"

# -----------------------------------------------------------------

plotting_libraries = [mpl, bokeh]

# -----------------------------------------------------------------

class Plot(object):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # The plot path
        self.path = kwargs.pop("path", None)

        # The title
        self.title = kwargs.pop("title", None)

    # -----------------------------------------------------------------

    @abstractmethod
    def show(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check if path is defined
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

# -----------------------------------------------------------------

# stylesheet and script can also be accessed with https!
bokeh_page_template = """<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <title>Bokeh Scatter Plots</title>

        <link rel="stylesheet" href="http://cdn.pydata.org/bokeh/release/bokeh-0.12.6.min.css" type="text/css" />
        <script type="text/javascript" src="http://cdn.pydata.org/bokeh/release/bokeh-0.12.6.min.js"></script>

        <!-- COPY/PASTE SCRIPT HERE -->

    </head>
    <body>
        <!-- INSERT DIVS HERE -->
    </body>
</html>
"""

# -----------------------------------------------------------------

class BokehPlot(Plot):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(BokehPlot, self).__init__(**kwargs)

        from bokeh.plotting import figure

        #
        self.figure = figure()
        self.figure.circle([1, 2], [3, 4])

    # -----------------------------------------------------------------

    def to_html(self):

        """
        Thisf unction ...
        :return:
        """

        from bokeh.resources import CDN
        from bokeh.embed import file_html

        # Generate the html
        html = file_html(self.figure, CDN, self.title)

        # Return the HTML
        return html

    # -----------------------------------------------------------------

    def to_html_components(self):

        """
        This function ...
        :return:
        """

        from bokeh.embed import components

        script, div = components(self.figure)
        return script, div

    # -----------------------------------------------------------------

    def write_html(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        from ..tools import filesystem as fs
        html = self.to_html()
        fs.write_text(path, html)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        from ..tools import filesystem as fs
        from ..tools import introspection
        from ..tools import time
        from ..tools import browser

        # Determine temporary filepath
        filename = time.unique_name("bokeh_plot") + ".html"
        filepath = fs.join(introspection.pts_temp_dir, filename)

        # Write the html
        self.write_html(filepath)

        # Open the file in a web browser
        browser.open_path(filepath)

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        pass

# -----------------------------------------------------------------

class MPLPlot(Plot):
        
    """
    This class ...
    """

    def __init__(self, size=(10,6), **kwargs):

        """
        This function ...
        :param size:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MPLPlot, self).__init__(**kwargs)

        # Setup the figure
        self.figure = plt.figure(figsize=size)
        plt.clf()
        self.ax = self.figure.gca()

        # Properties
        self.add_borders = False
        self.transparent = False
        self.format = None

    # -----------------------------------------------------------------

    def set_x_log_scale(self):

        """
        This function ...
        :return:
        """

        plt.xscale("log")

    # -----------------------------------------------------------------

    def set_y_log_scale(self):

        """
        This function ...
        :return:
        """

        plt.yscale("log")

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

    def set_labels(self, x_label, y_label, fontsize=18):

        """
        This function ...
        :param x_label:
        :param y_label:
        :param fontsize:
        :return:
        """

        # Labels
        plt.xlabel(x_label, fontsize=fontsize)
        plt.ylabel(y_label, fontsize=fontsize)

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

    def add_curve(self, curve, label=None):

        """
        This function ...
        :param curve:
        :param label:
        :return:
        """

        self.ax.plot(curve.x_data, curve.y_data, label=label)

    # -----------------------------------------------------------------

    def finish(self, out=None):

        """
        This function ...
        :param out:
        :return:
        """

        # Set tight layout
        plt.tight_layout()

        # Save the plot
        if out is None: self.show()
        else: self.saveto(out)

        # Close
        self.close()

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
