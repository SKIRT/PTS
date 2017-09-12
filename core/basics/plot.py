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
    This fucntion ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, **kwargs):

        """
        Thisf unction ...
        :param kwargs:
        """

        # Set the plot
        self._plot = kwargs.pop("plot", None)

    # -----------------------------------------------------------------

    @property
    def initialized(self):

        """
        This function ...
        :return:
        """

        return self._plot is not None

# -----------------------------------------------------------------

class Figure(object):

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
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(BokehPlot, self).__init__(**kwargs)

        # Initialize the plot
        if not self.initialized: self.initialize(**kwargs)

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        This function ...
        :return:
        """

        from bokeh.plotting import figure

        # The Bokeh plot
        #self._plot = figure(plot_width=250, plot_height=250, title=None)

        if "share_x" in kwargs and "share_y" in kwargs: self._plot = figure(x_range=kwargs.pop("share_x")._plot.x_range, y_range=kwargs.pop("share_y")._plot.y_range)
        elif "share_x" in kwargs: self._plot = figure(x_range=kwargs.pop("share_x")._plot.x_range)
        elif "share_y" in kwargs: self._plot = figure(y_range=kwargs.pop("share_y")._plot.y_range)
        else: self._plot = figure()

    # -----------------------------------------------------------------

    def plot(self, x, y, color="red", point_kwargs=None, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param color:
        :param point_kwargs:
        :param kwargs:
        :return:
        """

        if point_kwargs is None: point_kwargs = {}

        # Plot data points
        self._plot.circle(x, y, color=color, **point_kwargs)

    # -----------------------------------------------------------------

    def errorbar(self, x, y, xerr=None, yerr=None, color='red', point_kwargs=None, error_kwargs=None, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param xerr:
        :param yerr:
        :param color:
        :param point_kwargs:
        :param error_kwargs:
        :param kwargs:
        :return:
        """

        if point_kwargs is None: point_kwargs = {}
        if error_kwargs is None: error_kwargs = {}

        # Plot data points
        self._plot.circle(x, y, color=color, **point_kwargs)

        from ..tools import types

        if types.is_real_type(x): x = [x]
        if types.is_real_type(y): y = [y]
        if types.is_real_type(xerr): xerr = [xerr]
        if types.is_real_type(yerr): yerr = [yerr]

        # X error bars
        if xerr is not None:

            x_err_x = []
            x_err_y = []

            for px, py, err in zip(x, y, xerr):
                x_err_x.append((px - err, px + err))
                x_err_y.append((py, py))

            self._plot.multi_line(x_err_x, x_err_y, color=color, **error_kwargs)

        # Y error bars
        if yerr is not None:

            y_err_x = []
            y_err_y = []

            for px, py, err in zip(x, y, yerr):
                y_err_x.append((px, px))
                y_err_y.append((py - err, py + err))

            self._plot.multi_line(y_err_x, y_err_y, color=color, **error_kwargs)

    # -----------------------------------------------------------------

    def legend(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def add_artist(self, *args, **kwargs):

        """
        Thisn function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def axhline(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def axvline(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_xlim(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_ylim(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_xticks(self, ticks, fontsize=None):

        """
        This function ...
        :param ticks:
        :param fontsize:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_yticks(self, ticks=None, fontsize=None):

        """
        This function ...
        :param ticks:
        :param fontsize:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_xscale(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_yscale(self, *args, **kwargs):

        """
        Thisf unction ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_xlabel(self, *args, **kwargs):

        """
        Thisf ucntion ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_ylabel(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        pass

# -----------------------------------------------------------------

class BokehFigure(Figure):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(BokehFigure, self).__init__(**kwargs)

        # Rows and columns
        self.rows = []
        self.columns = []

    # -----------------------------------------------------------------

    def create_column(self, size, share_axis=False, height_ratios=None):

        """
        This function ...
        :param size:
        :param share_axis:
        :param height_ratios:
        :return:
        """

        plots = []

        if share_axis:

            first_plot = BokehPlot()
            plots.append(first_plot)
            for i in range(1, size):
                next_plot = BokehPlot(share_x=first_plot)
                plots.append(next_plot)

        else:

            for i in range(size):
                plot = BokehPlot()
                plots.append(plot)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_row(self, size, share_axis=False, width_ratios=None):

        """
        This function ...
        :param size:
        :param share_axis:
        :param width_ratios:
        :return:
        """

        plots = []

        if share_axis:

            first_plot = BokehPlot()
            plots.append(first_plot)
            for i in range(1, size):
                next_plot = BokehPlot(share_y=first_plot)
                plots.append(next_plot)

        else:

            for i in range(size):
                plot = BokehPlot()
                plots.append(plot)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    @property
    def nrows(self):

        """
        This function ...
        :return:
        """

        return len(self.rows)

    # -----------------------------------------------------------------

    @property
    def has_rows(self):

        """
        Thisf unction ...
        :return:
        """

        return self.nrows > 0

    # -----------------------------------------------------------------

    @property
    def row_size(self):

        """
        This function ...
        :return:
        """

        return len(self.rows[0])

    # -----------------------------------------------------------------

    def add_row(self, *figures):

        """
        Tihs function ...
        :return:
        """

        # Check
        if self.has_columns: raise ValueError("Cannot add rows in column appending mode")

        # Check size
        if self.has_rows and len(figures) != self.row_size: raise ValueError("Invalid number of values: must be " + str(self.row_size))

    # -----------------------------------------------------------------

    @property
    def ncolumns(self):

        """
        This function ...
        :return:
        """

        return len(self.columns)

    # -----------------------------------------------------------------

    @property
    def has_columns(self):

        """
        This function ...
        :return:
        """

        return self.ncolumns > 0

    # -----------------------------------------------------------------

    @property
    def column_size(self):

        """
        This function ...
        :return:
        """

        return len(self.columns[0])

    # -----------------------------------------------------------------

    def add_column(self, *figures):

        """
        This function ...
        :param figures:
        :return:
        """

        # Check
        if self.has_rows: raise ValueError("Cannot add rows in row appending mode")

        # Check size
        if self.has_columns and len(figures) != self.column_size: raise ValueError("Invalid number of values: must be " + str(self.column_size))

    # -----------------------------------------------------------------

    @property
    def columns_to_rows(self):

        """
        This function ...
        :return:
        """

        rows = []
        for row_index in range(self.column_size):

            # Take the row_index'th value of each column
            row = []
            for column_index in range(self.ncolumns):
                row.append(self.columns[column_index][row_index])

            # Add the row
            rows.append(row)

        # Return the rows
        return rows

    # -----------------------------------------------------------------

    @property
    def grid(self):

        """
        This function ...
        :return:
        """

        from bokeh.layouts import gridplot

        if self.has_rows: return gridplot(self.rows)
        elif self.has_columns: return gridplot(self.columns_to_rows)
        else: raise ValueError("No rows or columns")

    # -----------------------------------------------------------------

    def to_html(self):

        """
        Thisf unction ...
        :return:
        """

        from bokeh.resources import CDN
        from bokeh.embed import file_html

        # Generate the html
        html = file_html(self.grid, CDN, self.title)

        # Return the HTML
        return html

    # -----------------------------------------------------------------

    def to_html_components(self):

        """
        This function ...
        :return:
        """

        from bokeh.embed import components

        script, div = components(self.grid)
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
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MPLPlot, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    def plot(self, x, y, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        :return:
        """

        self._plot.plot(x, y, **kwargs)

    # -----------------------------------------------------------------

    def errorbar(self, x, y, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        :return:
        """

        # Plot
        self._plot.errorbar(x, y, **kwargs)

    # -----------------------------------------------------------------

    def legend(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return self._plot.legend(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_artist(self, *args, **kwargs):

        """
        Thisn function ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.add_artist(*args, **kwargs)

    # -----------------------------------------------------------------

    def axhline(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Set linestyle and limit for axis2
        self._plot.axhline(*args, **kwargs)

    # -----------------------------------------------------------------

    def axvline(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.axvline(*args, **kwargs)

    # -----------------------------------------------------------------

    def set_xlim(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.set_xlim(*args, **kwargs)

    # -----------------------------------------------------------------

    def set_ylim(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.set_ylim(*args, **kwargs)

    # -----------------------------------------------------------------

    def set_xticks(self, ticks, fontsize=None):

        """
        This function ...
        :param ticks:
        :param fontsize:
        :return:
        """

        if ticks is not None:
            self._plot.set_xticks(ticks)
            self._plot.set_xticklabels(ticks)

        # Format
        self._plot.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        plt.setp(self._plot.get_xticklabels(), rotation='horizontal', fontsize=fontsize)

    # -----------------------------------------------------------------

    def set_yticks(self, ticks=None, fontsize=None):

        """
        This function ...
        :param ticks:
        :param fontsize:
        :return:
        """

        if ticks is not None:
            self._plot.set_yticks(ticks)
            self._plot.set_yticklabels(ticks)

        # Format
        self._plot.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        plt.setp(self._plot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)

    # -----------------------------------------------------------------

    def set_xscale(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Set axis label
        self._plot.set_xscale(*args, **kwargs)

    # -----------------------------------------------------------------

    def set_yscale(self, *args, **kwargs):

        """
        Thisf unction ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.set_yscale(*args, **kwargs)

    # -----------------------------------------------------------------

    def set_xlabel(self, *args, **kwargs):

        """
        Thisf ucntion ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.set_xlabel(*args, **kwargs)

    # -----------------------------------------------------------------

    def set_ylabel(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.set_ylabel(*args, **kwargs)

# -----------------------------------------------------------------

class MPLFigure(Figure):
        
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
        super(MPLFigure, self).__init__(**kwargs)

        # Setup the figure
        self.figure = plt.figure(figsize=size)
        plt.clf()
        self.ax = self.figure.gca()

        # Properties
        self.add_borders = False
        self.transparent = False
        self.format = None

    # -----------------------------------------------------------------

    def create_column(self, size, share_axis=False, height_ratios=None):

        """
        This function ...
        :param size:
        :param share_axis:
        :param height_ratios:
        :return:
        """

        # Make the grid
        gs = gridspec.GridSpec(size, 1, height_ratios=height_ratios)

        # Create the (sub)plots
        plots = []
        if share_axis:

            first_plot = plt.subplot(gs[0])
            plots.append(MPLPlot(plot=first_plot))
            for index in range(1, size):
                next_plot = plt.subplot(gs[index], sharex=first_plot)
                plots.append(MPLPlot(plot=next_plot))
        else:

            for index in range(size):
                plot = plt.subplot(gs[index])
                plots.append(MPLPlot(plot=plot))

        # Return the plots
        return plots

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
