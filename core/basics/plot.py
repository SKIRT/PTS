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
from textwrap import wrap
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import LinearLocator, LogLocator, AutoMinorLocator, AutoLocator, NullLocator
from matplotlib.ticker import ScalarFormatter, NullFormatter, LogFormatter, PercentFormatter, EngFormatter, LogFormatterMathtext, LogFormatterSciNotation
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import cbook
from matplotlib.legend import Legend
from matplotlib.colorbar import colorbar_factory
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LogNorm, Normalize
#from matplotlib.ticker import LogFormatterMathtext, LogLocator

# Import the relevant PTS classes and modules
from ..basics.log import log

# -----------------------------------------------------------------

#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8')

# -----------------------------------------------------------------

class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self, vmin, vmax):  # Override function that finds format to use.
        #print(vmin, vmax)
        self.format = "%1.1f"  # Give format here

# -----------------------------------------------------------------

uniform_colormaps = ["viridis", "plasma", "inferno", "magma"]
misc_colormaps = ["flag", "prism", "ocean", "gist_earth", "terrain", "gist_stern", "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "brg", "hsv", "gist_rainbow", "rainbow", "jet", "nipy_spectral", "gist_ncar"]
sequential_colormaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd',
                        'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg',
                        'gist_gray', 'gray', 'bone', 'pink', 'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                        'hot', 'afmhot', 'gist_heat', 'copper']
diverging_colormaps = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
qualitative_colormaps = ['Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c']

# -----------------------------------------------------------------

normal_colormaps = sequential_colormaps + uniform_colormaps + misc_colormaps
all_colormaps = uniform_colormaps + sequential_colormaps + misc_colormaps + diverging_colormaps + qualitative_colormaps

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

#pretty_colors_no_yellows = [x for x in pretty_colors if x != 'yellow']
dark_pretty_colors = ["dodgerblue", "r", "purple", "darkorange", "darkblue", "teal", "darkgreen",
                 "lightcoral", "crimson", "saddlebrown", "mediumslateblue", "lightslategrey", "firebrick",
                 "orange", "darkcyan", "hotpink", "indianred", "aqua"]

# -----------------------------------------------------------------

mpl = "matplotlib"
bokeh = "bokeh"

# -----------------------------------------------------------------

plotting_libraries = [mpl, bokeh]

# -----------------------------------------------------------------

pdf = "pdf"
png = "png"
svg = "svg"
html = "html"

# -----------------------------------------------------------------

plotting_formats = [pdf, png, svg, html]

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

# BOKEH MARKERS:
# Asterisk
# Circle
# CircleCross
# CircleX
# Cross
# Diamond
# DiamondCross
# InvertedTriangle
# Square
# SquareCross
# SquareX
# Triangle
# X

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

        options = dict()
        if "share_x" in kwargs: options["x_range"] = kwargs.pop("share_x")._plot.x_range
        if "share_y" in kwargs: options["y_range"] = kwargs.pop("share_y")._plot.y_range
        if "x_axis_location" in kwargs: options["x_axis_location"] = kwargs.pop("x_axis_location")
        if "y_axis_location" in kwargs: options["y_axis_location"] = kwargs.pop("y_axis_location")
        if "width" in kwargs: options["plot_width"] = kwargs.pop("width")
        if "height" in kwargs: options["plot_height"] = kwargs.pop("height")

        if "x_range" in kwargs: options["x_range"] = kwargs.pop("x_range")
        if "y_range" in kwargs: options["y_range"] = kwargs.pop("y_range")

        # Number of minor ticks between adjacent x-axis major ticks
        # Number of minor ticks between adjacent y-axis major ticks
        if "x_minor_ticks" in kwargs: options["x_minor_ticks"] = kwargs.pop("x_minor_ticks")
        if "y_minor_ticks" in kwargs: options["y_minor_ticks"] = kwargs.pop("y_minor_ticks")

        if "x_axis_label" in kwargs: options["x_axis_label"] = kwargs.pop("x_axis_label")
        if "y_axis_label" in kwargs: options["y_axis_label"] = kwargs.pop("y_axis_label")

        # "linear", "log", "datetime"
        if "x_axis_type" in kwargs: options["x_axis_type"] = kwargs.pop("x_axis_type")
        if "y_axis_type" in kwargs: options["y_axis_type"] = kwargs.pop("y_axis_type")

        default_tools = ["pan", "wheel_zoom", "box_zoom", "save", "reset", "help"]
        options["tools"] = kwargs.pop("tools", default_tools)
        #if "tools" in kwargs:

        # Specify a drag tool to be active when the plot is displayed.
        # box_select, box_zoom, lasso_select, 'pan', 'xpan', 'ypan',
        if "box_zoom" in options["tools"]: options["active_drag"] = "box_zoom"

        # Specify an inspection tool or sequence of inspection tools to be active when the plot is displayed.
        # 'crosshair', 'hover',
        if "hover" in options["tools"]: options["active_inspect"] = "hover"

        # Specify a scroll/pinch tool to be active when the plot is displayed.
        # 'wheel_zoom', 'xwheel_zoom', 'ywheel_zoom', 'xwheel_pan', 'ywheel_pan'
        if "wheel_zoom" in options["tools"]: options["active_scroll"] = "wheel_zoom"

        # Specify a tap/click tool to be active when the plot is displayed.
        # 'poly_select', 'tap',
        if "tap" in options["tools"]: options["active_tap"] = "tap"

        # above, below, left, right
        if "toolbar_location" in kwargs: options["toolbar_location"] = kwargs.pop("toolbar_location")
        if "toolbar_sticky" in kwargs: options["toolbar_sticky"] = kwargs.pop("toolbar_sticky") # default = True

        # Create the plot figure
        self._plot = figure(**options)

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

    def scatter(self, x, y, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        :return:
        """

        raise NotImplementedError("Not implemented for Bokeh")

    # -----------------------------------------------------------------

    def text(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        raise NotImplementedError("Not yet implemented")

    # -----------------------------------------------------------------

    def horizontal_arrow(self, min_x, max_x, y):

        """
        This function ...
        :param min_x:
        :param max_x:
        :param y:
        :return:
        """

        raise NotImplementedError("Not yet implemented")

    # -----------------------------------------------------------------

    def vertical_arrow(self, x, min_y, max_y):

        """
        This function ...
        :param x:
        :param min_y:
        :param max_y:
        :return:
        """

        raise NotImplementedError("Not yet implemented")

    def vlines(self, x, ymin, ymax, **kwargs):

        """
        This function ...
        :param x:
        :param ymin:
        :param ymax:
        :param kwargs:
        :return:
        """

        raise NotImplementedError("Not yet implemented")

    # -----------------------------------------------------------------

    def fill(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        raise NotImplementedError("Not yet implemented")

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
            x_err_color = []

            for px, py, err in zip(x, y, xerr):
                x_err_x.append((px - err, px + err))
                x_err_y.append((py, py))
                x_err_color.append((color, color))

            self._plot.multi_line(xs=x_err_x, ys=x_err_y, color=x_err_color, **error_kwargs)

        # Y error bars
        if yerr is not None:

            y_err_x = []
            y_err_y = []
            y_err_color = []

            for px, py, err in zip(x, y, yerr):
                y_err_x.append((px, px))
                y_err_y.append((py - err, py + err))
                y_err_color.append((color, color))

            self._plot.multi_line(xs=y_err_x, ys=y_err_y, color=y_err_color, **error_kwargs)

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

    # def create_legends(self, nlegends):
    #
    #     """
    #     This function ...
    #     :param nlegends:
    #     :return:
    #     """
    #
    #     return None

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

        from bokeh.plotting.helpers import _get_scale

        scale = args[0]
        self._plot.x_scale = _get_scale(self._plot.x_range, scale)

    # -----------------------------------------------------------------

    def set_yscale(self, *args, **kwargs):

        """
        Thisf unction ...
        :param args:
        :param kwargs:
        :return:
        """

        from bokeh.plotting.helpers import _get_scale

        scale = args[0]
        self._plot.y_scale = _get_scale(self._plot.y_range, scale)

    # -----------------------------------------------------------------

    def set_xlabel(self, *args, **kwargs):

        """
        Thisf ucntion ...
        :param args:
        :param kwargs:
        :return:
        """

        label = args[0]
        self._plot.xaxis.axis_label = label

    # -----------------------------------------------------------------

    def set_ylabel(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        label = args[0]
        self._plot.yaxis.axis_label = label

    # -----------------------------------------------------------------

    def hide_axes(self):

        """
        This function ...
        :return:
        """

        self._plot.axis.visible = False

    # -----------------------------------------------------------------

    def hide_xaxis(self):

        """
        This function ...
        :return:
        """

        self._plot.xaxis.visible = False

    # -----------------------------------------------------------------

    def hide_yaxis(self):

        """
        This function ...
        :return:
        """

        self._plot.yaxis.visible = False

    # -----------------------------------------------------------------

    def hide_grid(self):

        """
        This function ...
        """

        self.hide_xgrid()
        self.hide_ygrid()

    # -----------------------------------------------------------------

    def hide_xgrid(self):

        """
        This function ...
        :return:
        """

        self._plot.xgrid.grid_line_color = None

    # -----------------------------------------------------------------

    def hide_ygrid(self):

        """
        This function ...
        :return:
        """

        self._plot.ygrid.grid_line_color = None

    # -----------------------------------------------------------------

    def set_xgrid(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Examples:
        self._plot.xgrid.grid_line_alpha = 0.5
        self._plot.xgrid.grid_line_dash = [6, 4]
        self._plot.xgrid.minor_grid_line_color = 'navy'
        self._plot.xgrid.bounds = (2, 4)

    # -----------------------------------------------------------------

    @property
    def xaxis(self):

        """
        This function ...
        :return:
        """

        return self._plot.xaxis

    # -----------------------------------------------------------------

    @property
    def yaxis(self):

        """
        This function ...
        :return:
        """

        return self._plot.yaxis

    # -----------------------------------------------------------------

    def set_xaxis_location(self, location):

        """
        Thisn function ...
        :param location:
        :return:
        """

        if location == "below": self._plot.below.append(self.xaxis)
        elif location == "above": self._plot.above.append(self.xaxis)
        else: raise ValueError("")

    # -----------------------------------------------------------------

    def set_yaxis_location(self, location):

        """
        This function ...
        :param location:
        :return:
        """

        if location == "left": self._plot.left.append(self.yaxis)
        elif location == "right": self._plot.right.append(self.yaxis)
        else: raise ValueError("")

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
                if i == size - 1: next_plot = BokehPlot(share_x=first_plot, x_axis_location="below")
                else: next_plot = BokehPlot(share_x=first_plot)
                #next_plot.hide_xaxis()
                #next_plot.set_xaxis_location("below")
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
                #next_plot.hide_yaxis()
                next_plot.set_yaxis_location("right")
                plots.append(next_plot)

        else:

            for i in range(size):
                plot = BokehPlot()
                plots.append(plot)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_one_plot(self):

        """
        Thisf unction ...
        :return:
        """

        plot = BokehPlot()
        return plot

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

    def add_row(self, *plots):

        """
        Tihs function ...
        :param plots:
        :return:
        """

        # Check
        if self.has_columns: raise ValueError("Cannot add rows in column appending mode")

        # Check size
        if self.has_rows and len(plots) != self.row_size: raise ValueError("Invalid number of values: must be " + str(self.row_size))

        # Add row
        bokeh_plots = [plot._plot for plot in plots]
        self.rows.append(bokeh_plots)

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

    def add_column(self, *plots):

        """
        This function ...
        :param plots:
        :return:
        """

        # Check
        if self.has_rows: raise ValueError("Cannot add rows in row appending mode")

        # Check size
        if self.has_columns and len(plots) != self.column_size: raise ValueError("Invalid number of values: must be " + str(self.column_size))

        # Add
        bokeh_plots = [plot._plot for plot in plots]
        self.columns.append(bokeh_plots)

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

    def set_grid(self, config, which="major"):

        """
        Thisf unction ...
        :param config:
        :param which:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_borders(self, *args):

        """
        This function ...
        :param args:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        pass

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

    def to_html_components(self, wrap_script=True):

        """
        This function ...
        :param wrap_script:
        :return:
        """

        from bokeh.embed import components

        script, div = components(self.grid, wrap_script=wrap_script)
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

        if path.endswith("html"): self.saveto_html(path)
        elif path.endswith("png"): self.saveto_png(path)
        elif path.endswith("svg"): self.saveto_svg(path)

    # -----------------------------------------------------------------

    def saveto_html(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.write_html(path)

    # -----------------------------------------------------------------

    def saveto_png(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # FOR transparent background:
        # Plot.background_fill_color and Plot.border_fill_color properties to None.

        from bokeh.io import export_png
        export_png(self.grid, filename=path)

    # -----------------------------------------------------------------

    def saveto_svg(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        from bokeh.io import export_svgs

        # FOR EACH PLOT SEPERATELY!
        #plot.output_backend = "svg"
        #export_svgs(plot, filename="plot.svg")

    # -----------------------------------------------------------------

    def close(self):

        """
        This function ...
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

    @property
    def bounding_box(self):
        return self.axes.get_position()

    # -----------------------------------------------------------------

    @property
    def bounds(self):
        return self.bounding_box.bounds

    # -----------------------------------------------------------------

    @property
    def x_min(self):
        return self.bounds[0]

    # -----------------------------------------------------------------

    @property
    def y_min(self):
        return self.bounds[1]

    # -----------------------------------------------------------------

    @property
    def width(self):
        return self.bounds[2]

    # -----------------------------------------------------------------

    @property
    def height(self):
        return self.bounds[3]

    # -----------------------------------------------------------------

    @property
    def x_max(self):
        return self.x_min + self.width

    # -----------------------------------------------------------------

    @property
    def y_max(self):
        return self.y_min + self.height

    # -----------------------------------------------------------------

    def add_patch(self, patch):
        self.axes.add_patch(patch)

    # -----------------------------------------------------------------

    def set_background_color(self, color):
        self.axes.set_facecolor(color)

    # -----------------------------------------------------------------

    def plot(self, x, y, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        :return:
        """

        # Plot
        return self._plot.plot(x, y, **kwargs)

    # -----------------------------------------------------------------

    def contour(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Plot
        return self._plot.contour(*args, **kwargs)

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
        return self._plot.errorbar(x, y, **kwargs)

    # -----------------------------------------------------------------

    def scatter(self, x, y, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        :return:
        """

        return self._plot.scatter(x, y, **kwargs)

    # -----------------------------------------------------------------

    def imshow(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return self._plot.imshow(*args, **kwargs)

    # -----------------------------------------------------------------

    def text(self, x, y, s, *args, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param s:
        :param args:
        :param kwargs:
        :return:
        """

        return self._plot.text(x, y, s, *args, **kwargs)

    # -----------------------------------------------------------------

    def add_text(self, text, vertical_position="top", horizontal_position="center", x_shift=0, y_shift=0, fontsize=None):

        """
        This function ...
        :param text:
        :param vertical_position:
        :param horizontal_position:
        :param x_shift:
        :param y_shift:
        :param fontsize:
        :return:
        """

        # Set x position
        if horizontal_position == "left":
            x = 0.05 + x_shift
            halignment = "left"
        elif horizontal_position == "right":
            x = 0.95 + x_shift
            halignment = "right"
        elif horizontal_position == "center":
            x = 0.5 + x_shift
            halignment = "center"
        else: raise ValueError("Invalid horizontal position: '" + horizontal_position + "'")

        # Set y position
        if vertical_position == "top":
            #y = self.y_min + 0.9 * self.height
            y = 0.95 + y_shift
            valignment = "top"
        elif vertical_position == "bottom":
            #y = self.y_min + 0.1 * self.height
            y = 0.05 + y_shift
            valignment = "bottom"
        elif vertical_position == "center":
            y = 0.5 + y_shift
            valignment = "center"
        else: raise ValueError("Invalid vertical position: '" + vertical_position + "'")

        # Create
        if fontsize is not None: return self.text(x, y, text, horizontalalignment=halignment, verticalalignment=valignment, transform=self.axes.transAxes, fontsize=fontsize)
        else: return self.text(x, y, text, horizontalalignment=halignment, verticalalignment=valignment, transform=self.axes.transAxes)

    # -----------------------------------------------------------------

    def horizontal_arrow(self, x_min, x_max, y, **kwargs):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y:
        :param kwargs:
        :return:
        """

        return self._plot.annotate(s='', xy=(x_max, y), xytext=(x_min, y), arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0), **kwargs)

    # -----------------------------------------------------------------

    def vertical_arrow(self, x, y_min, y_max, **kwargs):

        """
        This function ...
        :param x:
        :param y_min:
        :param y_max:
        :param kwargs:
        :return:
        """

        return self._plot.annotate(s='', xy=(x, y_max), xytext=(x, y_min), arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0), **kwargs)

    # -----------------------------------------------------------------

    def vlines(self, x, ymin, ymax, **kwargs):

        """
        This function ...
        :param x:
        :param ymin:
        :param ymax:
        :param kwargs:
        :return:
        """

        # Plot
        return self._plot.vlines(x, ymin, ymax, **kwargs)

    # -----------------------------------------------------------------

    def bar(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return self._plot.bar(*args, **kwargs)

    # -----------------------------------------------------------------

    def fill(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return self._plot.fill(*args, **kwargs)

    # -----------------------------------------------------------------

    def create_legend(self, handles, labels, **kwargs):

        """
        Thisf unction ...
        :param handles:
        :param labels:
        :param kwargs:
        :return:
        """

        return Legend(self.axes, handles, labels, **kwargs)

    # -----------------------------------------------------------------

    @property
    def legend_handles(self):
        # changed with matplotlib version?
        #return [handle for handle in self.axes._get_legend_handles() if true_and_not_startswith(handle.get_label(), "_")]
        return self.axes.get_legend_handles_labels()[0]

    # -----------------------------------------------------------------

    @property
    def legend_labels(self):
        # changed with matplotlib version?
        #return [handle.get_label() for handle in self.legend_handles]
        return self.axes.get_legend_handles_labels()[1]

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
        return self._plot.axhline(*args, **kwargs)

    # -----------------------------------------------------------------

    def axvline(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        return self._plot.axvline(*args, **kwargs)

    # -----------------------------------------------------------------

    def disable_auto_scaling(self):

        """
        This function ...
        :return:
        """

        self.axes.set_autoscale_on(False)

    # -----------------------------------------------------------------

    def enable_auto_scaling(self):

        """
        This function ...
        :return:
        """

        self.axes.set_autoscale_on(True)

    # -----------------------------------------------------------------

    def disable_x_auto_scaling(self):

        """
        This function ...
        :return:
        """

        self.axes.set_autoscalex_on(False)

    # -----------------------------------------------------------------

    def enable_x_auto_scaling(self):

        """
        This function ...
        :return:
        """

        self.axes.set_autoscalex_on(True)

    # -----------------------------------------------------------------

    def disable_y_auto_scaling(self):

        """
        This function ...
        :return:
        """

        self.axes.set_autoscaley_on(False)

    # -----------------------------------------------------------------

    def enable_y_auto_scaling(self):

        """
        This function ...
        :return:
        """

        self.axes.set_autoscaley_on(True)

    # -----------------------------------------------------------------

    def set_xlim(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.set_xlim(*args, **kwargs)
        self.disable_x_auto_scaling()

    # -----------------------------------------------------------------

    def set_ylim(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        self._plot.set_ylim(*args, **kwargs)
        self.disable_y_auto_scaling()

    # -----------------------------------------------------------------

    @property
    def axes(self):
        return self._plot.axes

    # -----------------------------------------------------------------

    @property
    def xaxis(self):
        return self.axes.get_xaxis()

    # -----------------------------------------------------------------

    @property
    def yaxis(self):
        return self.axes.get_yaxis()

    # -----------------------------------------------------------------

    @property
    def xticklabels(self):
        return self._plot.get_xticklabels()

    # -----------------------------------------------------------------

    def set_xticks(self, ticks=None, tick_labels=None, fontsize=None, major_locator=None, minor_locator=None, major_formatter=None,
                   minor_formatter=None, log_subs=(1., 2., 5.), formatter=None, minor=True, magnitude=None):

        """
        This function ...
        :param ticks:
        :param tick_labels:
        :param fontsize:
        :param major_locator:
        :param minor_locator:
        :param major_formatter:
        :param minor_formatter:
        :param log_subs:
        :param formatter:
        :param minor:
        :param magnitude:
        :return:
        """

        # Set major formatter automatically
        if major_formatter is None:
            if formatter is None: major_formatter = FormatStrFormatter('%g') # default
            elif self.log_xscale:
                if formatter == "scalar": major_formatter = ScalarFormatter()
                elif formatter == "log": major_formatter = LogFormatter()
                elif formatter == "math": major_formatter = LogFormatterMathtext()
                elif formatter == "sci": major_formatter = LogFormatterSciNotation()
                elif formatter == "eng": major_formatter = EngFormatter()
                else: raise ValueError("Invalid formatter for log scale: '" + formatter + "'")
            elif self.linear_xscale:
                if formatter == "scalar": major_formatter = ScalarFormatter()
                elif formatter == "eng": major_formatter = EngFormatter()
                elif formatter == "percent": major_formatter = PercentFormatter()
                else: raise ValueError("Invalid formatter for linear scale : '" + formatter + "'")
            else: raise ValueError("Unknown xscale")

        # Set minor formatter automatically
        if minor_formatter is None: minor_formatter = NullFormatter()

        # Set major formatter
        self.xaxis.set_major_formatter(major_formatter)

        # Set minor formatter
        self.xaxis.set_minor_formatter(minor_formatter)

        # Set the ticks
        if ticks is not None:

            # Check
            if major_locator is not None: raise ValueError("Cannot specify ticks and pass major locator")
            if minor_locator is not None: raise ValueError("Cannot specify ticks and pass minor locator")

            # Set
            if tick_labels is None: tick_labels = ticks
            # print(tick_labels)
            self._plot.set_xticks(ticks)
            self._plot.set_xticklabels(tick_labels)

        # Ticks are not passed
        else:

            # Set locator automatically
            if major_locator is None:

                if self.log_xscale: major_locator = LogLocator(subs=log_subs)
                elif self.linear_xscale: major_locator = LinearLocator()
                else: raise ValueError("Unknown xscale")

            if minor_locator is None:

                if minor:  # have minor ticks
                    if self.log_xscale: minor_locator = None
                    elif self.linear_xscale: minor_locator = AutoMinorLocator()
                    else: raise ValueError("Unknown xscale")
                else: minor_locator = NullLocator()

            # Set the locators
            self.xaxis.set_major_locator(major_locator)
            if minor_locator is not None: self.xaxis.set_minor_locator(minor_locator)

        # Accomodate for the magnitude
        if magnitude is not None:
            if tick_labels is not None: raise ValueError("Tick labels cannot be passed when magnitude is also defined")
            ticks = self.xaxis.get_majorticklocs() # returns numpy array
            tick_labels = ticks / 10**magnitude #[tick / 10 ** magnitude for tick in ticks]
            self._plot.set_xticklabels(tick_labels)

        # Set cosmetic properties of the tick labels
        if fontsize is not None: plt.setp(self.xticklabels, rotation='horizontal', fontsize=fontsize)
        else: plt.setp(self.xticklabels, rotation='horizontal')

    # -----------------------------------------------------------------

    def set_xtick_labels(self, labels):

        """
        This function ...
        :param labels:
        :return:
        """

        self.axes.set_xticklabels(labels)

    # -----------------------------------------------------------------

    @property
    def yticklabels(self):
        return self._plot.get_yticklabels()

    # -----------------------------------------------------------------

    def set_yticks(self, ticks=None, fontsize=None, major_locator=None, minor_locator=None, major_formatter=None,
                   minor_formatter=None, log_subs=(1., 2., 5.), formatter=None, minor=True):

        """
        This function ...
        :param ticks:
        :param fontsize:
        :param major_locator:
        :param minor_locator:
        :param major_formatter:
        :param minor_formatter:
        :param log_subs:
        :param formatter:
        :param minor:
        :return:
        """

        # Set the ticks
        if ticks is not None:

            # Check
            if major_locator is not None: raise ValueError("Cannot specify ticks and pass major locator")
            if minor_locator is not None: raise ValueError("Cannot specify ticks and pass minor locator")

            # Set
            self._plot.set_yticks(ticks)
            self._plot.set_yticklabels(ticks)

        # Set locator automatically
        if major_locator is None:
            if self.log_yscale: major_locator = LogLocator(subs=log_subs)
            elif self.linear_yscale: major_locator = AutoLocator()
            else: raise ValueError("Unknown yscale")
        if minor_locator is None:
            if minor: # have minor ticks
                if self.log_yscale: minor_locator = None # auto
                elif self.linear_yscale: minor_locator = AutoMinorLocator() # auto
                else: raise ValueError("Unknown yscale")
            else: minor_locator = NullLocator()

        # Set the locators
        self.yaxis.set_major_locator(major_locator)
        if minor_locator is not None: self.yaxis.set_minor_locator(minor_locator)

        # Set major formatter automatically
        if major_formatter is None:
            if formatter is None: major_formatter = FormatStrFormatter('%g')  # default
            elif self.log_yscale:
                if formatter == "scalar": major_formatter = ScalarFormatter()
                elif formatter == "log": major_formatter = LogFormatter()
                elif formatter == "math": major_formatter = LogFormatterMathtext()
                elif formatter == "sci": major_formatter = LogFormatterSciNotation()
                elif formatter == "eng": major_formatter = EngFormatter()
                else: raise ValueError("Invalid formatter for log scale: '" + formatter + "'")
            elif self.linear_yscale:
                if formatter == "scalar": major_formatter = ScalarFormatter()
                elif formatter == "eng": major_formatter = EngFormatter()
                elif formatter == "percent": major_formatter = PercentFormatter()
                else: raise ValueError("Invalid formatter for linear scale : '" + formatter + "'")
            else: raise ValueError("Unknown yscale")

        # Set minor formatter automatically
        if minor_formatter is None: minor_formatter = NullFormatter()

        # Set major formatter
        self.yaxis.set_major_formatter(major_formatter)

        # Set minor formatter
        self.yaxis.set_minor_formatter(minor_formatter)

        # Set fontsize
        if fontsize is not None: plt.setp(self.yticklabels, rotation='horizontal', fontsize=fontsize)
        else: plt.setp(self.yticklabels, rotation='horizontal')

    # -----------------------------------------------------------------

    def set_ytick_labels(self, labels):

        """
        This function ...
        :param labels:
        :return:
        """

        self.axes.set_yticklabels(labels)

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

    @property
    def xscale(self):
        return self.xaxis.get_scale()

    # -----------------------------------------------------------------

    @property
    def linear_xscale(self):
        return self.xscale == "linear"

    # -----------------------------------------------------------------

    @property
    def log_xscale(self):
        return self.xscale == "log"

    # -----------------------------------------------------------------

    def set_yscale(self, *args, **kwargs):
        self._plot.set_yscale(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def yscale(self):
        return self.yaxis.get_scale()

    # -----------------------------------------------------------------

    @property
    def linear_yscale(self):
        return self.yscale == "linear"

    # -----------------------------------------------------------------

    @property
    def log_yscale(self):
        return self.yscale == "log"

    # -----------------------------------------------------------------

    def _clean_label(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # label = label.decode("utf8").replace("_", "\_").replace(u'\xa0', u' ').replace("&", "\&")
        label = label.decode("utf8").replace(u'\xa0', u' ').replace("&", "\&")
        if "$" in label:
            if label.count("$") != 2: log.warning("Cannot handle label '" + label + "'")
            else:
                before, inside, after = label.split("$")
                before = before.replace("_", "\_")
                after = after.replace("_", "\_")
                label = before + "$" + inside + "$" + after
        else: label = label.replace("_", "\_")
        return label

    # -----------------------------------------------------------------

    def set_xlabel(self, label, **kwargs):

        """
        This fucntion ...
        :param label:
        :param kwargs:
        :return:
        """

        label = self._clean_label(label)
        self._plot.set_xlabel(label, **kwargs)

    # -----------------------------------------------------------------

    def set_ylabel(self, label, **kwargs):

        """
        This function ...
        :param label:
        :param kwargs:
        :return:
        """

        label = self._clean_label(label)
        self._plot.set_ylabel(label, **kwargs)

    # -----------------------------------------------------------------

    def set_x_scalar(self):

        """
        This function ...
        :return:
        """

        self.xaxis.set_major_formatter(ScalarFormatterForceFormat(useOffset=True, useMathText=True))

    # -----------------------------------------------------------------

    def set_y_scalar(self):

        """
        This function ...
        :return:
        """

        self.yaxis.set_major_formatter(ScalarFormatterForceFormat(useOffset=True, useMathText=True))

    # -----------------------------------------------------------------

    def hide_xticks(self):

        """
        This function ...
        :return:
        """

        # No ticks
        self.xaxis.set_ticks([])

    # -----------------------------------------------------------------

    def hide_yticks(self):

        """
        This function ...
        :return:
        """

        # No ticks
        self.yaxis.set_ticks([])

    # -----------------------------------------------------------------

    def hide_xtick_labels(self, shared_axis=False):

        """
        This function ...
        :param shared_axis:
        :return:
        """

        # Shared axis required other method, or the other axis also loses its labels
        if shared_axis:
            plt.setp(self.axes.get_xticklabels(), visible=False)

        # Not shared
        else:

            # Set labels to empty strings
            labels = [item.get_text() for item in self.axes.get_xticklabels()]
            empty_string_labels = [''] * len(labels)
            self.set_xtick_labels(empty_string_labels)

    # -----------------------------------------------------------------

    def hide_ytick_labels(self, shared_axis=False):

        """
        This function ...
        :param shared_axis:
        :return:
        """

        # Shared axis requires other method, or the other axis also loses its labels
        if shared_axis:
            plt.setp(self.axes.get_yticklabels(), visible=False)

        # Not shared
        else:

            # Set labels to empty strings
            labels = [item.get_text() for item in self.axes.get_yticklabels()]
            empty_string_labels = [''] * len(labels)
            self.set_ytick_labels(empty_string_labels)

    # -----------------------------------------------------------------

    def hide_first_xtick_label(self):

        """
        This function ...
        :return:
        """

        # Get labels
        #labels = self.xtick_labels

        # Hide first
        #labels[0] = ""
        #print("XLABELS", labels)
        #self.set_xtick_labels(labels)

        # NEW
        a = self.axes.get_xticks().tolist()
        a[0] = " " # SPACE AND NOT EMPTY!!! '' is default and will be replaced by tick value when drawing!!
        self.axes.set_xticklabels(a)

    # -----------------------------------------------------------------

    def hide_first_ytick_label(self):

        """
        This function ...
        :return:
        """

        # Get labels
        #labels = self.ytick_labels

        # Hide first
        #labels[0] = ""
        #self.set_ytick_labels(labels)

        # NEW
        a = self.axes.get_yticks().tolist()
        a[0] = " "  # SPACE AND NOT EMPTY!!! '' is default and will be replaced by tick value when drawing!!
        self.axes.set_yticklabels(a)

    # -----------------------------------------------------------------

    def hide_last_xtick_label(self):

        """
        This function ...
        :return:
        """

        # Get labels
        #labels = self.xtick_labels

        # Hide last
        #labels[-1] = ""
        #print("XLABELS", labels)
        #self.set_xtick_labels(labels)

        # NEW
        a = self.axes.get_xticks().tolist()
        a[-1] = " "  # SPACE AND NOT EMPTY!!! '' is default and will be replaced by tick value when drawing!!
        self.axes.set_xticklabels(a)

    # -----------------------------------------------------------------

    def hide_last_ytick_label(self):

        """
        This function ...
        :return:
        """

        # Get labels
        #labels = self.ytick_labels

        # Hide last
        #labels[-1] = ""
        #self.set_ytick_labels(labels)

        # NEW
        a = self.axes.get_yticks().tolist()
        a[-1] = " "  # SPACE AND NOT EMPTY!!! '' is default and will be replaced by tick value when drawing!!
        self.axes.set_yticklabels(a)

    # -----------------------------------------------------------------

    def hide_xaxis(self):

        """
        This function ...
        :return:
        """

        # Not visible
        self.xaxis.set_visible(False)

        # No ticks
        self.xaxis.set_ticks([])

    # -----------------------------------------------------------------

    def hide_yaxis(self):

        """
        This function ...
        :return:
        """

        # Not visible
        self.yaxis.set_visible(False)

        # No ticks
        self.yaxis.set_ticks([])

    # -----------------------------------------------------------------

    def set_grid(self, config, which="major"):

        """
        This function ...
        :param config:
        :param which:
        :return:
        """

        if config.add_grid: self.axes.grid(linewidth=config.grid_linewidth, linestyle=config.grid_linestyle, color=config.grid_color, which=which)

    # -----------------------------------------------------------------

    @property
    def xticks(self):
        return self.axes.get_xticklabels()

    # -----------------------------------------------------------------

    @property
    def xtick_labels(self):
        return [tick.get_text() for tick in self.xticks]

    # -----------------------------------------------------------------

    @property
    def yticks(self):
        return self.axes.get_yticklabels()

    # -----------------------------------------------------------------

    @property
    def ytick_labels(self):
        return [tick.get_text() for tick in self.yticks]

    # -----------------------------------------------------------------

    def set_xaxis_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        self.xaxis.set_ticks_position(position)
        self.xaxis.set_label_position(position)

    # -----------------------------------------------------------------

    def set_xaxis_top(self):
        self.set_xaxis_position("top")

    # -----------------------------------------------------------------

    def set_xaxis_bottom(self):
        self.set_xaxis_position("bottom")

    # -----------------------------------------------------------------

    def set_yaxis_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        self.yaxis.set_ticks_position(position)
        self.yaxis.set_label_position(position)

    # -----------------------------------------------------------------

    def set_yaxis_left(self):
        self.set_yaxis_position("left")

    # -----------------------------------------------------------------

    def set_yaxis_right(self):
        self.set_yaxis_position("right")

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

        # Properties
        self.transparent = False
        self.format = None

    # -----------------------------------------------------------------

    @property
    def ax(self):
        return self.figure.gca()

    # -----------------------------------------------------------------

    def create_one_plot(self, projection=None):

        """
        Thisf unction ...
        :param projection:
        :return:
        """

        # Create axes
        if projection is not None: axes = self.figure.add_subplot(111, projection=projection)
        else: axes = self.ax

        # Create
        return MPLPlot(plot=axes)

    # -----------------------------------------------------------------

    def create_twin_plots(self, share="x"):

        """
        This function ...
        :param share:
        :return:
        """

        # Create first plot
        plot0 = self.create_one_plot()

        # Create second plot
        if share == "x": ax1 = plot0.axes.twinx()
        elif share == "y": ax1 = plot0.axes.twiny()
        else: raise ValueError("Invalid value for 'share'")
        plot1 = MPLPlot(plot=ax1)

        # Return plots
        return plot0, plot1

    # -----------------------------------------------------------------

    def add_plot(self, left, bottom, width, height):

        """
        This function ...
        :param left:
        :param bottom:
        :param width:
        :param height:
        :return:
        """

        axes = self.figure.add_axes([left, bottom, width, height])
        return MPLPlot(plot=axes)

    # -----------------------------------------------------------------

    def create_axes(self, left, bottom, width, height):
        return self.figure.add_axes([left, bottom, width, height])

    # -----------------------------------------------------------------

    def add_colorbar(self, left, bottom, width, height, cmap, orientation, interval, logscale=False, ticks=None):

        """
        This function ...
        :param left:
        :param bottom:
        :param width:
        :param height:
        :param cmap:
        :param orientation:
        :param interval:
        :param logscale:
        :param ticks:
        :return:
        """

        # Create axes
        axes = self.create_axes(left, bottom, width, height)

        # Create cb
        return self.create_colorbar(axes, cmap, orientation, interval, logscale=logscale, ticks=ticks)

    # -----------------------------------------------------------------

    def create_colorbar(self, axes, cmap, orientation, interval, logscale=False, ticks=None):

        """
        This function ...
        :param axes:
        :param cmap:
        :param orientation:
        :param interval:
        :param logscale:
        :param ticks:
        :return:
        """

        # Create norm
        if logscale: norm = LogNorm(vmin=interval[0], vmax=interval[1])
        else: norm = Normalize(vmin=interval[0], vmax=interval[1])

        # Create colorbar
        cb = ColorbarBase(axes, cmap=cmap, norm=norm, orientation=orientation)

        # Set ticks
        if ticks is not None: cb.set_ticks(ticks)

        # Return
        return cb

    # -----------------------------------------------------------------

    def create_vertical_colorbar_in_plot(self, plot, plotted=None, position="right", relwidth=0.05, relheight=0.8,
                                         cmap=None, norm=None, ticks=None, label=None, logarithmic=False, x_shift=0,
                                         y_shift=0):

        """
        This function ...
        :param plot:
        :param plotted:
        :param position:
        :param relwidth:
        :param relheight:
        :param cmap:
        :param norm:
        :param ticks:
        :param label:
        :param logarithmic:
        :param x_shift:
        :param y_shift:
        :return:
        """

        # Centered in height
        bottomtop_margin = 0.5 * (1 - relheight)
        bottom = plot.y_min + (bottomtop_margin + y_shift) * plot.height

        # Width
        if position == "left": left = plot.x_min + (0.1 + x_shift) * plot.width
        elif position == "right": left = plot.x_min + (0.9 - relwidth + x_shift) * plot.width
        else: raise ValueError("Invalid position: '" + position + "'")

        # Set size
        width = plot.width * relwidth
        height = plot.height * relheight

        # Create axes for colorbar
        cb_axes = self.create_axes(left, bottom, width, height)

        # Create
        if plotted is not None:
            if logarithmic:  # something went wrong when from logarithmic auxilary axis scatter density plot (Astrofrog implementation)
                norm = LogNorm(vmin=plotted.norm.vmin, vmax=plotted.norm.vmax)
                cb = ColorbarBase(cb_axes, cmap=plotted.cmap, norm=norm, orientation="vertical")
            else: cb = colorbar_factory(cb_axes, plotted, orientation="vertical")
        else:
            if cmap is None: raise ValueError("If plotted is not passed, colormap must be passed")
            if norm is None: raise ValueError("If plotted is not passed, normalization must be passed")
            cb = ColorbarBase(cb_axes, cmap=cmap, norm=norm, orientation="vertical")

        # Add label
        if label is not None:
            x_shift = 0.12
            plot.add_text(label, vertical_position=position, x_shift=x_shift, rotation=-45)

        # Set ticks
        if ticks is not None: cb.set_ticks(ticks)

        # Return
        return cb

    # -----------------------------------------------------------------

    def create_horizontal_colorbar_in_plot(self, plot, plotted=None, position="bottom", relwidth=0.8, relheight=0.05,
                                            cmap=None, norm=None, ticks=None, label=None, logarithmic=False, label_position=None,
                                            limits=None):

        """
        This function ...
        :param plot:
        :param plotted:
        :param position:
        :param relwidth:
        :param relheight:
        :param cmap:
        :param norm:
        :param ticks:
        :param label:
        :param logarithmic:
        :param label_position:
        :param limits:
        :return:
        """

        ## Set anchor

        # Centered in width
        leftright_margin = 0.5 * (1 - relwidth)
        left = plot.x_min + leftright_margin * plot.width

        # Height
        if position == "bottom": bottom = plot.y_min + 0.1 * plot.height
        elif position == "top": bottom = plot.y_min + (0.9 - relheight) * plot.height
        else: raise ValueError("Invalid position: '" + position + "'")

        # Set size
        width = plot.width * relwidth
        height = plot.height * relheight

        # Create axes for colorbar
        cb_axes = self.create_axes(left, bottom, width, height) #self.add_colorbar(left, bottom, width, height, "RdBu", "horizontal", (-1, 1), ticks=[-1, -0.5, 0, 0.5, 1])

        # Create
        if plotted is not None:
            if logarithmic: # something went wrong when from logarithmic auxilary axis scatter density plot (Astrofrog implementation)
                #cmap = "inferno"
                #print(plotted.cmap) if hasattr(plotted, "cmap") else print("nope")
                print("log colorbar", plotted.norm.vmin, plotted.norm.vmax)
                norm = LogNorm(vmin=plotted.norm.vmin, vmax=plotted.norm.vmax)
                cb = ColorbarBase(cb_axes, cmap=plotted.cmap, norm=norm, orientation="horizontal")
            elif plotted.norm.vmin is not None and plotted.norm.vmax is not None:
                print("linear colorbar", plotted.norm.vmin, plotted.norm.vmax)
                cb = colorbar_factory(cb_axes, plotted, orientation="horizontal")
                #norm = Normalize(vmin=plotted.norm.vmin, vmax=plotted.norm.vmax)
                #cb = ColorbarBase(cb_axes, cmap=plotted.cmap, norm=norm, orientation="horizontal")
            elif limits is not None:
                print("LIMITS", limits)
                norm = Normalize(vmin=limits[0], vmax=limits[1])
                cb = ColorbarBase(cb_axes, cmap=plotted.cmap, norm=norm, orientation="horizontal")
            else:
                log.warning("Limits of linear scale are not defined")
                cb = colorbar_factory(cb_axes, plotted, orientation="horizontal")
        else:
            if cmap is None: raise ValueError("If plotted is not passed, colormap must be passed")
            if norm is None: raise ValueError("If plotted is not passed, normalization must be passed")
            print("colorbar", norm.vmin, norm.vmax)
            cb = ColorbarBase(cb_axes, cmap=cmap, norm=norm, orientation="horizontal")

        # Add label
        if label is not None:
            y_shift = 0.12
            if label_position is None: plot.add_text(label, vertical_position=position, y_shift=y_shift)
            elif label_position == "right":
                x = left + width*2
                if position == "top":
                    y = 0.95 + y_shift
                    alignment = "top"
                elif position == "bottom":
                    y = 0.05 + y_shift
                    alignment = "bottom"
                else: raise ValueError("Invalid position: '" + position + "'")
                print("LABEL X (right)", x, left, width)
                plot.text(x, y, label, horizontalalignment='right', verticalalignment=alignment, transform=plot.axes.transAxes)
            elif label_position == "left":
                x = left
                if position == "top":
                    y = 0.95 + y_shift
                    alignment = "top"
                elif position == "bottom":
                    y = 0.05 + y_shift
                    alignment = "bottom"
                else: raise ValueError("Invalid position: '" + position + "'")
                print("LABEL X (left)", x, left)
                plot.text(x, y, label, horizontalalignment='left', verticalalignment=alignment, transform=plot.axes.transAxes)
            else: raise ValueError("Invalid value for 'label_position'")

        # Set scale
        #if logarithmic:
        #    print(plotted.norm)
        #    print(plotted.norm.vmin)
        #    print(plotted.norm.vmax)
        #    #cb.locator = LogLocator()
        #    #cb.formatter = LogFormatterMathtext()

        # Set ticks
        if ticks is not None: cb.set_ticks(ticks)

        # Return
        return cb

    # -----------------------------------------------------------------

    def create_column(self, size, share_axis=False, height_ratios=None, x_label=None, x_label_fontsize="small",
                      x_labels=None, y_labels=None, y_label_fontsize="small", x_scale="linear", x_scales=None, y_scales=None,
                      x_limits=None, y_limits=None, x_log_scalar=False, y_log_scalar=False, projections=None, hspace=None):

        """
        This function ...
        :param size:
        :param share_axis:
        :param height_ratios:
        :param x_label:
        :param x_label_fontsize:
        :param x_labels:
        :param y_labels:
        :param y_label_fontsize:
        :param x_scale:
        :param x_scales
        :param y_scales:
        :param x_limits:
        :param y_limits:
        :param x_log_scalar:
        :param y_log_scalar:
        :param projections:
        :param hspace:
        :return:
        """

        from ..tools import types

        # Define hspace
        if hspace is None:
            if share_axis: hspace = 0.0
            else: hspace = 0.05

        # Make the grid
        gs = GridSpec(size, 1, height_ratios=height_ratios, hspace=hspace)

        # Set projections: same for all plots?
        if projections is not None and not types.is_sequence(projections): projections = [projections for _ in range(size)]

        # Create the (sub)plots
        plots = []
        if share_axis:

            if x_labels is not None: raise ValueError("Cannot specify different x labels when sharing axis")
            if x_scales is not None: raise ValueError("Cannot specify different x scales when sharing axis")

            # Get first projection
            if projections is not None: projection = projections[0]
            else: projection = None

            # Create first plot
            first_mpl_plot = self.figure.add_subplot(gs[0], projection=projection)
            first_plot = MPLPlot(plot=first_mpl_plot)

            # Hide x axis for all but last plot (if sharing)
            # NO: hide x tick labels only
            first_plot.hide_xtick_labels()

            if y_labels is not None:
                label = y_labels[0]
                if label is not None: first_plot.set_ylabel(label, fontsize=y_label_fontsize)
                else: first_plot.hide_yaxis()

            if y_scales is not None:
                scale = y_scales[0]
                first_plot.set_yscale(scale)

            # Scalar?
            if y_log_scalar: first_plot.set_y_scalar()

            # Set shared x scale
            first_plot.set_xscale(x_scale)

            # Scalar?
            if x_log_scalar: first_plot.set_x_scalar()

            if x_limits is not None:
                lower, upper = x_limits
                first_plot.set_xlim(lower, upper)

            if y_limits is not None:
                limits = y_limits[0]
                first_plot.set_ylim(*limits)

            # Hide x axis
            #first_plot.hide_xticks()

            # Add the first plot
            plots.append(first_plot)

            # Create next plots
            for index in range(1, size):
                is_last = index == size - 1

                # Get projection
                if projections is not None: projection = projections[index]
                else: projection = None

                # Create next plot
                next_plot = self.figure.add_subplot(gs[index], sharex=first_mpl_plot, projection=projection)
                next_plot = MPLPlot(plot=next_plot)

                # Hide x axis for all but last plot (if sharing)
                # NO: hide x tick labels only
                #if not is_last: next_plot.hide_xaxis()
                if not is_last: next_plot.hide_xtick_labels()

                # Set y label
                if y_labels is not None:
                    label = y_labels[index]
                    if label is not None: next_plot.set_ylabel(label, fontsize=y_label_fontsize)
                    else: next_plot.hide_yaxis()

                # Set y scale
                if y_scales is not None:
                    scale = y_scales[index]
                    next_plot.set_yscale(scale)

                # Scalar?
                if y_log_scalar: next_plot.set_y_scalar()

                # Set y limits
                if y_limits is not None:
                    limits = y_limits[index]
                    next_plot.set_ylim(*limits)

                # Is this not last plot? -> hide x labels
                if index != size - 1:
                    next_plot.hide_xticks()

                # Add the plot
                plots.append(next_plot)

            last_plot = plots[-1]

            # Set shared x label
            if x_label is not None:
                last_plot.set_xlabel(x_label, fontsize=x_label_fontsize)

            # Set shared x scale
            #last_plot.set_xscale(x_scale)

        # Not sharing
        else:

            if x_label is not None: raise ValueError("Cannot specify one x label when not sharing axis")
            #if x_scale is not None:

            # Create plots
            for index in range(size):

                # Get projection
                if projections is not None: projection = projections[index]
                else: projection = None

                # Get plot
                plot = self.figure.add_subplot(gs[index], projection=projection)
                plot = MPLPlot(plot=plot)

                if x_labels is not None:
                    label = x_labels[index]
                    if label is not None: plot.set_xlabel(label, fontsize=x_label_fontsize)
                    else: plot.hide_xaxis()

                if y_labels is not None:
                    label = y_labels[index]
                    if label is not None: plot.set_ylabel(label, fontsize=y_label_fontsize)
                    else: plot.hide_yaxis()

                if x_scales is not None:
                    scale = x_scales[index]
                    plot.set_xscale(scale)

                # Scalar?
                if x_log_scalar: plot.set_x_scalar()

                if y_scales is not None:
                    scale = y_scales[index]
                    plot.set_yscale(scale)

                # Scalar?
                if y_log_scalar: plot.set_y_scalar()

                if x_limits is not None:
                    limits = x_limits[index]
                    plot.set_xlim(*limits)

                if y_limits is not None:
                    limits = y_limits[index]
                    plot.set_ylim(*limits)

                # Add the plot
                plots.append(plot)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_sed_plots(self, nresiduals=1):

        """
        This function ...
        :param nresiduals:
        :return:
        """

        # No residual panel?
        if nresiduals == 0: return self.create_one_plot(), []

        # Set subplot height ratios
        nplots = nresiduals + 1
        height_ratios = [4] + [1] * nresiduals

        # Create plots
        plots = self.create_column(nplots, share_axis=True, height_ratios=height_ratios)

        # Split
        main_plot = plots[0]
        residual_plots = plots[1:]

        # Return
        return main_plot, residual_plots

    # -----------------------------------------------------------------

    def create_row_of_sed_plots(self, npanels, nresiduals=1):

        """
        This function ...
        :param npanels:
        :param nresiduals:
        :return:
        """

        # Main plots
        main_plots = []
        residual_plots = []

        # The number of residual panels is the same for each SED panel
        if isinstance(nresiduals, int):

            # Set subplot height ratios
            nplots = nresiduals + 1 # nplots per column (nrows)
            height_ratios = [4] + [1] * nresiduals

            # Create plots
            plots = self.create_grid(nplots, npanels, height_ratios=height_ratios, sharex=True, sharey=True)

            # Set main plots and residual plots
            for index in range(npanels):
                main_plots.append(plots[0][index])
                residual_plots_panel = [row[index] for row in plots[1:]]
                residual_plots.append(residual_plots_panel)

        # Different number of residual panels for each SED panel? -> nresiduals is a list of integers
        else:

            max_nresiduals = max(nresiduals)
            nrows = max_nresiduals + 1
            height_ratios = [4] + [1] * max_nresiduals

            # Define space
            wspace = 0.0
            hspace = 0.0

            # Create grid spec
            grid = self.create_gridspec(nrows, npanels, wspace=wspace, hspace=hspace, height_ratios=height_ratios)

            # Main plots
            main_plots = []
            residual_plots = []

            # Loop over the SED panels
            for index in range(npanels):

                # Get the number of residual panels for this SED panel
                nres = nresiduals[index]
                nrows_for_main = nrows - nres

                # Create the main plot
                #rect = grid[:-nres, index]
                rect = grid[:nrows_for_main, index]

                # Create plot
                #main_plot = self._create_plot_not_shared(rect)
                if index == 0: main_plot = self._create_plot_not_shared(rect)
                else: main_plot = self._create_plot_shared_y(rect, main_plots[0].axes)
                main_plots.append(main_plot)

                # Create the residual plots
                res_plots = []
                for j in range(nres):
                    resrect = grid[nrows_for_main+j, index]
                    #resplot = self._create_plot_not_shared(resrect)
                    resplot = self._create_plot_shared_x(resrect, main_plot.axes)
                    #resplot = self._create_plot_not_shared(resrect)
                    res_plots.append(resplot)
                residual_plots.append(res_plots)

        # Return
        return main_plots, residual_plots

    # -----------------------------------------------------------------

    def create_column_of_sed_plots(self):
        # not implemented yet
        pass

    # -----------------------------------------------------------------

    def create_row(self, size, share_axis=False, width_ratios=None, y_label=None, y_label_fontsize="small",
                   y_labels=None, x_labels=None, x_label_fontsize="small", y_scale="linear", y_scales=None,
                   x_scales=None, y_limits=None, x_limits=None, y_log_scalar=False, x_log_scalar=False,
                   projections=None, wspace=None):

        """
        This function ...
        :param size:
        :param share_axis:
        :param width_ratios:
        :param y_label:
        :param y_label_fontsize:
        :param y_labels:
        :param x_labels:
        :param x_label_fontsize:
        :param y_scale:
        :param y_scales:
        :param x_scales:
        :param y_limits:
        :param x_limits:
        :param y_log_scalar:
        :param x_log_scalar:
        :param projections:
        :param wspace:
        :return:
        """

        from ..tools import types

        # Define wspace
        if wspace is None:
            if share_axis: wspace = 0.0
            else: wspace = 0.05

        # Make the grid
        gs = GridSpec(1, size, width_ratios=width_ratios, wspace=wspace)

        # Set projections: same for all plots?
        if projections is not None and not types.is_sequence(projections): projections = [projections for _ in range(size)]

        # Create the (sub)plots
        plots = []
        if share_axis:

            if y_labels is not None: raise ValueError("Cannot specify different y labels when sharing axis")
            if y_scales is not None: raise ValueError("Cannot specify different y scales when sharing axis")

            # Get the first projection
            if projections is not None: projection = projections[0]
            else: projection = None

            # Create the first plot
            first_mpl_plot = self.figure.add_subplot(gs[0], projection=projection)
            first_plot = MPLPlot(plot=first_mpl_plot)

            if x_labels is not None:
                label = x_labels[0]
                if label is not None: first_plot.set_xlabel(label, fontsize=x_label_fontsize)
                else: first_plot.hide_xaxis()

            if x_scales is not None:
                scale = x_scales[0]
                first_plot.set_xscale(scale)

            # Scalar?
            if y_log_scalar: first_plot.set_y_scalar()

            # Set shared y scale
            first_plot.set_yscale(y_scale)

            # Set shared y label
            # if y_label is not None: last_plot.set_ylabel(y_label, fontsize=y_label_fontsize)
            if y_label is not None: first_plot.set_ylabel(y_label, fontsize=y_label_fontsize)

            # Scalar?
            if y_log_scalar: first_plot.set_y_scalar()

            if y_limits is not None:
                lower, upper = y_limits
                first_plot.set_ylim(lower, upper)

            if x_limits is not None:
                limits = x_limits[0]
                first_plot.set_xlim(*limits)

            # Add the first plot
            plots.append(first_plot)

            # Create next plots
            for index in range(1, size):

                # Get the next projection
                if projections is not None: projection = projections[index]
                else: projection = None

                # Create the next plot
                next_plot = self.figure.add_subplot(gs[index], sharey=first_mpl_plot, projection=projection)
                next_plot = MPLPlot(plot=next_plot)

                # Hide y tick labels for all but first plot (if sharing)
                next_plot.hide_ytick_labels()

                # Set x label
                if x_labels is not None:
                    label = x_labels[index]
                    if label is not None: next_plot.set_xlabel(label, fontsize=x_label_fontsize)
                    else: next_plot.hide_xaxis()

                # Set x scale
                if x_scales is not None:
                    scale = x_scales[index]
                    next_plot.set_xscale(scale)

                # Scalar?
                if x_log_scalar: next_plot.set_x_scalar()

                # Set x limits
                if x_limits is not None:
                    limits = x_limits[index]
                    next_plot.set_xlim(*limits)

                # Add the next plot
                plots.append(next_plot)

            last_plot = plots[-1]

            # Set shared y label
            #if y_label is not None: last_plot.set_ylabel(y_label, fontsize=y_label_fontsize)

        # Not sharing axis
        else:

            if y_label is not None: raise ValueError("Cannot specify one y label when not sharing axis")

            # Create plots
            for index in range(size):

                # Get projection
                if projections is not None: projection = projections[index]
                else: projection = None

                # Get plot
                plot = self.figure.add_subplot(gs[index], projection=projection)
                plot = MPLPlot(plot=plot)

                if y_labels is not None:
                    label = y_labels[index]
                    if label is not None: plot.set_ylabel(label, fontsize=y_label_fontsize)
                    else: plot.hide_yaxis()

                if x_labels is not None:
                    label = x_labels[index]
                    if label is not None: plot.set_xlabel(label, fontsize=x_label_fontsize)
                    else: plot.hide_xaxis()

                if y_scales is not None:
                    scale = y_scales[index]
                    plot.set_yscale(scale)

                # Scalar?
                if y_log_scalar: plot.set_y_scalar()

                if x_scales is not None:
                    scale = x_scales[index]
                    plot.set_xscale(scale)

                # Scalar?
                if x_log_scalar: plot.set_x_scalar()

                if y_limits is not None:
                    limits = y_limits[index]
                    plot.set_ylim(*limits)

                if x_limits is not None:
                    limits = x_limits[index]
                    plot.set_xlim(*limits)

                # Add the plot
                plots.append(plot)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_row_of_image_grids(self, nrows, ncols, ngrids, wspace=0.0, hspace=0.0, colorbar_relsize=0.05,
                                  return_colorbars=False, share_colorbars=True, adjust_grid=False):

        """
        This function ...
        :param nrows:
        :param ncols:
        :param ngrids:
        :param wspace:
        :param hspace:
        :param colorbar_relsize:
        :param return_colorbars:
        :param share_colorbars:
        :param adjust_grid:
        :return:
        """

        # Set colorbar size in percentage
        cbar_size = str(colorbar_relsize * 100) + "%"

        # Set the colorbar mode
        if share_colorbars: cbar_mode = "single"
        else: cbar_mode = "edge"

        axes_pad = (wspace, hspace)

        axes_class = None
        label_mode = "L"

        # No axes for the main figure
        self.ax.set_axis_off()

        # The list of grids
        grids = []

        # IMPORTANT
        share_all = True

        # SET ASPECT
        if adjust_grid: aspect = True # default
        else: aspect = False # don't change the grid based on the shape of the images

        cbar_set_cax = False
        #cbar_set_cax = True

        # Create the grids
        for index in range(ngrids):

            # Determine subplot specified
            rect = "1" + str(ngrids) + str(index+1)

            # Create the grid
            grid = ImageGrid(self.figure, rect, nrows_ncols=(nrows, ncols), axes_pad=axes_pad, aspect=aspect,
                             cbar_mode=cbar_mode, add_all=True, cbar_set_cax=cbar_set_cax, cbar_size=cbar_size,
                             axes_class=axes_class, label_mode=label_mode, share_all=share_all)
            grids.append(grid)

        # Initialize structure to contain the plots
        plots = [[[None for i in range(ncols)] for j in range(nrows)] for k in range(ngrids)]

        # Initialize structure to contain the colorbar axes
        colorbars = [[[None for i in range(ncols)] for j in range(nrows)] for k in range(ngrids)]

        # Loop over the images
        for i in range(ngrids):

            index = 0

            grid = grids[i]

            # Create the plots
            for row in range(nrows):
                for col in range(ncols):

                    # Get axes, create subplot?
                    ax = grid[index]
                    plot = ax

                    # last column
                    if col == ncols-1:
                        ax._sharex = None
                        #ax._shared_x_axes = cbook.Grouper()
                        #ax._adjustable = 'box'
                        ax._sharey = None
                        ax._shared_y_axes = cbook.Grouper()
                        ax._adjustable = 'box'

                    # Create plot
                    plot = MPLPlot(plot=plot)

                    # Add the plot
                    plots[i][row][col] = plot

                    # Get the colorbar axes
                    colorbar = grid.cbar_axes[index]

                    # Add the colorbar
                    colorbars[i][row][col] = colorbar

                    # Increment the index
                    index += 1

        # Return the plots (and errorbars)
        if return_colorbars: return plots, colorbars
        else: return plots

    # -----------------------------------------------------------------

    def create_image_grid(self, nrows, ncols, wspace=0.0, hspace=0.0, return_colorbar=False, colorbar_relsize=0.05,
                          edgecolor=None, projection=None, adjust_grid=False):

        """
        This function ...
        :param nrows:
        :param ncols:
        :param wspace:
        :param hspace:
        :param return_colorbar:
        :param colorbar_relsize:
        :param edgecolor:
        :param projection:
        :param adjust_grid:
        :return:
        """

        # Set colorbar size in percentage
        cbar_size = str(colorbar_relsize*100) + "%"

        cbar_mode = "single"
        axes_pad = (wspace, hspace)

        # No axes for the main figure
        self.ax.set_axis_off()

        #if edgecolor is not None:
        #    self.ax.spines['bottom'].set_color("white")
        #    self.ax.spines['top'].set_color("white")
        #    self.ax.spines['left'].set_color("white")
        #    self.ax.spines['right'].set_color("white")

        #axes_class, axes_class_args = axes_class

        # Set axes class and projection
        #axes_class = ImageGrid._defaultLocatableAxesClass
        #axes_kwargs = {}
        #if projection is not None: axes_kwargs["projection"] = projection
        #axes_class = (axes_class, axes_kwargs)

        if projection is not None:

            from mpl_toolkits.axes_grid1.axes_divider import locatable_axes_factory
            projection_class, extra_kwargs = projection._as_mpl_axes()
            loc_axes_class = locatable_axes_factory(projection_class)
            axes_class = (loc_axes_class, extra_kwargs)

            label_mode = None #### I ADAPTED THE IMAGEGRID CLASS SO TO AVOID A CRASH WITH WCS AXES

        else:
            axes_class = None
            label_mode = "L"

        # IMPORTANT
        share_all = True

        # SET ASPECT
        if adjust_grid: aspect = True  # default
        else: aspect = False  # don't change the grid based on the shape of the images

        # Create the grid
        grid = ImageGrid(self.figure, 111, nrows_ncols=(nrows, ncols), axes_pad=axes_pad, aspect=aspect,
                         cbar_mode=cbar_mode, add_all=True, cbar_set_cax=False, cbar_size=cbar_size,
                         axes_class=axes_class, label_mode=label_mode, share_all=share_all)

        if projection is not None:

            # SET LABEL MODE, FOR WCS PROJECTIONS!!! BUG???
            def _tick_only(ax, bottom_on, left_on):

                bottom_off = not bottom_on
                left_off = not left_on
                # [l.set_visible(bottom_off) for l in ax.get_xticklabels()]
                # [l.set_visible(left_off) for l in ax.get_yticklabels()]
                # ax.xaxis.label.set_visible(bottom_off)
                # ax.yaxis.label.set_visible(left_off)

                #xaxis = ax.xaxis
                #yaxis = ax.yaxis

                #print(type(xaxis), vars(xaxis))
                #for label in vars(xaxis):
                    #if "label" in label: print(label, getattr(xaxis, label))

                # OUTPUT:
                # _label
                # isDefault_label True
                # label Text(0.5,0,u'')
                # labelpad 4.0
                # _autolabelpos True
                # label_position bottom

                try:
                    # Original, also works
                    # ax.axis["bottom"].toggle(ticklabels=bottom_off, label=bottom_off)
                    # ax.axis["left"].toggle(ticklabels=left_off, label=left_off)
                    # Try
                    ax.xaxis.toggle(ticklabels=bottom_off, label=bottom_off)
                    ax.yaxis.toggle(ticklabels=left_off, label=left_off)
                except:
                    xaxis = ax.xaxis
                    yaxis = ax.yaxis
                    xaxis.label = None
                    xaxis.label_position = None
                    yaxis.label = None
                    yaxis.label_position = None

            # left-most axes
            for ax in grid.axes_column[0][:-1]:
                _tick_only(ax, bottom_on=True, left_on=False)

            # lower-left axes
            ax = grid.axes_column[0][-1]
            _tick_only(ax, bottom_on=False, left_on=False)

            for col in grid.axes_column[1:]:
                # axes with no labels
                for ax in col[:-1]: _tick_only(ax, bottom_on=True, left_on=True)

                # bottom
                ax = col[-1]
                _tick_only(ax, bottom_on=False, left_on=True)

        #print(grid._divider)
        #subplotspec = grid._divider._subplotspec
        #for label in vars(subplotspec): print(label, vars(subplotspec)[label])
        #gs = subplotspec._gridspec
        #print(gs[0, 0])
        #print(gs[1, 0])
        #print(gs[0, 1])

        # Set the color map axes
        colorbar = grid.cbar_axes[0]
        #grid[0].cax = colorbar # set colorbar axes explicitly

        # Initialize structure to contain the plots
        plots = [[None for i in range(ncols)] for j in range(nrows)]

        # Loop over the images
        index = 0
        for row in range(nrows):
            for col in range(ncols):

                # Get axes, create subplot?
                ax = grid[index]
                plot = ax

                # Create plot
                plot = MPLPlot(plot=plot)

                # Add the plot
                plots[row][col] = plot

                index += 1

        # Return the plots
        if return_colorbar: return plots, colorbar
        else: return plots

    # -----------------------------------------------------------------

    def create_row_of_grids(self, nrows, ncols, ngrids, width_ratios=None, height_ratios=None, sharex=False, sharey=False,
                            projections=None, first_row_not_shared_x=None, first_row_not_shared_y=None, last_row_not_shared_x=None,
                            last_row_not_shared_y=None, first_column_not_shared_x=None, first_column_not_shared_y=None,
                            last_column_not_shared_x=None, last_column_not_shared_y=None, rows_shared_x=None, rows_shared_y=None,
                            columns_shared_x=None, columns_shared_y=None, share_per_row=True, share_per_column=True):

        """
        This function ...
        :param nrows:
        :param ncols:
        :param ngrids:
        :param width_ratios:
        :param height_ratios:
        :param sharex:
        :param sharey:
        :param projections:
        :param first_row_not_shared_x:
        :param first_row_not_shared_y:
        :param last_row_not_shared_x:
        :param last_row_not_shared_y:
        :param first_column_not_shared_x:
        :param first_column_not_shared_y:
        :param last_column_not_shared_x:
        :param last_column_not_shared_y:
        :param rows_shared_x:
        :param rows_shared_y:
        :param columns_shared_x:
        :param columns_shared_y:
        :param share_per_row:
        :param share_per_column:
        :return:
        """

        # The main grid spec
        main = GridSpec(1, ngrids)

        # Initialize a list for the plots
        plots = []

        # Create the grids
        for index in range(ngrids):

            # Get the projections
            if projections is not None: projections_grid = projections[index]
            else: projections_grid = None

            # Create the grid
            grid_plots = self.create_grid(nrows, ncols, projections=projections_grid, subplotspec=main[index],
                                          width_ratios=width_ratios, height_ratios=height_ratios, sharex=sharex, sharey=sharey,
                                          first_row_not_shared_x=first_row_not_shared_x, first_row_not_shared_y=first_row_not_shared_y,
                                          last_row_not_shared_x=last_row_not_shared_x, last_row_not_shared_y=last_row_not_shared_y,
                                          first_column_not_shared_x=first_column_not_shared_x, first_column_not_shared_y=first_column_not_shared_y,
                                          last_column_not_shared_x=last_column_not_shared_x, last_column_not_shared_y=last_column_not_shared_y,
                                          rows_shared_x=rows_shared_x, rows_shared_y=rows_shared_y, columns_shared_x=columns_shared_x,
                                          columns_shared_y=columns_shared_y, share_per_row=share_per_row, share_per_column=share_per_column)

            # Add the plots
            plots.append(grid_plots)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_column_of_grids(self, nrows, ncols, ngrids, width_ratios=None, height_ratios=None, sharex=False, sharey=False,
                            projections=None, first_row_not_shared_x=None, first_row_not_shared_y=None, last_row_not_shared_x=None,
                            last_row_not_shared_y=None, first_column_not_shared_x=None, first_column_not_shared_y=None,
                            last_column_not_shared_x=None, last_column_not_shared_y=None, rows_shared_x=None, rows_shared_y=None,
                            columns_shared_x=None, columns_shared_y=None, share_per_row=True, share_per_column=True):

        """
        This function ...
        :param nrows:
        :param ncols:
        :param ngrids:
        :param width_ratios:
        :param height_ratios:
        :param sharex:
        :param sharey:
        :param projections:
        :param first_row_not_shared_x:
        :param first_row_not_shared_y:
        :param last_row_not_shared_x:
        :param last_row_not_shared_y:
        :param first_column_not_shared_x:
        :param first_column_not_shared_y:
        :param last_column_not_shared_x:
        :param last_column_not_shared_y:
        :param rows_shared_x:
        :param rows_shared_y:
        :param columns_shared_x:
        :param columns_shared_y:
        :param share_per_row:
        :param share_per_column:
        :return:
        """

        # The main grid spec
        main = GridSpec(ngrids, 1)

        # Initialize a list for the plots
        plots = []

        # Create the grids
        for index in range(ngrids):

            # Get the projections
            if projections is not None: projections_grid = projections[index]
            else: projections_grid = None

            # Create the grid
            grid_plots = self.create_grid(nrows, ncols, projections=projections_grid, subplotspec=main[index],
                                          width_ratios=width_ratios, height_ratios=height_ratios, sharex=sharex, sharey=sharey,
                                          first_row_not_shared_x=first_row_not_shared_x, first_row_not_shared_y=first_row_not_shared_y,
                                          last_row_not_shared_x=last_row_not_shared_x, last_row_not_shared_y=last_row_not_shared_y,
                                          first_column_not_shared_x=first_column_not_shared_x, first_column_not_shared_y=first_column_not_shared_y,
                                          last_column_not_shared_x=last_column_not_shared_x, last_column_not_shared_y=last_column_not_shared_y,
                                          rows_shared_x=rows_shared_x, rows_shared_y=rows_shared_y, columns_shared_x=columns_shared_x,
                                          columns_shared_y=columns_shared_y, share_per_row=share_per_row, share_per_column=share_per_column)

            # Add the plots
            plots.append(grid_plots)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_grid(self, nrows, ncols, width_ratios=None, height_ratios=None, sharex=False, sharey=False,
                    projections=None, first_row_not_shared_x=None, first_row_not_shared_y=None,
                    last_row_not_shared_x=None, last_row_not_shared_y=None, first_column_not_shared_x=None, first_column_not_shared_y=None,
                    last_column_not_shared_x=None, last_column_not_shared_y=None, rows_shared_x=None, rows_shared_y=None,
                    columns_shared_x=None, columns_shared_y=None, subplotspec=None, share_per_row=True, share_per_column=True,
                    wspace=None, hspace=None, empty=None):

        """
        This function ...
        :param nrows:
        :param ncols:
        :param width_ratios:
        :param height_ratios:
        :param sharex:
        :param sharey:
        :param projections:
        :param first_row_not_shared_x: flag
        :param first_row_not_shared_y: flag
        :param last_row_not_shared_x: flag
        :param last_row_not_shared_y: flag
        :param first_column_not_shared_x: flag
        :param first_column_not_shared_y: flag
        :param last_column_not_shared_x: flag
        :param last_column_not_shared_y: flag
        :param rows_shared_x:
        :param rows_shared_y:
        :param columns_shared_x:
        :param columns_shared_y:
        :param subplotspec:
        :param share_per_row: share only the y axes per row
        :param share_per_column: share only the x axes per column
        :param wspace:
        :param hspace:
        :param empty:
        :return:
        """

        # Define wspace
        if wspace is None:
            if sharey: wspace = 0.0
            else: wspace = 0.05

        # Define hspace
        if hspace is None:
            if sharex: hspace = 0.0
            else: hspace = 0.05

        # Set width and height ratios
        if width_ratios is None: width_ratios = [1] * ncols
        if height_ratios is None: height_ratios = [1] * nrows

        # Set list of shared x flags for rows
        if rows_shared_x is None:
            rows_shared_x = []
            for row in range(nrows):
                if row == 0 and first_row_not_shared_x: shared = False
                elif row == nrows - 1 and last_row_not_shared_x: shared = False
                else: shared = sharex
                rows_shared_x.append(shared)
        elif first_row_not_shared_x is not None or last_row_not_shared_x is not None: raise ValueError("Cannot specify first_row_not_shared_x or last_row_not_shared_x")

        # Set list of shared y flags for rows
        if rows_shared_y is None:
            rows_shared_y = []
            for row in range(nrows):
                if row == 0 and first_row_not_shared_y: shared = False
                elif row == nrows - 1 and last_row_not_shared_y: shared = False
                else: shared = sharey
                rows_shared_y.append(shared)
        elif first_row_not_shared_y is not None or last_row_not_shared_y is not None: raise ValueError("Cannot specify first_row_not_shared_y or last_row_not_shared_y")

        # Set list of shared x flags for columns
        if columns_shared_x is None:
            columns_shared_x = []
            for col in range(ncols):
                if col == 0 and first_column_not_shared_x: shared = False
                elif col == ncols - 1 and last_column_not_shared_x: shared = False
                else: shared = sharex
                columns_shared_x.append(shared)
        elif first_column_not_shared_x is not None or last_column_not_shared_x is not None: raise ValueError("Cannot specify first_column_not_shared_x or last_column_not_shared_x")

        # Set list of shared y flags for columns
        if columns_shared_y is None:
            columns_shared_y = []
            for col in range(ncols):
                if col == 0 and first_column_not_shared_y: shared = False
                elif col == ncols - 1 and last_column_not_shared_y: shared = False
                else: shared = sharey
                columns_shared_y.append(shared)
        elif first_column_not_shared_y is not None or last_column_not_shared_y is not None: raise ValueError("Cannot specify first_column_not_shared_y or last_column_not_shared_y")

        # Create grid spec
        grid = self.create_gridspec(nrows, ncols, wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios, subplotspec=subplotspec)

        # Both axes are shared
        if sharex and sharey: plots = self._create_grid_shared(grid, nrows, ncols, projections=projections,
                                                               rows_shared_x=rows_shared_x, rows_shared_y=rows_shared_y,
                                                               columns_shared_x=columns_shared_x, columns_shared_y=columns_shared_y,
                                                               share_per_row=share_per_row, share_per_column=share_per_column, empty=empty)

        # Only x axis is shared
        elif sharex: plots = self._create_grid_shared_x(grid, nrows, ncols, projections=projections, empty=empty)

        # Only y axis is shared
        elif sharey: plots = self._create_grid_shared_y(grid, nrows, ncols, projections=projections, empty=empty)

        # No axes are shared
        else: plots = self._create_grid_not_shared(grid, nrows, ncols, projections=projections, empty=empty)

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_gridspec(self, nrows, ncols, wspace=None, hspace=None, width_ratios=None, height_ratios=None, subplotspec=None):

        """
        This function ...
        :param nrows:
        :param ncols:
        :param wspace:
        :param hspace:
        :param width_ratios:
        :param height_ratios:
        :param subplotspec:
        :return:
        """

        # Set width and height ratios
        if width_ratios is None: width_ratios = [1] * ncols
        if height_ratios is None: height_ratios = [1] * nrows

        # Create grid spec
        if subplotspec is not None: grid = GridSpecFromSubplotSpec(nrows, ncols, wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios, subplot_spec=subplotspec)
        else: grid = GridSpec(nrows, ncols, wspace=wspace, hspace=hspace, width_ratios=width_ratios, height_ratios=height_ratios)

        # Return the grid
        return grid

    # -----------------------------------------------------------------

    def _create_grid_plot_shared(self, grid, row, col, sharex, sharey, projections=None):

        """
        This function ...
        :param grid:
        :param row:
        :param col:
        :param sharex:
        :param sharey:
        :param projections:
        :return:
        """

        from ..tools import types

        # Get the grid spec
        rect = grid[row, col]

        # Get the projection
        if projections is not None:
            if types.is_sequence(projections) or types.is_dictionary(projections): projection = projections[row][col]
            else: projection = projections  # assume one projection for all plots
        else: projection = None

        # Create and return the plot
        return self._create_plot_shared(rect, sharex, sharey, projection=projection)

    # -----------------------------------------------------------------

    def _create_plot_shared(self, rect, sharex, sharey, projection=None):

        """
        This function ...
        :param rect:
        :param sharex:
        :param sharey:
        :param projection:
        :return:
        """

        # Get the axes
        plot = self.figure.add_subplot(rect, projection=projection, sharex=sharex, sharey=sharey)

        # Create plot
        plot = MPLPlot(plot=plot)

        # Return the plot
        return plot

    # -----------------------------------------------------------------

    def _create_grid_plot_shared_x(self, grid, row, col, sharex, projections=None):

        """
        This function ...
        :param grid:
        :param row:
        :param col:
        :param sharex:
        :param projections:
        :return:
        """

        from ..tools import types

        # Get the grid spec
        rect = grid[row, col]

        # Get the projection
        if projections is not None:
            if types.is_sequence(projections) or types.is_dictionary(projections): projection = projections[row][col]
            else: projection = projections  # assume one projection for all plots
        else: projection = None

        # Create and return the plot
        return self._create_plot_shared_x(rect, sharex, projection=projection)

    # -----------------------------------------------------------------

    def _create_plot_shared_x(self, rect, sharex, projection=None):

        """
        This function ...
        :param rect:
        :param sharex:
        :param projection:
        :return:
        """

        # Get the axes
        plot = self.figure.add_subplot(rect, projection=projection, sharex=sharex)

        # Create plot
        plot = MPLPlot(plot=plot)

        # Return the plot
        return plot

    # -----------------------------------------------------------------

    def _create_grid_plot_shared_y(self, grid, row, col, sharey, projections=None):

        """
        This function ...
        :param grid:
        :param row:
        :param col:
        :param sharey:
        :param projections:
        :return:
        """

        from ..tools import types

        # Get the grid spec
        rect = grid[row, col]

        # Get the projection
        if projections is not None:
            if types.is_sequence(projections) or types.is_dictionary(projections): projection = projections[row][col]
            else: projection = projections  # assume one projection for all plots
        else: projection = None

        # Create and return the plot
        return self._create_plot_shared_y(rect, sharey, projection=projection)

    # -----------------------------------------------------------------

    def _create_plot_shared_y(self, rect, sharey, projection=None):

        """
        This function ...
        :param rect:
        :param sharey:
        :param projection:
        :return:
        """

        # Get the axes
        plot = self.figure.add_subplot(rect, projection=projection, sharey=sharey)

        # Create plot
        plot = MPLPlot(plot=plot)

        # Return the plot
        return plot

    # -----------------------------------------------------------------

    def _create_grid_plot_not_shared(self, grid, row, col, projections=None):

        """
        This function ...
        :param grid:
        :param row:
        :param col:
        :param projections:
        :return:
        """

        from ..tools import types

        # Get the grid spec
        rect = grid[row, col]

        # Get the projection
        if projections is not None:
            if types.is_sequence(projections) or types.is_dictionary(projections): projection = projections[row][col]
            else: projection = projections  # assume one projection for all plots
        else: projection = None

        # Create and return the plot
        return self._create_plot_not_shared(rect, projection=projection)

    # -----------------------------------------------------------------

    def _create_plot_not_shared(self, rect, projection=None):

        """
        This function ...
        :param rect:
        :param projection:
        :return:
        """

        # Get the axes
        plot = self.figure.add_subplot(rect, projection=projection)

        # Create plot
        plot = MPLPlot(plot=plot)

        # Return the plot
        return plot

    # -----------------------------------------------------------------

    def _create_grid_shared(self, grid, nrows, ncols, projections=None, rows_shared_x=None, rows_shared_y=None,
                            columns_shared_x=None, columns_shared_y=None, share_per_row=True, share_per_column=True,
                            empty=None):

        """
        This function ...
        :param grid:
        :param nrows:
        :param ncols:
        :param projections:
        :param rows_shared_x:
        :param rows_shared_y:
        :param columns_shared_x:
        :param columns_shared_y:
        :param share_per_row:
        :param share_per_column:
        :param empty:
        :return:
        """

        # Initialize structure to contain the plots
        plots = [[None for i in range(ncols)] for j in range(nrows)]

        # Set indices
        #first_row = 0
        #first_column = 0
        #last_row = nrows - 1
        #last_column = ncols - 1

        # Initialize lists of reference rows and columns
        reference_rows = [] # reference row (for sharing x) per column
        reference_columns = [] # reference column (for sharing y) per row

        # Set reference rows (per column)
        for col in range(ncols):
            reference_row = None
            for row in range(nrows):
                shared_x = rows_shared_x[row] and columns_shared_x[col]
                if shared_x:
                    reference_row = row
                    break
            if reference_row is None: raise RuntimeError("Something went wrong")
            reference_rows.append(reference_row)

        # Set reference columns (per row)
        for row in range(nrows):
            reference_column = None
            for col in range(ncols):
                shared_y = rows_shared_y[row] and columns_shared_y[col]
                if shared_y:
                    reference_column = col
                    break
            if reference_column is None: raise RuntimeError("Something went wrong")
            reference_columns.append(reference_column)

        #print("reference rows:", reference_rows)
        #print("reference columns:", reference_columns)

        # Create the reference plots
        for col in reference_rows:
            row = reference_rows[col]
            if plots[row][col] is not None: continue # already created
            plot = self._create_grid_plot_not_shared(grid, row, col, projections=projections)
            #print("reference row plotted", row, col)
            plots[row][col] = plot
        for row in reference_columns:
            col = reference_columns[row]
            if plots[row][col] is not None: continue # already created
            plot = self._create_grid_plot_not_shared(grid, row, col, projections=projections)
            #print("reference col plotted", row, col)
            plots[row][col] = plot

        # Create the next rows
        for row in range(nrows):

            # Loop over the columns
            for col in range(ncols):

                # Already plotted
                if plots[row][col] is not None: continue

                # No plot necessary
                if empty is not None and (row,col,) in empty: continue

                #print("Creating plot " + str(row) + " " + str(col) + " ...")

                # Shared?
                shared_x = rows_shared_x[row] and columns_shared_x[col]
                shared_y = rows_shared_y[row] and columns_shared_y[col]

                #print("sharing x: " + str(shared_x))
                #print("sharing y: " + str(shared_y))

                # X and Y shared with other row and column
                if shared_x and shared_y:

                    reference_row = reference_rows[col]
                    reference_col = reference_columns[row]
                    #print("reference x", reference_row, col)
                    #print("reference y", row, reference_col)

                    # Share y only per row or for every row?
                    if share_per_row: reference_y = (row, reference_col)
                    else: reference_y = (reference_row, reference_col)

                    # Share x only per column or for every column?
                    if share_per_column: reference_x = (reference_row, col)
                    else: reference_x = (reference_row, reference_col)

                    # Is this plot it's own reference for either x or y?
                    is_reference_x = reference_x == (row, col)
                    is_reference_y = reference_y == (row, col)

                    #print("reference x: ", reference_x, is_reference_x)
                    #print("reference y: ", reference_y, is_reference_y)

                    # This plot is its own reference for both x and y: don't share, just create the plot
                    if is_reference_x and is_reference_y: plot = self._create_grid_plot_not_shared(grid, row, col, projections=projections)

                    # This plot is its own reference for x: share y
                    elif is_reference_x:

                        reference_plot = plots[reference_y[0]][reference_y[1]]
                        plot = self._create_grid_plot_shared_y(grid, row, col, reference_plot.axes, projections=projections)

                    # This plot is its own reference for y: share x
                    elif is_reference_y:

                        reference_plot = plots[reference_x[0]][reference_x[1]]
                        plot = self._create_grid_plot_shared_x(grid, row, col, reference_plot.axes, projections=projections)

                    # This plot is not a reference for either x or y
                    else:

                        reference_x_plot = plots[reference_x[0]][reference_x[1]]
                        reference_y_plot = plots[reference_y[0]][reference_y[1]]
                        plot = self._create_grid_plot_shared(grid, row, col, reference_x_plot.axes, reference_y_plot.axes, projections=projections)

                # X shared with other row
                elif shared_x:

                    reference_row = reference_rows[col]

                    #if share_per_column: reference = (reference_row, col)
                    #else: reference = (reference_row, reference_col)
                    reference = (reference_row, col)
                    is_reference = reference == (row, col)

                    # This plot is its own reference
                    if is_reference: plot = self._create_grid_plot_not_shared(grid, row, col, projections=projections)

                    # This plot is not its own reference
                    else:

                        reference_plot = plots[reference[0]][reference[1]]
                        plot = self._create_grid_plot_shared_x(grid, row, col, reference_plot.axes, projections=projections)

                # Y shared with other columns
                elif shared_y:

                    reference_col = reference_columns[row]

                    reference = (row, reference_col)
                    is_reference = reference == (row, col)

                    # This plot is its own reference
                    if is_reference: plot = self._create_grid_plot_not_shared(grid, row, col, projections=projections)

                    # This plot is not its own reference
                    else:

                        reference_plot = plots[reference[0]][reference[1]]
                        plot = self._create_grid_plot_shared_y(grid, row, col, reference_plot.axes, projections=projections)

                # No axis shared
                else: plot = self._create_grid_plot_not_shared(grid, row, col, projections=projections)

                # Create the plot and add it
                plots[row][col] = plot

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def _create_grid_shared_x(self, grid, nrows, ncols, projections=None, empty=None):

        """
        This function ...
        :param grid:
        :param nrows:
        :param ncols:
        :param projections:
        :param empty:
        :return:
        """

        # Initialize structure to contain the plots
        plots = [[None for i in range(ncols)] for j in range(nrows)]

        # Not yet implemented
        raise NotImplementedError("Not yet implemented")

    # -----------------------------------------------------------------

    def _create_grid_shared_y(self, grid, nrows, ncols, projections=None, empty=None):

        """
        This function ...
        :param grid:
        :param nrows:
        :param ncols:
        :param projections:
        :param empty:
        :return:
        """

        # Initialize structure to contain the plots
        plots = [[None for i in range(ncols)] for j in range(nrows)]

        # Not yet implemented
        raise NotImplementedError("Not yet implemented")

    # -----------------------------------------------------------------

    def _create_grid_not_shared(self, grid, nrows, ncols, projections=None, empty=None):

        """
        This function ...
        :param gs:
        :param nrows:
        :param ncols:
        :param projections:
        :param empty:
        :return:
        """

        from ..tools import types

        # Initialize structure to contain the plots
        plots = [[None for i in range(ncols)] for j in range(nrows)]

        # Loop over the rows
        for row in range(nrows):

            # Loop over the columns
            for col in range(ncols):

                # No plot necessary
                if empty is not None and (row, col,) in empty: continue

                # Get sub plot specification
                rect = grid[row, col]

                # Get the projection
                if projections is not None:
                    if types.is_sequence(projections) or types.is_dictionary(projections): projection = projections[row][col]
                    else: projection = projections  # same projection for each plot
                else: projection = None

                # Get the axes
                plot = self.figure.add_subplot(rect, projection=projection)

                # Create plot
                plot = MPLPlot(plot=plot)

                # Add the plot
                plots[row][col] = plot

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def set_xscale(self, scale):

        """
        This function ...
        :param scale:
        :return:
        """

        #plt.xscale(scale)
        self.ax.set_xscale(scale)

    # -----------------------------------------------------------------

    def set_yscale(self, scale):

        """
        This function ...
        :param scale:
        :return:
        """

        #plt.yscale(scale)
        self.ax.set_yscale(scale)

    # -----------------------------------------------------------------

    def set_x_log_scale(self):

        """
        This function ...
        :return:
        """

        #plt.xscale("log")
        self.set_xscale("log")

    # -----------------------------------------------------------------

    def set_y_log_scale(self):

        """
        This function ...
        :return:
        """

        #plt.yscale("log")
        self.set_yscale("log")

    # -----------------------------------------------------------------

    def set_borders(self, config):

        """
        This function ...
        :return:
        """

        # Set border width
        if config.add_border:
            #print("Adding borders ...")
            [i.set_linewidth(config.borderwidth) for i in self.ax.spines.itervalues()]

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
        ax.set_position([box.x0, box.y0 + box.height * percentage / 100., box.width, box.height * (100 - percentage) / 100.])

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

    def tight_layout(self):

        """
        This function ...
        :return:
        """

        self.figure.tight_layout()

    # -----------------------------------------------------------------

    def set_xlabel(self, label, fontsize=18):

        """
        This function ...
        :param label:
        :param fontsize:
        :return:
        """

        self.ax.set_xlabel(label, fontsize=fontsize)

    # -----------------------------------------------------------------

    def set_ylabel(self, label, fontsize=18):

        """
        This function ...
        :param label:
        :param fontsize:
        :return:
        """

        self.ax.set_ylabel(label, fontsize=fontsize)

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
        #plt.xlabel(x_label, fontsize=fontsize)
        #plt.ylabel(y_label, fontsize=fontsize)

        self.set_xlabel(x_label, fontsize=fontsize)
        self.set_ylabel(y_label, fontsize=fontsize)

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
        ticks = x_range.log(nxticks)
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

    def set_title(self, title, width=60):

        """
        This function ...
        :param title:
        :param width:
        :return:
        """

        title = title.replace("_", "\_").replace("&", "\&")
        self.figure.suptitle("\n".join(wrap(title, width)))

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

    def draw(self):

        """
        This function ...
        :return:
        """

        self.figure.canvas.draw()

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

    def show(self, block=True):

        """
        This function ...
        :param block
        :return:
        """

        # Inform the user
        log.info("Showing the plot ...")

        # Show the figure
        if block: plt.show(self.figure) #self.figure.show() #block=block)
        else: self.figure.draw()

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

    def to_file(self, path, tight=True):

        """
        This function ...
        :param path:
        :param tight:
        :return:
        """

        # Inform the user
        log.info("Saving the plot to " + str(path) + " ...")

        # Save the figure
        if tight: bbox_inches = "tight"
        else: bbox_inches = None
        #plt.savefig(path, pad_inches=0.25, transparent=self.transparent, format=self.format, bbox_inches=bbox_inches)
        self.figure.savefig(path, pad_inches=0.25, transparent=self.transparent, format=self.format, bbox_inches=bbox_inches)

    # -----------------------------------------------------------------

    def to_buffer(self, buf, tight=True):

        """
        This function ...
        :param buf:
        :param tight:
        :return:
        """

        # Inform the user
        log.info("Saving the plot to a buffer ...")

        # Save to buffer
        if tight: bbox_inches = "tight"
        else: bbox_inches = None
        #plt.savefig(buf, pad_inches=0.25, transparent=self.transparent, format=self.format, bbox_inches=bbox_inches)
        self.figure.savefig(buf, pad_inches=0.25, transparent=self.transparent, format=self.format, bbox_inches=bbox_inches)

    # -----------------------------------------------------------------

    def close(self):

        """
        This function ...
        :return:
        """

        plt.close(self.figure)

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

# JS_CODE = """
# import {Label, LabelView} from "models/annotations/label"
#
# export class LatexLabelView extends LabelView
#   render: () ->
#
#     #--- Start of copied section from ``Label.render`` implementation
#
#     ctx = @plot_view.canvas_view.ctx
#
#     # Here because AngleSpec does units tranform and label doesn't support specs
#     switch @model.angle_units
#       when "rad" then angle = -1 * @model.angle
#       when "deg" then angle = -1 * @model.angle * Math.PI/180.0
#
#     if @model.x_units == "data"
#       vx = @xscale.compute(@model.x)
#     else
#       vx = @model.x
#     sx = @canvas.vx_to_sx(vx)
#
#     if @model.y_units == "data"
#       vy = @yscale.compute(@model.y)
#     else
#       vy = @model.y
#     sy = @canvas.vy_to_sy(vy)
#
#     if @model.panel?
#       panel_offset = @_get_panel_offset()
#       sx += panel_offset.x
#       sy += panel_offset.y
#
#     #--- End of copied section from ``Label.render`` implementation
#
#     # Must render as superpositioned div (not on canvas) so that KaTex
#     # css can properly style the text
#     @_css_text(ctx, "", sx + @model.x_offset, sy - @model.y_offset, angle)
#
#     # ``katex`` is loaded into the global window at runtime
#     # katex.renderToString returns a html ``span`` element
#     katex.render(@model.text, @el, {displayMode: true})
#
# export class LatexLabel extends Label
#   type: 'LatexLabel'
#   default_view: LatexLabelView
# """
#
# from bokeh.models import Label
#
# class LatexLabel(Label):
#
#     """
#     A subclass of the Bokeh built-in `Label` that supports rendering
#     LaTex using the KaTex typesetting library.
#
#     Only the render method of LabelView is overloaded to perform the
#     text -> latex (via katex) conversion. Note: ``render_mode="canvas``
#     isn't supported and certain DOM manipulation happens in the Label
#     superclass implementation that requires explicitly setting
#     `render_mode='css'`).
#     """
#
#     __javascript__ = ["https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.6.0/katex.min.js"]
#     __css__ = ["https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.6.0/katex.min.css"]
#     __implementation__ = JS_CODE

# -----------------------------------------------------------------

def true_and_not_startswith(x, pattern):

    """
    This function ...
    :param x:
    :param pattern:
    :return:
    """

    return x and not x.startswith(pattern)

# -----------------------------------------------------------------
