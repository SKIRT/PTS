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
import matplotlib.gridspec as gridspec
from matplotlib.ticker import LinearLocator, LogLocator, AutoMinorLocator, AutoLocator, NullLocator
from matplotlib.ticker import ScalarFormatter, NullFormatter, LogFormatter, PercentFormatter, EngFormatter, LogFormatterMathtext, LogFormatterSciNotation

# Import the relevant PTS classes and modules
from ..basics.log import log

# -----------------------------------------------------------------

class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self, vmin, vmax):  # Override function that finds format to use.
        #print(vmin, vmax)
        self.format = "%1.1f"  # Give format here

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

    #

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

        """
        This function ...
        :return:
        """

        return self._plot.axes

    # -----------------------------------------------------------------

    @property
    def xaxis(self):

        """
        Thisfunction ...
        :return:
        """

        return self.axes.get_xaxis()

    # -----------------------------------------------------------------

    @property
    def yaxis(self):

        """
        This function ...
        :return:
        """

        return self.axes.get_yaxis()

    # -----------------------------------------------------------------

    @property
    def xticklabels(self):

        """
        This function ...
        :return:
        """

        return self._plot.get_xticklabels()

    # -----------------------------------------------------------------

    def set_xticks(self, ticks=None, fontsize=None, major_locator=None, minor_locator=None, major_formatter=None,
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
            self._plot.set_xticks(ticks)
            self._plot.set_xticklabels(ticks)

        # Set locator automatically
        if major_locator is None:
            if self.log_xscale: major_locator = LogLocator(subs=log_subs)
            elif self.linear_xscale: major_locator = LinearLocator()
            else: raise ValueError("Unknown xscale")
        if minor_locator is None:
            if minor: # have minor ticks
                if self.log_xscale: minor_locator = None
                elif self.linear_xscale: minor_locator = AutoMinorLocator()
                else: raise ValueError("Unknown xscale")
            else: minor_locator = NullLocator()

        # Set the locators
        self.xaxis.set_major_locator(major_locator)
        if minor_locator is not None: self.xaxis.set_minor_locator(minor_locator)

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

        # Set fontsize
        plt.setp(self.xticklabels, rotation='horizontal', fontsize=fontsize)

    # -----------------------------------------------------------------

    @property
    def yticklabels(self):

        """
        Thisn function ...
        :return:
        """

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
        plt.setp(self.yticklabels, rotation='horizontal', fontsize=fontsize)

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

        """
        This function ...
        :return:
        """

        return self.xaxis.get_scale()

    # -----------------------------------------------------------------

    @property
    def linear_xscale(self):

        """
        This function ...
        :return:
        """

        return self.xscale == "linear"

    # -----------------------------------------------------------------

    @property
    def log_xscale(self):

        """
        This function ...
        :return:
        """

        return self.xscale == "log"

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

    @property
    def yscale(self):

        """
        This function ...
        :return:
        """

        return self.yaxis.get_scale()

    # -----------------------------------------------------------------

    @property
    def linear_yscale(self):

        """
        This function ...
        :return:
        """

        return self.yscale == "linear"

    # -----------------------------------------------------------------

    @property
    def log_yscale(self):

        """
        This function ...
        :return:
        """

        return self.yscale == "log"

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
        self.transparent = False
        self.format = None

    # -----------------------------------------------------------------

    def create_column(self, size, share_axis=False, height_ratios=None, x_label=None, x_label_fontsize="small",
                      x_labels=None, y_labels=None, y_label_fontsize="small", x_scale="linear", x_scales=None, y_scales=None,
                      x_limits=None, y_limits=None, x_log_scalar=False, y_log_scalar=False):

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
        :return:
        """

        # Define hspace
        wspace = 0.0
        if share_axis: hspace = 0.0
        else: hspace = 0.05

        # Make the grid
        gs = gridspec.GridSpec(size, 1, height_ratios=height_ratios, hspace=hspace)

        # Create the (sub)plots
        plots = []
        if share_axis:

            if x_labels is not None: raise ValueError("Cannot specify different x labels when sharing axis")
            if x_scales is not None: raise ValueError("Cannot specify different x scales when sharing axis")

            first_mpl_plot = plt.subplot(gs[0])
            first_plot = MPLPlot(plot=first_mpl_plot)

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

                next_plot = plt.subplot(gs[index], sharex=first_mpl_plot)
                next_plot = MPLPlot(plot=next_plot)

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

                plots.append(next_plot)

            last_plot = plots[-1]

            # Set shared x label
            if x_label is not None:
                last_plot.set_xlabel(x_label, fontsize=x_label_fontsize)

            # Set shared x scale
            #last_plot.set_xscale(x_scale)

        else:

            if x_label is not None: raise ValueError("Cannot specify one x label when not sharing axis")
            #if x_scale is not None:

            # Create plots
            for index in range(size):

                # Get plot
                plot = plt.subplot(gs[index])

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
                plots.append(MPLPlot(plot=plot))

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_row(self, size, share_axis=False, width_ratios=None, y_label=None, y_label_fontsize="small",
                   y_labels=None, x_labels=None, x_label_fontsize="small", y_scale="linear", y_scales=None,
                   x_scales=None, y_limits=None, x_limits=None, y_log_scalar=False, x_log_scalar=False):

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
        :return:
        """

        # Define hspace
        if share_axis: wspace = 0.0
        else: wspace = 0.05

        # Make the grid
        gs = gridspec.GridSpec(1, size, width_ratios=width_ratios, wspace=wspace)

        # Create the (sub)plots
        plots = []
        if share_axis:

            if y_labels is not None: raise ValueError("Cannot specify different y labels when sharing axis")
            if y_scales is not None: raise ValueError("Cannot specify different y scales when sharing axis")

            first_mpl_plot = plt.subplot(gs[0])
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

                next_plot = plt.subplot(gs[index], sharey=first_mpl_plot)
                next_plot = MPLPlot(plot=next_plot)

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

                # Get plot
                plot = plt.subplot(gs[index])

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
                plots.append(MPLPlot(plot=plot))

        # Return the plots
        return plots

    # -----------------------------------------------------------------

    def create_one_plot(self):

        """
        Thisf unction ...
        :return:
        """

        plot = MPLPlot(plot=self.ax)
        return plot

    # -----------------------------------------------------------------

    def set_xscale(self, scale):

        """
        This function ...
        :param scale:
        :return:
        """

        plt.xscale(scale)

    # -----------------------------------------------------------------

    def set_yscale(self, scale):

        """
        This function ...
        :param scale:
        :return:
        """

        plt.yscale(scale)

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

        title = title.replace("_", "\_")

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
