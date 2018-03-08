#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add options
definition.add_optional("xsize", "integer", "figure x size", 10)
definition.add_optional("ysize", "integer", "figure y size", 6)
definition.add_optional("main_relsize", "percentage", "main panel relative size (wrt. x or y size)", "80", convert_default=True)
definition.add_optional("linewidth", "real", "line width", default=2.0)
definition.add_optional("borderwidth", "real", "border width", default=1.0)
definition.add_flag("add_titles", "add titles to the plots", True)
definition.add_optional("legend_fontsize", "positive_integer", "fontsize of legend text", default=10)
definition.add_optional("legend_title_fontsize", "positive_integer", "fontsize of legend titles", default=14)
definition.add_flag("legend_below", "place legends below the plots")
definition.add_optional("legend_borderwidth", "real", "border width of legend", default=2.0)
definition.add_optional("legend_bordercolor", "string", "border color of legend", default="k")
definition.add_flag("add_grid", "add a grid", True)
definition.add_optional("grid_linewidth", "real", "linewidth of grid", default=1.0)
definition.add_optional("grid_linestyle", "string", "linestyle of grid", default="dotted", choices=["solid", "dashed", "dashdot", "dotted", "offset", "on-off-dash-seq", '-', '--', '-.', ':', 'None', ' ', ''])
definition.add_optional("grid_color", "string", "color of the grid lines", default="0.70")
definition.add_flag("add_border", "add border to the plot", True)
definition.add_optional("label_fontsize", "positive_integer", "fontsize of the axes labels", default=18)
definition.add_flag("add_legend_border", "add border to legend", False)
definition.add_flag("transparent_background", "transparent background", True)
definition.add_optional("title_fontsize", "positive_integer", "fontsize of the plot title", default=20)
definition.add_optional("ticks_fontsize", "positive_integer", "fontsize of the axes ticks", default=12)
definition.add_optional("markersize", "positive_integer", "size of the data point markers", default=12)
definition.add_flag("connect_points", "connect data points with lines", True)

# -----------------------------------------------------------------
