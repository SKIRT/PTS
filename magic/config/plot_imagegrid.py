# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.plot import normal_colormaps
from pts.magic.plot.imagegrid import default_cmap, light_theme, themes

# ------------------------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# From directory
definition.add_optional("from_directory", "directory_path", "load images from directory")

# Colormap
definition.add_optional("cmap", "string", "colormap", default_cmap, choices=normal_colormaps)

# Options
definition.add_optional("axes_label_size", "positive_integer", "axes label size", 14)
definition.add_optional("ticks_label_size", "positive_integer", "ticks label size", 8)
definition.add_optional("legend_fontsize", "positive_integer", "legend fontsize", 14)
definition.add_optional("legend_markers_cale", "positive_integer", "legend marker scale", 0)
definition.add_optional("lines_marker_size", "positive_real", "lines marker size", 2.5)
definition.add_optional("linewidth", "positive_integer", "linewidth", 1)

# Dark or light theme?
definition.add_optional("theme", "string", "theme for the plot", light_theme, choices=themes)

# Pixelscale
#definition.add_optional("pixelscale", "angle", "pixelscale for the images", "5 arcsec", convert_default=True)
definition.add_optional("radius", "angle", "radius of the field of view", "0.1 deg", convert_default=True)

# Interval
definition.add_optional("interval", "string", "interval", "pts")
definition.add_optional("scale", "string", "scale", "log")

# Path for the plot file
definition.add_optional("path", "string", "output filepath")

# Fontsizes
definition.add_optional("label_fontsize", "positive_real", "fontsize for image labels", 18)

# Colors
definition.add_optional("background_color", "string", "background color in plots (None means colormap background)")
definition.add_optional("text_color_in", "string", "color of all text in plots", "white")

# Center
definition.add_optional("center", "skycoordinate", "center coordinate")

# Spacing
definition.add_optional("spacing", "angle", "spacing of the ticks", "1 arcsec", convert_default=True)

# For saving plot file
definition.add_optional("dpi", "positive_integer", "dots per inch", 300)

# ------------------------------------------------------------------------------

# Vmin and vmax?
definition.add_optional("vmin", "real", "vmin")
definition.add_optional("vmax", "real", "vmax")

# ------------------------------------------------------------------------------

# Showing
definition.add_flag("show", "show", True)

# ------------------------------------------------------------------------------

# Writing
definition.add_flag("write", "do writing", True)

# ------------------------------------------------------------------------------

definition.add_optional("major_tick_length", "positive_integer", "length of major ticks", 7)
definition.add_optional("minor_tick_length", "positive_integer", "length of minor ticks", 4)

# ------------------------------------------------------------------------------
