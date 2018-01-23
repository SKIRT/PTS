#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.config.plot import definition as plot_definition
from pts.core.basics.plot import plotting_libraries, mpl, all_colormaps

# -----------------------------------------------------------------

styles = ["dark", "light"]
formats = ["pdf", "png"]
default_colormap = "magma"

# -----------------------------------------------------------------

scales = ["log", "sqrt"]
default_scale = "log"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

definition.add_optional("contains", "string_list", "only load images containing this string in their name")
definition.add_optional("not_contains", "string_list", "don't load images containing this string in their name")
definition.add_optional("exact_name", "string_list", "only load images with this exact string as their name")
definition.add_optional("exact_not_name", "string_list", "don't load images with this exact string as their name")
definition.add_optional("startswith", "string_list", "only load images whose name starts with this string")
definition.add_optional("endswith", "string_list", "only load images whose name starts with this string")

# -----------------------------------------------------------------

# Layout
definition.add_optional("max_nrows", "positive_integer", "maximum number of rows (per grid)", 10)
definition.add_optional("ngrids", "positive_integer", "number of grids", 1)

# -----------------------------------------------------------------

# Load from data
definition.add_flag("from_data", "load from data written out by the residuals image grid plotter")

# -----------------------------------------------------------------

# Sort on filters
definition.add_flag("sort_filters", "sort the frames on filter", True)

# -----------------------------------------------------------------

# Downsampling
definition.add_flag("downsample", "perform downsampling")
definition.add_optional("max_npixels", "positive_integer", "maximum number of pixels to enabled downsampling", 400)

# -----------------------------------------------------------------

# Output path
definition.add_optional("output", "directory_path", "output path")

# -----------------------------------------------------------------

# Add plotting options
definition.import_section("plot", "plotting options", plot_definition)

# Change defaults: sizes per grid panel!
definition.sections["plot"].optional["xsize"].default = 3
definition.sections["plot"].optional["ysize"].default = 3

# -----------------------------------------------------------------

# The plotting library to use
definition.add_optional("library", "string", "plotting library", mpl, plotting_libraries)

# -----------------------------------------------------------------

# Plotting options
definition.add_optional("style", "string", "plotting style", "dark", choices=styles)
definition.add_flag("transparent", "transparency", True)
definition.add_optional("format", "string", "plotting format", "pdf", choices=formats)
definition.add_optional("colormap", "string", "color map", default_colormap, choices=all_colormaps)
definition.add_optional("residuals_colormap", "string", "color map for residual frames (none means equal to image colormap)", choices=all_colormaps)

# -----------------------------------------------------------------

definition.add_optional("scale", "string", "scaling", default_scale, scales)
definition.add_optional("interval", "string", "interval", "pts")
definition.add_optional("alpha", "positive_real", "alpha of the images", 1)
definition.add_flag("background", "plot a background", True)

# -----------------------------------------------------------------

definition.add_flag("weighed", "plot weighed residuals", None)
definition.add_flag("distributions", "plot the residual distributions", False)
definition.add_flag("relative", "show relative residuals", True)
definition.add_flag("absolute", "show the residuals as absolute values", False)

# -----------------------------------------------------------------

# Extra flags
definition.add_flag("normalize", "normalize the images")
definition.add_flag("share_scale", "share the scales of the images")
definition.add_optional("scale_reference", "string", "name of the image to determine the scale for to use for the other images")
definition.add_flag("same_residuals_scale", "use the same scale for the residuals as for the observation and models")

# -----------------------------------------------------------------

definition.add_flag("write", "write out the processed frames, masks, regions, residual maps and distributions", False)
definition.add_flag("show", "show the plot (default is automatic)", None)

# -----------------------------------------------------------------

# Add coordinates?
definition.add_flag("coordinates", "show the coordinates", False)

# -----------------------------------------------------------------
