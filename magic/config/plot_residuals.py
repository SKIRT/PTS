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
definition.add_flag("absolute", "plot absolute residuals", False)

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
definition.add_flag("coordinates", "show the coordinates", True)

# -----------------------------------------------------------------
