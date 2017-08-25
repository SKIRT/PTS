#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.magic.view.html import scales, colormaps, zooms
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

default_colour = "jet"
default_colormap = "viridis"
default_scale = "log"
default_zoom = "toFit"

# -----------------------------------------------------------------

# Create the configuration
definition.add_flag("show", "show the page", False)

# VIEW
definition.add_optional("colormap", "string", "color map", default_colormap, choices=colormaps)
definition.add_optional("scale", "string", "image scaling", default=default_scale, choices=scales)
definition.add_optional("zoom", "string", "zoom function", default_zoom, choices=zooms)

# Scales for different kinds of maps
definition.add_optional("colour_scale", "string", "scaling for colour maps", "histeq")
definition.add_optional("ssfr_scale", "string", "scaling for ssfr maps", "log")
definition.add_optional("tir_scale", "string", "scaling for TIR maps", "log")
definition.add_optional("attenuation_scale", "string", "scaling for attenuation maps", "log")
definition.add_optional("old_scale", "string", "scaling for old stellar maps", "log")
definition.add_optional("young_scale", "string", "scaling for young stellar maps", "log")
definition.add_optional("ionizing_scale", "string", "scaling for ionizing stellar maps", "log")
definition.add_optional("dust_scale", "string", "scaling for dust maps", "log")

# For PNG
definition.add_optional("colours", "string", "colour or colour map for plotting", default=default_colour)

definition.add_optional("softening_start", "real", "relative radius for softening to start (relative to truncation ellipse)", 0.75)
definition.add_flag("view_png", "use the pngs for viewing instead of the original data", False)
definition.add_optional("cropping_factor", "positive_real", "multiply the cropping box with this factor", 1.3)

# Exclusively for the views
definition.add_flag("menubar", "add menubars", True)
definition.add_flag("colorbar", "add colorbars", True)

# Flags
definition.add_flag("replot", "replot already existing figures", True)
definition.add_flag("info", "add info about the images", True)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)
definition.add_flag("replace_nans", "replace NaNs (and Infs) by zero (Infs are converted to NaN by default)", False)

# -----------------------------------------------------------------
