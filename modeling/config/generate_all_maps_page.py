#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.view.html import scales, colormaps, zooms

# -----------------------------------------------------------------

default_colormap = "jet"
default_scale = "log"
default_zoom = "toFit"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")
definition.add_flag("show", "show the page", False)

# Plot options
definition.add_optional("colours", "string", "colour or colour map", default=default_colormap, choices=colormaps)
definition.add_optional("scale", "string", "image scaling", default=default_scale, choices=scales)
definition.add_optional("softening_start", "real", "relative radius for softening to start (relative to truncation ellipse)", 0.75)

# Flags
definition.add_flag("replot", "replot already existing figures", True)
definition.add_flag("info", "add info about the images", True)

# -----------------------------------------------------------------
