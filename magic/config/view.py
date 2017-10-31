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
from pts.magic.view.html import scales, colormaps, zooms

# -----------------------------------------------------------------

default_scale = "log"
default_colormap = "sls"
default_zoom = "toFit"

# -----------------------------------------------------------------

# Cretae configuration definition
definition = ConfigurationDefinition()

# Page
definition.add_optional("page_width", "positive_integer", "page width (in pixels)", 600)
definition.add_flag("menubar", "show menu bar", True)
definition.add_flag("colorbar", "show color bar", True)
definition.add_flag("resize", "allow resize", True)
definition.add_flag("scrolling", "allow scrolling", True)
definition.add_flag("onload", "display the masks and regions when the image is loaded", False) # doesn't work yet

# Additional features
definition.add_flag("panner", "add a panner", True)
definition.add_flag("magnifier", "add a magnifier", True)
definition.add_flag("combine_panner_and_magnifier", "combine panner and manifier next to each other", True)

# Regions settings
definition.add_flag("movable", "movable regions", True)
definition.add_flag("rotatable", "rotatable regions", True)
definition.add_flag("removable", "removable regions", True)
definition.add_flag("resizable", "resizable regions", True)

# View settings
definition.add_optional("scale", "string", "scale", default_scale, choices=scales)
definition.add_optional("colormap", "string", "color map", default_colormap, choices=colormaps)
definition.add_optional("zoom", "string", "zoom function", default_zoom, choices=zooms)

# -----------------------------------------------------------------

# Make page?
definition.add_flag("page", "generate the page", True)

# Show?
definition.add_flag("show", "show the view", True)

# -----------------------------------------------------------------

# Regions settings
definition.add_section("regions", "regions settings")

# Color of the regions
definition.sections["regions"].add_optional("color", "string", "color to render all regions (None means default of the regions)")

# Other settings
definition.sections["regions"].add_flag("changeable", "changeable", False)
definition.sections["regions"].add_flag("movable", "movable", False)
definition.sections["regions"].add_flag("rotatable", "rotatable", False)
definition.sections["regions"].add_flag("removable", "removable", False)
definition.sections["regions"].add_flag("resizable", "resizable", True)

# -----------------------------------------------------------------
