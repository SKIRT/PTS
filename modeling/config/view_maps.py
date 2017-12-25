#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.magic.view.html import scales, colormaps, zooms
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

types = ["colours", "ssfr", "tir", "attenuation", "old", "dust", "young", "ionizing", "components"]

# -----------------------------------------------------------------

definition = definition.copy()

# Which maps?
definition.add_required("which", "string", "types of map to be plotted", choices=types)

# Methods
definition.add_positional_optional("method", "string", "map making method")
definition.add_positional_optional("startswith", "string", "map filename should start with this string")
definition.add_optional("factors", "real_list", "factors")

# Add the truncation ellipse?
definition.add_flag("truncation_ellipse", "add truncation ellipses", True)

# -----------------------------------------------------------------

default_zoom = "toFit;x2"

# -----------------------------------------------------------------

definition.add_optional("scale", "string", "scale", choices=scales) # default is auto
definition.add_optional("colormap", "string", "color map", choices=colormaps) # default is auto
definition.add_optional("zoom", "string", "zoom function", default_zoom, choices=zooms)

# -----------------------------------------------------------------

# ADVANCED
definition.add_flag("preload_all", "preload all images", False)
definition.add_optional("preload", "string_list", "names for which to preload the image")
definition.add_flag("dynamic", "create the viewers dynamically", False)

# -----------------------------------------------------------------
