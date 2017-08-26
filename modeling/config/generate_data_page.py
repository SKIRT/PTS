#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.generate_page import definition
from pts.modeling.core.environment import verify_modeling_cwd
from pts.magic.core.rgba import alpha_methods

# -----------------------------------------------------------------

default_alpha_method = "combined"
default_scale = "log"
default_color = "jet"
default_mask_color = "black"

scales = ["log", "sqrt"]
default_colour = "jet"
default_interval = "pts"

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Flags
definition.add_flag("replot", "replot", False)
definition.add_flag("use_session", "use remote python session to create the images", False)

# Group
definition.add_flag("group_observatories", "group the images with observatories", False)

# Flags
definition.add_flag("thumbnails", "add map thumbnails", True)
definition.add_optional("thumbnail_height", "positive_integer", "height of the thumbnails (in pixels)", 50)
definition.add_flag("previews", "add previews of the maps when hovering over the thumbnails", True)

# Maxium numberof pixels
definition.add_optional("max_npixels", "positive_integer", "maximum number of pixels for the plots", 1000)

# Flags
definition.add_flag("fetch_missing", "fetch missing images from the DustPedia archive", False)

# For PNG
definition.add_optional("colours", "string", "colour or colour map for plotting", default=default_color)
definition.add_optional("scale", "string", "scaling", default_scale, scales)
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("alpha_method", "string", "alpha method", default_alpha_method, suggestions=alpha_methods)
definition.add_optional("peak_alpha", "real", "alpha of peak value", 2.)

# -----------------------------------------------------------------
