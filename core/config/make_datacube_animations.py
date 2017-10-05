#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
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

formats = ["avi", "gif", "apng"]
default_format = "gif"

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Instruments for which to create animations (None means all)
definition.add_positional_optional("instruments", "string_list", "instruments for which to create the images")

# -----------------------------------------------------------------

# Add optional
definition.add_optional("output", "string", "output directory")

# Write frames?
definition.add_flag("write_frames", "write the frames as seperate images in a directory", False)

# -----------------------------------------------------------------

# For PNG
definition.add_optional("colours", "string", "colour or colour map for plotting", default=default_color)
definition.add_optional("scale", "string", "scaling", default_scale, scales)
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("alpha_method", "string", "alpha method", default_alpha_method, suggestions=alpha_methods)
definition.add_optional("peak_alpha", "real", "alpha of peak value", 1.5)

# -----------------------------------------------------------------

# For animation
definition.add_optional("fps", "positive_integer", "frames per second", 10)
definition.add_optional("format", "string", "output format", default_format, choices=formats)

# -----------------------------------------------------------------
