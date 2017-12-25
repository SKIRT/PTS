#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

sigma_levels = [0.5, 0.75, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8., 9., 10., 12., 15., 20., 25., 30.]
default_sigma_level = 1.

# -----------------------------------------------------------------

default_color = "jet"
default_mask_color = "black"
default_dark_mask_color = "white"

scales = ["log", "sqrt"]
default_colour = "jet"
default_interval = "pts"

# -----------------------------------------------------------------

definition = definition.copy()

# Create the configuration
definition.add_flag("show", "show the page", False)

# Sigma levels
definition.add_positional_optional("sigma_levels", "real_list", "different sigma levels for which to generate significance masks", sigma_levels)
definition.add_optional("default_level", "real", "default sigma level", default_sigma_level)

# Flags
definition.add_flag("replot", "replot already existing figures", False) # default False because prep data doesn't really change (normally)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# Image
definition.add_optional("image_width", "positive_integer", "width of the image")
definition.add_optional("image_height", "positive_integer", "height of the image", 400)

# CROP
definition.add_optional("cropping_factor", "real", "cropping factor relative to truncation box", 1.3)

# For clip mask
definition.add_optional("min_npixels", "positive_integer", "minimum number of pixels", 1)
definition.add_optional("connectivity", "positive_integer", "connectiviy", 4)

# For PNG
definition.add_optional("colours", "string", "colour or colour map for plotting", default=default_color)
definition.add_optional("scale", "string", "scaling", "log", scales)
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("alpha_method", "string", "alpha method", "combined")
definition.add_optional("peak_alpha", "real", "alpha of peak value", 2.)

# For masks
definition.add_optional("mask_colour", "string", "colour for the mask", default=default_mask_color)
definition.add_optional("dark_mask_colour", "string", "colour for the dark mask", default=default_dark_mask_color)
definition.add_flag("mask_alpha", "use alpha for the mask", True)
definition.add_flag("fuzzy_mask", "use fuzzy masks", True)
definition.add_optional("fuzziness", "percentage", "relative fuzziness edge width", "50", convert_default=True)
definition.add_flag("dark_masks", "plot masks for dark theme", False)
definition.add_flag("remove_dark_masks", "remove the dark mask files", False)

definition.add_optional("normalization_ellipse_factor", "real", "normalize the image within the truncation ellipse scaled with this factor", 0.7)

# -----------------------------------------------------------------
