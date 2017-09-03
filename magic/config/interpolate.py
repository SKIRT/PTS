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
from pts.magic.core.cutout import interpolation_methods

# -----------------------------------------------------------------

default_method = "pts"
default_interpolation_method = "pts"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_required("image", "string", "name of the input image")
definition.add_required("regions", "string", "name of the regions file over which to interpolate")
definition.add_flag("write_mask", "write out the mask")
definition.add_optional("color", "string", "only interpolate over the shapes with this color")
definition.add_optional("ignore_color", "string", "ignore shapes with this particular color")
definition.add_optional("shapes", "string_list", "only interpolate over these kinds of shapes")
definition.add_optional("input", "string", "name of the input directory", letter="i")
definition.add_optional("output", "string", "name of the output directory", letter="o")
definition.add_optional("method", "string", "interpolation method to use", default=default_method)
definition.add_optional("interpolation_method", "string", "interpolation method", default_interpolation_method, choices=interpolation_methods)
definition.add_flag("sigma_clip", "apply sigma clipping before interpolation", True)
definition.add_optional("source_outer_factor", "real", "outer factor", 1.4)
definition.add_flag("plot", "plot after interpolation", False)
definition.add_flag("replace", "allow the original image to be replaced", False)
definition.add_flag("backup", "backup if replaced", True)
definition.add_optional("backup_suffix", "string", "backup suffix", "backup")

# -----------------------------------------------------------------
