#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Required
definition.add_positional_optional("image", "file_path", "name/path of the input image")

# Input and output
definition.add_optional("input", "directory_path", "input directory path", letter="i")
definition.add_optional("output", "directory_path", "output directory path", letter="o")

# Regions
definition.add_optional("special_region", "file_path", "region indicating areas that require special attention")
definition.add_optional("ignore_region", "file_path", "region indicating areas that should be ignored")
definition.add_optional("bad", "file_path", "region specifying areas that have to be added to the mask of bad pixels")

definition.add_flag("animation", "make an animation of the extraction procedure")

# Detailed settings
definition.add_optional("interpolation_method", "string", "interpolation method", "pts")
definition.add_flag("sigma_clip", "perform sigma-clipping when interpolating", True)
definition.add_optional("source_outer_factor", "real", "outer factor", 1.4)

definition.add_flag("dilate_saturation", "dilate saturation")
definition.add_optional("saturation_dilation_factor", "real", "saturation dilation factor", 2.0)

definition.add_flag("dilate_other", "dilate other sources")
definition.add_optional("other_dilation_factor", "real", "dilation factor for other sources", 2.0)

definition.add_flag("only_foreground", "only interpolate over the stars that are in the foreground of the galaxy", False)

definition.add_flag("write", "do writing", True)

# Flags
definition.add_flag("remove_companions", "remove companion galaxies", False)

# -----------------------------------------------------------------
