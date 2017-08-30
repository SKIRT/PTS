#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

methods = ["skirt", "pts"]
default_method = "pts"
default_downsample_factor = 2.

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Method
definition.add_optional("method", "string", "method for deprojection", default_method, choices=methods)

# Downsample
definition.add_optional("downsample_factor", "positive_real", "downsample factor of the deprojected map w.r.t the original map", default_downsample_factor)

# Writing
definition.add_section("writing", "writing options")
definition.sections["writing"].add_flag("deprojections", "write the deprojections", True)
definition.sections["writing"].add_flag("maps", "write the maps (or don't clear them)", True)

# -----------------------------------------------------------------
