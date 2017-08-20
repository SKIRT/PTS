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

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Method
definition.add_optional("method", "string", "method for deprojection", default_method, choices=methods)

# Writing
definition.add_section("writing", "writing options")
definition.sections["writing"].add_flag("deprojections", "write the deprojections", True)

# -----------------------------------------------------------------
