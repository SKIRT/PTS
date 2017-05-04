#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

scales = ["linear", "logarithmic"]

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add optional
definition.add_optional("ngenerations", "positive_integer", "number of generations to run in one run (ngenerations > 1 is only allowed for local execution)", 1)
definition.add_optional("nmodels", "positive_integer", "number of models per generation", 80)

definition.add_flag("manual_initial_generation", "generate intiial popultion manually", False)
definition.add_optional("manual_initial_generation_scale", "string", "gegegeg", default="logarithmic", choices=scales)

# -----------------------------------------------------------------
