#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

default_scale = "logarithmic"
scales = ["linear", "logarithmic"]

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The number of models
definition.add_optional("nmodels", "positive_integer", "number of models per generation", 80)

# Scales
definition.add_optional("scales", "string_string_dictionary", "scales (linear/logarithmic) for the different free parameters")

# Make animations?
definition.add_flag("animate", "make animations", False)

# -----------------------------------------------------------------
