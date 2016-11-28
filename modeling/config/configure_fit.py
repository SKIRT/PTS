#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add optional
definition.add_optional("parameters", "string_list", "parameters to be used as free parameters during the fitting")
definition.add_section("ranges", "parameter ranges")
definition.add_optional("filters", "string_list", "fit to the observed data of these filters")

# -----------------------------------------------------------------
