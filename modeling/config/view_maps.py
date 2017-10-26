#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

types = ["colours", "ssfr", "tir", "attenuation", "old", "dust", "young", "ionizing", "components"]

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Which maps?
definition.add_required("which", "string_list", "types of maps to be plotted", choices=types)

# -----------------------------------------------------------------
