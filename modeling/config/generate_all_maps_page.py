#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

default_cmap = "jet"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")
definition.add_flag("show", "show the page", False)

# Colour map
definition.add_optional("colours", "string", "colour or colour map", default=default_cmap)

# Flags
definition.add_flag("replot", "replot already existing figures", True)
definition.add_flag("info", "add info about the images", True)

# -----------------------------------------------------------------
