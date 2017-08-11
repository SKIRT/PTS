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
definition.add_flag("show", "show the page", False)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# Flags
definition.add_flag("thumbnails", "add map thumbnails", True)

definition.add_optional("thumbnail_height", "positive_integer", "height of the thumbnails (in pixels)", 50)

# -----------------------------------------------------------------
