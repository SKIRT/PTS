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

# Flags
definition.add_flag("recursive", "look for simulation in directories recursively", True)

# Output
definition.add_optional("output", "directory_path", "output directory", letter="o")

# Flags
definition.add_flag("other", "also plot the 'other' phases")

# -----------------------------------------------------------------
