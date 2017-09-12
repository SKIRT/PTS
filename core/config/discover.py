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

# Positional optional
definition.add_positional_optional("directories", "directorypath_list", "input directories to search in")

# Flags
definition.add_flag("recursive", "look for simulations in directories recursively", True)
definition.add_flag("list", "list the found simulations", True)

# NEW
definition.add_flag("output", "list the simulation output", False)

# -----------------------------------------------------------------
