#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from ...core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add settings
definition.add_positional_optional("database", "file_path", "path to the database file")
definition.add_positional_optional("runs", "string_list", "plot these runs from the database (None = all)")

# Output
definition.add_optional("output", "directory_path", "output path", letter="o")

# Flags
definition.add_optional("minmax", "string", "whether the scores should be minimized or maximized", "max")

definition.add_optional("format", "string", "file format", "pdf")

# -----------------------------------------------------------------
