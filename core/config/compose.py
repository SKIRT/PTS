#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Ski file path
definition.add_positional_optional("skifile", "file_path", "ski file path")

# -----------------------------------------------------------------

# Commands to be run
definition.add_optional("commands", "string_list", "commands to be run in interactive mode")

# -----------------------------------------------------------------
