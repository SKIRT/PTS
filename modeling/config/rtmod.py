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

# Commands to be run
definition.add_positional_optional("commands", "string_list", "commands to be run in interactive mode")

# Interactive mode
definition.add_flag("interactive", "use interactive mode", True)

# -----------------------------------------------------------------

# Caching
definition.add_optional("cache_volume", "string", "name of the volume to be used for caching")

# -----------------------------------------------------------------

# Set plotting backend
definition.add_optional("mpl_backend", "string", "set the matplotlib backend")

# -----------------------------------------------------------------
