#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_positional_optional("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_positional_optional("matching", "string", "only adapt settings with a name matching this string", suggestions=["remote"])
definition.add_positional_optional("ids", "integer_list", "simulation IDs")
definition.add_optional("names", "string_list", "simulation names")
definition.add_flag("from_directories", "use directory names as simulation names")

# -----------------------------------------------------------------

# Select certain properties
definition.add_optional("contains", "string", "only adapt properties containing this string in their name")
definition.add_optional("not_contains", "string", "don't adapt properties containing this string in their name")
definition.add_optional("exact_name", "string", "only adapt properties with this exact string as their name")
definition.add_optional("exact_not_name", "string", "don't adapt properties with this exact string as their name")
definition.add_optional("startswith", "string", "only adapt properties whose name starts with this string")
definition.add_optional("endswith", "string", "only adapt properties whose name starts with this string")

# -----------------------------------------------------------------

definition.add_flag("update", "update the analysis options", True)

# -----------------------------------------------------------------
