#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# The path to the hosts configuration directory
hosts_config_path = fs.join(introspection.pts_config_dir("core"), "hosts")

# The names of the hosts which are pre-configured
preconfigured_names = fs.files_in_path(hosts_config_path, extension="cfg", not_contains="template", returns="name")

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition(write_config=False)

# Add required
definition.add_required("name", "nonempty_string_no_spaces", "name to give to the host")

# Add optional
definition.add_positional_optional("preconfigured", "string", "name of a preconfigured remote for which to adapt the user-specific settings", choices=preconfigured_names)

# -----------------------------------------------------------------
