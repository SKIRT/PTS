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
definition.add_positional_optional("host_ids", "string_list", "remote host IDs", choices=find_host_ids(), default=find_host_ids())
definition.add_optional("cluster_name", "string", "cluster name")

# Add flags
definition.add_flag("show", "show the loads", True)

# -----------------------------------------------------------------
