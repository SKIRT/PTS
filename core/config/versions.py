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
definition.add_positional_optional("host_ids", "string_list", "remote host(s)", choices=find_host_ids(), default=find_host_ids())

# Add optional
definition.add_optional("not_remotes", "string_list", "skip these remote hosts", choices=find_host_ids())

# Add flag
definition.add_flag("show", "show the versions on the terminal", True)

# -----------------------------------------------------------------
