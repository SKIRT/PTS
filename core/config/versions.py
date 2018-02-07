#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids, find_hosts

# -----------------------------------------------------------------

all_host_ids = find_host_ids()
all_hosts = find_hosts()

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("host_ids", "host_list", "remote host(s)", choices=all_host_ids, default=all_hosts)
definition.add_positional_optional("clustername", "string", "cluster name (if one host is specified)")

# Add optional
definition.add_optional("not_remotes", "string_list", "skip these remote hosts", choices=find_host_ids())

# Add flag
definition.add_flag("show", "show the versions on the terminal", True)
definition.add_flag("one_attempt", "only perform one attempt at connecting to a remote")

# -----------------------------------------------------------------
