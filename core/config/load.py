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

# Remote hosts
definition.add_positional_optional("hosts", "host_list", "remote hosts", choices=all_host_ids, default=all_hosts)

# Add flags
definition.add_flag("show", "show the loads", True)
definition.add_flag("one_attempt", "only perform one attempt at connecting to a remote")

# -----------------------------------------------------------------
