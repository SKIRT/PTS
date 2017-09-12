#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.clear_sessions_and_temp Clear the complete PTS temporary directory and sessions on remote hosts.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_positional_optional("remotes", "string_list", "remote hosts on which to clear", choices=find_host_ids(), default=find_host_ids())

# Create setter
config = parse_arguments("clear_sessions_and_temp", definition)

# -----------------------------------------------------------------

# Loop over the remote hosts
for host_id in config.remotes:

    # Setup the remote
    remote = Remote()
    if not remote.setup(host_id):
        log.warning("Could not connect to remote host '" + host_id + "'")
        continue

    # Clear temporary directory and clear sessions
    remote.clear_temp_and_sessions()

# -----------------------------------------------------------------
