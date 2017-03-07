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
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools.logging import setup_log
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_positional_optional("remotes", "string_list", "remote hosts on which to clear", choices=find_host_ids(), default=find_host_ids())

# Create setter
setter = ArgumentConfigurationSetter("clear_sessions_and_temp")
config = setter.run(definition)

# -----------------------------------------------------------------

# Set logger
level = "DEBUG" if config.debug else "INFO"
log = setup_log(level)

# -----------------------------------------------------------------

# Loop over the remote hosts
for host_id in config.remotes:

    # Setup the remote
    remote = Remote()
    if not remote.setup(host_id): log.warning("Could not connect to remote host '" + host_id + "'")

    # Clear temporary directory
    remote.clear_pts_temp()

    # Clear sessions
    remote.close_all_screen_sessions()
    remote.close_all_tmux_sessions()

# -----------------------------------------------------------------
