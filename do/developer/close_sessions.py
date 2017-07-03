#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.sessions Check screen/tmux sessions on the remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools.logging import setup_log
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.remote import Remote
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_optional("remotes", "string_list", "remote host ID", choices=find_host_ids(), default=find_host_ids())

setter = ArgumentConfigurationSetter("sessions")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create logger
level = "DEBUG" if config.debug else "INFO"
log = setup_log(level)

# -----------------------------------------------------------------

# Loop over the remote host
for host_id in config.remotes:

    # Create and setup the remote
    remote = Remote()
    if not remote.setup(host_id): raise RuntimeError("The remote host '" + host_id + "' is not available at the moment")

    # Get screen names and tmux session names
    screen_names = remote.screen_names()
    tmux_names = remote.tmux_names()

    # Inform the user
    log.info("Closing sessions on remote host '" + host_id + "' ...")

    # Loop over the screens
    for name in screen_names:

        # Inform the user
        log.info("Closing screen session '" + name + "' ...")

        # Close
        remote.close_screen_session(name)

    # Loop over the tmux sessions
    for name in tmux_names:

        # Inform the user
        log.info("Closing tmux session '" + name + "' ...")

        # Close
        remote.close_tmux_session(name)

# -----------------------------------------------------------------
