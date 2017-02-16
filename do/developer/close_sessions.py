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
from pts.core.tools.logging import log
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host ID")

setter = ArgumentConfigurationSetter("sessions")
config = setter.run(definition)

# -----------------------------------------------------------------

remote = Remote()
if not remote.setup(config.remote): raise RuntimeError("The remote host '" + config.remote + "' is not available at the moment")

screen_names = remote.screen_names()

tmux_names = remote.tmux_names()

for name in screen_names:

    # Inform the user
    log.info("Closing screen session '" + name + "' ...")

    # Close
    remote.close_screen_session(name)

for name in tmux_names:

    # Inform the user
    log.info("Closing tmux session '" + name + "' ...")

    # Close
    remote.close_tmux_session(name)

# -----------------------------------------------------------------
