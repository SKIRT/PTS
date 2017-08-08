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
from pts.core.basics.log import log
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host ID")
setter = parse_arguments("sessions", definition)

# -----------------------------------------------------------------

remote = Remote()
if not remote.setup(config.remote): raise RuntimeError("The remote host '" + config.remote + "' is not available at the moment")

screen_names = remote.screen_names()

tmux_names = remote.tmux_names()

print("SCREEN")
print(screen_names)

print("")

print("TMUX")
print(tmux_names)

# -----------------------------------------------------------------
