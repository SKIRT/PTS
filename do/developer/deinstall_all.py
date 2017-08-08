#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.deinstall_all Uninstall SKIRT/PTS on all remotes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.prep.uninstaller import Uninstaller
from pts.core.remote.host import find_host_ids
from pts.core.basics.log import log
from pts.core.remote.remote import Remote
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_positional_optional("skirt_and_or_pts", "string_list", "SKIRT and/or PTS", default=["skirt", "pts"], choices=["skirt", "pts"])

# Add flags
definition.add_flag("conda", "also remove conda installation")
definition.add_flag("qt", "also remove Qt installation")
definition.add_flag("one_attempt", "only perform one attempt at connecting to a remote")

# Get the config
config = parse_arguments("deinstall_all", definition)

# -----------------------------------------------------------------

# Loop over the remote hosts
for host_id in find_host_ids():

    # Setup
    remote = Remote()
    if not remote.setup(host_id, one_attempt=config.one_attempt):
        log.warning("Remote host '" + host_id + "' is offline")
        continue

    # Create uninstaller
    uninstaller = Uninstaller()

    # Set options
    uninstaller.config.skirt_and_or_pts = config.skirt_and_or_pts
    uninstaller.config.conda = config.conda
    uninstaller.config.qt = config.qt

    # Inform the user
    log.info("Running the uninstaller for remote host '" + host_id + "' ...")

    # Uninstall
    uninstaller.run(remote=remote)

# -----------------------------------------------------------------
