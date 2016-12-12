#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.update_all Update SKIRT and PTS on all remotes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids, HostDownException
from pts.core.prep.update import SKIRTUpdater, PTSUpdater
from pts.core.tools.logging import log

# -----------------------------------------------------------------

# Loop over the different hosts
for host_id in find_host_ids():

    # Inform the user
    log.info("Checking and updating SKIRT and PTS on remote host '" + host_id + "' ...")

    # Create and configure the SKIRT updater
    skirt = SKIRTUpdater()
    skirt.config.remote = host_id

    # Inform the user
    log.info("Updating SKIRT ...")

    # Update SKIRT
    try: skirt.run()
    except HostDownException:
        log.warning("Could not connect to the remote host, skipping ...")
        continue

    # Inform the user
    log.info("Updating PTS ...")

    # Create and configure the PTS updater
    pts = PTSUpdater()
    pts.config.remote = host_id

    # Update PTS
    pts.run()

# -----------------------------------------------------------------
