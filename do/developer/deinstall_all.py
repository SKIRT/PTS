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
from pts.core.tools.logging import log
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

# Loop over the remote hosts
for host_id in find_host_ids():

    # Setup
    remote = Remote()
    if not remote.setup(host_id):
        log.warning("Remote host '" + host_id + "' is offline")
        continue

    # Run uninstaller
    uninstaller = Uninstaller()
    uninstaller.run(remote=remote)

# -----------------------------------------------------------------
