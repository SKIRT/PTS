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
from pts.core.remote.host import find_host_ids
from pts.core.prep.update import SKIRTUpdater, PTSUpdater

# -----------------------------------------------------------------

for host_id in find_host_ids():

    # Create and configure the SKIRT updater
    skirt = SKIRTUpdater()
    skirt.config.remote = host_id

    # Update SKIRT
    skirt.run()

    # Create and configure the PTS updater
    pts = PTSUpdater()
    pts.config.remote = host_id

    # Update PTS
    pts.run()

# -----------------------------------------------------------------
