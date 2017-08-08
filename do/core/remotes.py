#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.remotes Check the status of the remotes (whether they are up).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.tools import formatting as fmt
from pts.core.basics.log import setup_log
from pts.core.basics.configuration import parse_logging_arguments

# -----------------------------------------------------------------

config = parse_logging_arguments("hosts")

# -----------------------------------------------------------------

# Set log level in a special way
if config.debug: setup_log("DEBUG")
else: setup_log("ERROR")

# -----------------------------------------------------------------

# Loop over the hosts
for host_id in find_host_ids():

    remote = Remote()
    connected = remote.setup(host_id, login_timeout=15) # only try for 15 seconds
    if connected: print(fmt.green + host_id + ": up" + fmt.reset)
    else: print(fmt.red + host_id + ": down" + fmt.reset)

# -----------------------------------------------------------------
