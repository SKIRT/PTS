#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.hosts Print the properties of the hosts.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.remote.host import all_host_ids, load_host
from pts.core.basics.log import setup_log
from pts.core.basics.configuration import parse_logging_arguments
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

config = parse_logging_arguments("hosts")

# -----------------------------------------------------------------

# Set log level in a special way
if config.debug: setup_log("DEBUG")
else: setup_log("ERROR")

# -----------------------------------------------------------------

# Loop over the hosts
print("")
for host_id in all_host_ids():

    print(fmt.underlined + fmt.red + host_id + fmt.reset)
    print("")

    host = load_host(host_id)
    print(host)

    #if host_id == "hpc": print(host.clusters)

    print("")

# -----------------------------------------------------------------
