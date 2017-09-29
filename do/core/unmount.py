#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.unmount Unmount a remote configured in PTS into a local directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.mounter import RemoteMounter
from pts.core.remote.host import all_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to unmount", choices=all_host_ids())

# Read the command line arguments
config = parse_arguments("mount", definition, description="Unmount a remote mounted with PTS")

# -----------------------------------------------------------------

# Create the remote mounter
mounter = RemoteMounter()

# Unmount the remote
mounter.unmount(config.remote)

# -----------------------------------------------------------------
