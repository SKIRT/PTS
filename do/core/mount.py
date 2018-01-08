#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.mount Mount a remote configured in PTS into a local directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.mounter import mount_and_open
from pts.core.remote.host import all_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=all_host_ids())
definition.add_positional_optional("path", "directory_path", "path of directory in which to create the mount point")

# Read the command line arguments
config = parse_arguments("mount", definition, description="Mount a remote configured in PTS into the local filesystem")

# -----------------------------------------------------------------

# Mount and open
mount_and_open(config.remote, path=config.path)

# -----------------------------------------------------------------
