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
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.mounter import RemoteMounter
from pts.core.remote.host import find_host_ids
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=find_host_ids())
definition.add_positional_optional("path", "directory_path", "path of directory in which to create the mount point")

# Read the command line arguments
setter = ArgumentConfigurationSetter("mount", "Mount a remote configured in PTS into the local filesystem")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create the mounter
mounter = RemoteMounter()

# Mount and get the mount path
path = mounter.mount(config.remote, config.path)

# Open the directory
fs.open_directory(path)

# -----------------------------------------------------------------
