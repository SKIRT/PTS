#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.open Open a file on a remote host to view and/or edit it.

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
definition.add_required("remote", "string", "remote host", choices=find_host_ids())
definition.add_required("filename", "string", "file name")

# Read the command line arguments
setter = ArgumentConfigurationSetter("open", "Open a file on a remote host")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create the mounter
mounter = RemoteMounter()

# Mount and get the mount path
path = mounter.mount(config.remote)

# Determine the file path
filepath = fs.join(path, config.filename)

# Open the file
fs.open_file(filepath)

# -----------------------------------------------------------------
