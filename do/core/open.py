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
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_finish
from pts.core.remote.mounter import RemoteMounter
from pts.core.remote.host import find_host_ids
from pts.core.tools import filesystem as fs
from pts.core.remote.remote import Remote
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host", choices=find_host_ids())
definition.add_required("filename", "string", "file name (or path)")

# Read the command line arguments
config = parse_arguments("open", definition, description="Open a file on a remote host")

# -----------------------------------------------------------------

# Connect to remote
remote = Remote()
remote.setup(config.remote)

# Determine the absolute file path
filepath = remote.absolute_path(config.filename)

# Determine the file path relative to the home directory
relative_filepath = remote.relative_to_home(filepath)

# Check if the file exists; otherwise don't bother mounting
if not remote.is_file(filepath):
    log.error("The file does not exist")
    exit()

# Disconnect the remote
remote.logout()

# -----------------------------------------------------------------

# Create the mounter
mounter = RemoteMounter()

# Mount and get the mount path
mount_path = mounter.mount(config.remote)

# Determine the local file path
filepath = fs.join(mount_path, relative_filepath)

# Open the file
fs.open_file(filepath)

# Wait
#prompt_finish()

# -----------------------------------------------------------------
