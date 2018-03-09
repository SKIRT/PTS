#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.synchronize Synchronize a remote directory with a local one.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.remote import Remote
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Required
definition.add_required("local_path", "directory_path", "path or name of the local directory")
definition.add_required("remote", "string", "the remote host to send to", choices=find_host_ids())
definition.add_required("remote_path", "string", "path of remote directory to synchronize")
definition.add_flag("create", "create remote directory if necessary", False)

# Create the configuration
config = parse_arguments("show_datacubes", definition)

# -----------------------------------------------------------------

# Create remote
remote = Remote(host_id=config.remote)

# -----------------------------------------------------------------

# Create?
remote_path = remote.absolute_path(config.remote_path)
if not remote.is_directory(remote_path):
    if config.create: remote.create_directory(remote_path, recursive=True)
    else: raise IOError("Remote directory '" + remote_path + "' does not exist")

# -----------------------------------------------------------------

# Synchronize
remote.synchronize(config.local_path, remote_path, show_output=True)

# -----------------------------------------------------------------
