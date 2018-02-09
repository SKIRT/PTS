#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.send Send a file or directory to a remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Required
definition.add_required("local_path", "string", "path or name of the file or directory to send")
definition.add_required("remote", "string", "the remote host to send to", choices=find_host_ids())

# Remote path
definition.add_positional_optional("remote_path", "string", "path of the remote directory to send to")

# Create configuration
config = parse_arguments("send", definition, "Send a file or directory to a remote host")

# -----------------------------------------------------------------

# Create remote
remote = Remote(host_id=config.remote)

# -----------------------------------------------------------------

# Set full path of origin
origin = fs.absolute_or_in_cwd(config.local_path)
name = fs.name(origin)

# Set full path to the destination
if config.remote_path is None: destination = fs.join(remote.home_directory, name)
else: destination = remote.absolute_path(config.remote_path)

# Upload
remote.upload(origin, destination)

# -----------------------------------------------------------------
