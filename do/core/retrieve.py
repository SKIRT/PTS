#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.retrieve Retrieve a file or directory from a remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Required
definition.add_required("remote_path", "string", "remote path of the file or directory to retrieve")
definition.add_required("remote", "string", "remote host to retrieve from", choices=find_host_ids())

# Local path
definition.add_positional_optional("local_path", "string", "path of the local directory to store the file/directory")

# Create configuration
config = parse_arguments("retrieve", definition, "Retrieve a file or directory from a remote host")

# -----------------------------------------------------------------

# Create remote
remote = Remote(host_id=config.remote)

# -----------------------------------------------------------------

# Set full path of origin
origin = remote.absolute_path(config.remote_path)

# Set full path to the destination
if config.local_path is not None: destination = fs.absolute_or_in_cwd(config.local_path)
else: destination = fs.cwd()

# -----------------------------------------------------------------

# Debugging
log.debug("Origin: " + origin)
log.debug("Destination: " + destination)

# Upload
remote.download(origin, destination)

# -----------------------------------------------------------------
