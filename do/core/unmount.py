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

# Import standard modules
import subprocess

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.host import Host
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "the remote host to unmount")

# Read the command line arguments
setter = ArgumentConfigurationSetter("mount", "Mount a remote configured in PTS into a local directory")
config = setter.run(definition)

# -----------------------------------------------------------------

# Get host
host = Host(config.remote)

# PTS remotes directory
pts_remotes_path = fs.join(introspection.pts_root_dir, "remotes")
if not fs.is_directory(pts_remotes_path): fs.create_directory(pts_remotes_path)

# Create directory for remote
path = fs.join(pts_remotes_path, host.id)
if not fs.is_directory(path): fs.create_directory(path)

# If not yet mounted
if len(fs.files_in_path(path)) == 0: log.warning(host.id + " was not mounted")
else: subprocess.call(["umount", path])

# -----------------------------------------------------------------
