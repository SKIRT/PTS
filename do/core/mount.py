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

# Import standard modules
import platform
import subprocess

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.mounter import RemoteMounter
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=find_host_ids())

# Read the command line arguments
setter = ArgumentConfigurationSetter("mount", "Mount a remote configured in PTS into the local filesystem")
config = setter.run(definition)

# -----------------------------------------------------------------

def open_file(path):
    if platform.system() == "Darwin": subprocess.Popen(["open", path])
    else: subprocess.Popen(["xdg-open", path])

# -----------------------------------------------------------------

# Create the mounter
mounter = RemoteMounter()

# Mount and get the mount path
path = mounter.mount(config.remote)

# Open the directory
open_file(path)

# -----------------------------------------------------------------
