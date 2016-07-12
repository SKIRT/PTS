#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.compare_python_packages

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.remote import Remote
from pts.core.basics.configuration import Configuration
from pts.core.tools import introspection
from pts.core.tools.logging import log

# -----------------------------------------------------------------

# Create configuration instance
config = Configuration("compare_python_packages")

# Add required
config.add_required("remote", str, "the remote host ID")

# Read
config.read()

# Settings
settings = config.get_settings()

# -----------------------------------------------------------------

# Create the remote execution environment
remote = Remote()

# Log in
remote.setup(settings.remote)

# Get all python packages installed on the remote host
remote_packages = remote.python_packages

# Get local python package version
local_versions = introspection.get_pip_versions()

# Loop over all python packages necessary for PTS
for dependency in introspection.get_all_dependencies():

    # Skip standard modules
    if introspection.is_std_lib(dependency): continue

    # Check if present
    locally_present = introspection.is_present(dependency)

    # Get version
    #version = local_versions[dependency.lower()]

    # Check if present remotely
    remotely_present = dependency in remote_packages

    if locally_present and remotely_present: log.success(dependency + ": OK")
    elif not locally_present: log.error(dependency + ": not present on this system")
    elif not remotely_present: log.error(dependency + ": not present on " + settings.remote)

# -----------------------------------------------------------------
