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
from pts.core.remote.remote import Remote
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import introspection
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "the remote host ID")

# Add flags
definition.add_flag("versions", "compare versions", "v")

# Get configuration
config = parse_arguments("compare_python_packages", definition)

# -----------------------------------------------------------------

# Create the remote execution environment
remote = Remote()

# Log in
remote.setup(config.remote)

# Get remote python session
remote_python = remote.start_python_session(assume_pts=False)

# Get all python packages installed on the remote host
remote_packages = remote_python.installed_packages

# Get local python package version
local_packages = introspection.installed_python_packages()

# Loop over all python packages necessary for PTS
for dependency in introspection.get_all_dependencies():

    # Skip standard modules
    if introspection.is_std_lib(dependency): continue

    # Check present and version locally
    if dependency in local_packages:
        locally_present = True
        local_version = local_packages[dependency]
    else:
        # Check again for present by importing
        locally_present = introspection.is_present_package(dependency)
        local_version = None

    # Check present and version remotely
    if dependency in remote_packages:
        remotely_present = True
        remote_version = remote_packages[dependency]
    else:
        # Check again for present by importing
        remotely_present = remote_python.is_present_package(dependency)
        remote_version = None

    # If present both locally and remotely
    if locally_present and remotely_present:

        if config.versions:

            if local_version is None and remote_version is not None: log.warning(dependency + ": local version unknown")
            elif remote_version is None and local_version is not None: log.warning(dependency + ": remote version unknown")
            elif remote_version is None and local_version is None: log.warning(dependency + ": local and remote version unknown")
            elif local_version == remote_version: log.success(dependency + ": OK")
            else: log.warning(dependency + ": version " + local_version + " locally and version " + remote_version + " remotely")

        else: log.success(dependency + ": OK")

    # Not present on at least one system
    elif remotely_present and not locally_present: log.error(dependency + ": not present on this system")
    elif locally_present and not remotely_present: log.error(dependency + ": not present on remote '" + config.remote + "'")
    else: log.error(dependency + ": not present on either this sytem or remote '" + config.remote + "'")

# -----------------------------------------------------------------
