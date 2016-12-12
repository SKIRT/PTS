#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.versions Compare versions of PTS and SKIRT between local and remote installations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.host import find_host_ids
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("remote", "string", "remote host", choices=find_host_ids())

# Add optional
definition.add_optional("not_remotes", "string_list", "skip these remote hosts", choices=find_host_ids())

# Read the command line arguments
setter = ArgumentConfigurationSetter("open", "Open a file on a remote host")
config = setter.run(definition)

# -----------------------------------------------------------------

if config.remote is not None: host_ids = [config.remote]
else: host_ids = find_host_ids()

# -----------------------------------------------------------------

if config.not_remotes is not None:
    host_ids = [host_id for host_id in host_ids if host_id not in config.not_remotes]

skirt_versions = dict()
pts_versions = dict()

qt_versions = dict()
cpp_versions = dict()
mpi_versions = dict()

# Loop over the different hosts
for host_id in host_ids:

    # Setup the remote (login)
    remote = Remote()
    remote.setup(host_id)

    # Get SKIRT and PTS version
    skirt_version = remote.skirt_version
    pts_version = remote.pts_version

    if skirt_version is not None: skirt_versions[host_id] = skirt_version
    if pts_version is not None: pts_versions[host_id] = pts_version

    # Get Qt version


    # Get C++ compiler and MPI compiler version


print("")
print(fmt.bold + fmt.green + "SKIRT" + fmt.reset + ":")
print("")
print(" - local: " + introspection.skirt_version())
for host_id in host_ids:
    if host_id in skirt_versions: print(" - " + host_id + ": " + skirt_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

print(fmt.bold + fmt.green + "PTS" + fmt.reset + ":")
print("")
print(" - local: " + introspection.pts_version())
for host_id in host_ids:
    if host_id in pts_versions: print(" - " + host_id + ": " + pts_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

# -----------------------------------------------------------------
