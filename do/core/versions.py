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
from pts.core.remote.remote import Remote, HostDownException
from pts.core.tools.logging import log

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

python_versions = dict()

cpp_versions = dict()
mpi_versions = dict()
qt_versions = dict()

skirt_versions = dict()
pts_versions = dict()

down_hosts = []

# Loop over the different hosts
for host_id in host_ids:

    # Setup the remote (login)
    remote = Remote()
    try: remote.setup(host_id)
    except HostDownException:
        log.warning("Remote host '" + host_id + "' is down: skipping")
        down_hosts.append(host_id)
        continue

    # Loading compilers
    compiler_path = remote.find_and_load_cpp_compiler()
    mpi_compiler_path = remote.find_and_load_mpi_compiler()

    print(compiler_path)
    print(mpi_compiler_path)

    # Get C++ compiler and MPI compiler version
    compiler_version = remote.version_of(compiler_path) if compiler_path is not None else None
    mpi_compiler_version = remote.version_of(mpi_compiler_path) if mpi_compiler_path is not None else None

    # Load Qt module, find the qmake path
    qmake_path = remote.find_and_load_qmake()

    print(qmake_path)

    # Get qmake version
    qmake_version = remote.version_of(qmake_path) if qmake_path is not None else None

    # Get python version
    python_version = remote.version_of("python")

    # Get SKIRT and PTS version
    skirt_version = remote.skirt_version
    pts_version = remote.pts_version

    # Set versions
    if compiler_version is not None: cpp_versions[host_id] = compiler_version
    if mpi_compiler_version is not None: mpi_versions[host_id] = mpi_compiler_version
    if qmake_version is not None: qt_versions[host_id] = qmake_version
    if python_version is not None: python_versions[host_id] = python_version
    if skirt_version is not None: skirt_versions[host_id] = skirt_version
    if pts_version is not None: pts_versions[host_id] = pts_version

# -----------------------------------------------------------------

print("")
print(fmt.bold + fmt.green + "SKIRT" + fmt.reset + ":")
print("")
print(" - local: " + introspection.skirt_version())
for host_id in host_ids:
    if host_id in down_hosts: continue
    if host_id in skirt_versions: print(" - " + host_id + ": " + skirt_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

print(fmt.bold + fmt.green + "PTS" + fmt.reset + ":")
print("")
print(" - local: " + introspection.pts_version())
for host_id in host_ids:
    if host_id in down_hosts: continue
    if host_id in pts_versions: print(" - " + host_id + ": " + pts_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

print(fmt.bold + fmt.green + "C++ compiler" + fmt.reset + ":")
print("")
print(" - local: " + " ...")
for host_id in host_ids:
    if host_id in down_hosts: continue
    if host_id in cpp_versions: print(" - " + host_id + ": " + cpp_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

print(fmt.bold + fmt.green + "MPI compiler" + fmt.reset + ":")
print("")
print(" - local: " + "...")
for host_id in host_ids:
    if host_id in down_hosts: continue
    if host_id in mpi_versions: print(" - " + host_id + ": " + mpi_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

print(fmt.bold + fmt.green + "Qmake" + fmt.reset + ":")
print("")
print(" - local: " + "...")
for host_id in host_ids:
    if host_id in down_hosts: continue
    if host_id in qt_versions: print(" - " + host_id + ": " + qt_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

print(fmt.bold + fmt.green + "Python" + fmt.reset + ":")
print("")
print(" - local: " + introspection.python_version())
for host_id in host_ids:
    if host_id in down_hosts: continue
    if host_id in python_versions: print(" - " + host_id + ": " + python_versions[host_id])
    else: print(" - " + host_id + ": not found")
print("")

# -----------------------------------------------------------------
