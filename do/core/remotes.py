#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.remotes Check the status of the remotes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.tools import formatting as fmt
from pts.core.basics.log import setup_log
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
host_ids = find_host_ids()
definition.add_positional_optional("hosts", "string_list", "host IDs", default=host_ids, choices=host_ids)
config = parse_arguments("remotes", definition, "Check the status of the remotes")

# -----------------------------------------------------------------

# Set log level in a special way
if config.debug: setup_log("DEBUG")
else: setup_log("ERROR")

# -----------------------------------------------------------------

# Loop over the hosts
print("")
for host_id in config.hosts:

    # Create remote, try connecting
    remote = Remote()
    connected = remote.setup(host_id, login_timeout=15) # only try for 15 seconds

    # Connection succeeded?
    if connected: print(fmt.green + host_id + ": up" + fmt.reset)
    else:
        print(fmt.red + host_id + ": down" + fmt.reset)
        continue

    print("")
    print(" - " + fmt.bold + "uses scheduling system: " + fmt.reset + str(remote.scheduler))
    print(" - " + fmt.bold + "free memory: " + fmt.reset + str(remote.free_memory))
    print(" - " + fmt.bold + "total disk space: " + fmt.reset + str(remote.total_space))
    print(" - " + fmt.bold + "free disk space: " + fmt.reset + str(remote.free_space))
    print(" - " + fmt.bold + "used disk space: " + fmt.reset + str(remote.used_space))
    print(" - " + fmt.bold + "free disk space in home directory: " + fmt.reset + str(remote.free_space_home_directory))
    print(" - " + fmt.bold + "free cores: " + fmt.reset + str(remote.free_cores))
    print(" - " + fmt.bold + "free sockets: " + fmt.reset + str(remote.free_sockets))
    print(" - " + fmt.bold + "cpu_model: " + fmt.reset + str(remote.cpu_model))
    print(" - " + fmt.bold + "virtual_memory_per_node: " + fmt.reset + str(remote.virtual_memory_per_node))
    print(" - " + fmt.bold + "nodes: " + fmt.reset + str(remote.nodes))
    print(" - " + fmt.bold + "is_multinode: " + fmt.reset + str(remote.is_multinode))
    print(" - " + fmt.bold + "cores_per_node: " + fmt.reset + str(remote.cores_per_node))
    print(" - " + fmt.bold + "threads_per_core: " + fmt.reset + str(remote.threads_per_core))
    print(" - " + fmt.bold + "sockets_per_node: " + fmt.reset + str(remote.sockets_per_node))
    print(" - " + fmt.bold + "cores_per_socket: " + fmt.reset + str(remote.cores_per_socket))
    print(" - " + fmt.bold + "numa_domains: " + fmt.reset + str(remote.numa_domains))
    print(" - " + fmt.bold + "numa_cpus: " + fmt.reset + str(remote.numa_cpus))
    print(" - " + fmt.bold + "cpu_load: " + fmt.reset + str(remote.cpu_load))
    print(" - " + fmt.bold + "memory_load: " + fmt.reset + str(remote.memory_load))
    print(" - " + fmt.bold + "platform: " + fmt.reset + str(remote.platform))
    print(" - " + fmt.bold + "architecture: " + fmt.reset + str(remote.architecture))
    print("")

# -----------------------------------------------------------------
