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
from pts.core.remote.host import find_host_ids, find_hosts
from pts.core.remote.remote import Remote
from pts.core.tools import formatting as fmt
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.load import show_status
from pts.core.basics.log import log, setup_log

# -----------------------------------------------------------------

all_host_ids = find_host_ids()
all_hosts = find_hosts()

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Remote hosts
definition.add_positional_optional("hosts", "host_list", "remote hosts", default=all_hosts, choices=all_host_ids)

# Flags
definition.add_flag("clusters", "show the clusters")

# Create the configuration
config = parse_arguments("remotes", definition, "Check the status of the remotes")

# -----------------------------------------------------------------

# Set log level in a special way
if config.debug: setup_log("DEBUG")
else: setup_log("ERROR")

# -----------------------------------------------------------------

# Loop over the hosts
print("")
for host in config.hosts:

    # Debugging
    log.debug("Connecting to remote host '" + host.id + "' ...")

    # Create remote, try connecting
    remote = Remote()
    connected = remote.setup(host, login_timeout=10) # only try for 10 seconds

    # Connection succeeded?
    if connected: print(fmt.green + host.id + ": up" + fmt.reset)
    else:
        print(fmt.red + host.id + ": down" + fmt.reset)
        continue

    # Show clustername
    if remote.cluster_name is not None: print("Cluster name: " + remote.cluster_name)

    print("")
    print(" - " + fmt.bold + "uses scheduling system: " + fmt.reset + str(remote.scheduler))

    # Show status
    show_status(remote)

    # Disk space

    # print(" - " + fmt.bold + "free memory: " + fmt.reset + str(remote.free_memory))
    #
    # print(" - " + fmt.bold + "total disk space: " + fmt.reset + str(remote.total_space))
    # print(" - " + fmt.bold + "free disk space: " + fmt.reset + str(remote.free_space))
    # print(" - " + fmt.bold + "used disk space: " + fmt.reset + str(remote.used_space))
    # print(" - " + fmt.bold + "free disk space in home directory: " + fmt.reset + str(remote.free_space_home_directory))
    #
    # print(" - " + fmt.bold + "free cores: " + fmt.reset + str(remote.free_cores))
    # print(" - " + fmt.bold + "free sockets: " + fmt.reset + str(remote.free_sockets))
    #
    # print(" - " + fmt.bold + "CPU model: " + fmt.reset + str(remote.cpu_model))
    # print(" - " + fmt.bold + "virtual memory per node: " + fmt.reset + str(remote.virtual_memory_per_node))
    # print(" - " + fmt.bold + "nodes: " + fmt.reset + str(remote.nodes))
    #
    # print(" - " + fmt.bold + "multinode: " + fmt.reset + str(remote.is_multinode))
    # print(" - " + fmt.bold + "cores per node: " + fmt.reset + str(remote.cores_per_node))
    # print(" - " + fmt.bold + "threads per core: " + fmt.reset + str(remote.threads_per_core))
    # print(" - " + fmt.bold + "sockets per node: " + fmt.reset + str(remote.sockets_per_node))
    # print(" - " + fmt.bold + "cores per socket: " + fmt.reset + str(remote.cores_per_socket))
    #
    # print(" - " + fmt.bold + "numa domains: " + fmt.reset + str(remote.numa_domains))
    # print(" - " + fmt.bold + "numa cpus: " + fmt.reset + str(remote.numa_cpus))
    #
    # print(" - " + fmt.bold + "cpu load: " + fmt.reset + str(remote.cpu_load))
    # print(" - " + fmt.bold + "memory load: " + fmt.reset + str(remote.memory_load))
    # print(" - " + fmt.bold + "platform: " + fmt.reset + str(remote.platform))
    # print(" - " + fmt.bold + "architecture: " + fmt.reset + str(remote.architecture))
    # print("")

    # Show clusters
    if config.clusters and host.has_clusters:

        print(" - " + fmt.bold + "default cluster: " + fmt.reset + host.clusters.default)
        print(" - " + fmt.bold + "Other clusters:" + fmt.reset)
        print("")

        for cluster_name in host.clusters:

            if cluster_name == "default": continue
            if cluster_name == host.cluster_name: continue

            print(fmt.yellow + "    " + cluster_name + fmt.reset)
            print("")

            info = host.clusters[cluster_name]
            #print(info)

            numa_domains_per_node = info.numa_domains_per_node
            cores_per_socket = info.cores_per_socket
            sockets_per_node = info.sockets_per_node
            multi_node_communication = info.multi_node_communication
            memory = info.memory
            nodes = info.nodes
            threads_per_core = info.threads_per_core

            # Show
            print("    - " + fmt.bold + "number of NUMA domains per node: " + fmt.reset + str(numa_domains_per_node))
            print("    - " + fmt.bold + "number of nodes: " + fmt.reset + str(nodes))
            print("    - " + fmt.bold + "number of sockets per node: " + fmt.reset + str(sockets_per_node))
            print("    - " + fmt.bold + "number of cores per socket: " + fmt.reset + str(cores_per_socket))
            print("    - " + fmt.bold + "number of threads per core: " + fmt.reset + str(threads_per_core))
            print("    - " + fmt.bold + "suitable for multinode communication: " + fmt.reset + str(multi_node_communication))
            print("    - " + fmt.bold + "virtual memory per node: " + fmt.reset + str(memory))
            print("")

# -----------------------------------------------------------------

# #             # Report
#             print("")
#             print(fmt.bold + fmt.green + host_id + fmt.reset + ":")
#             print("")
#             print(" - Architecture:", self.architecture_dict[host_id])
#             print(" - Operating system:", self.os_dict[host_id])
#             print(" - CPU model:", self.models_dict[host_id])
#             print(" - Number of nodes:", self.nnodes_dict[host_id])
#             print(" - Number of sockets per node:", self.nsockets_dict[host_id])
#             print(" - Number of NUMA domains per node:", self.ndomains_dict[host_id])
#             print(" - Number of cores per socket:", self.ncores_dict[host_id])
#             print(" - Number of threads per core:", self.nthreads_dict[host_id])
#             print(" - Virtual memory per node:", self.memory_dict[host_id])
#             print(" - CPU load:", self.cpu_loads[host_id])
#             print(" - Memory load:", self.memory_loads[host_id])
#             print(" - Number of free nodes:", self.nfreenodes_dict[host_id])
#             print("")

# -----------------------------------------------------------------
