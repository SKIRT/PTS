#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.load Check the load of the remote hosts.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.host import find_host_ids
from pts.core.tools import formatting as fmt
from pts.core.remote.remote import Remote, HostDownException
from pts.core.tools.logging import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("host_id", "string", "remote host ID", choices=find_host_ids())
definition.add_optional("cluster_name", "string", "cluster name")

# Read the command line arguments
setter = ArgumentConfigurationSetter("load", "Check the load of the remote hosts")
config = setter.run(definition)

# -----------------------------------------------------------------

if config.remote is not None: host_ids = [config.host_id]
else: host_ids = find_host_ids()

# -----------------------------------------------------------------

down_hosts = []

# Architecture
architecture_dict = dict()
models_dict = dict()
nnodes_dict = dict()
nfreenodes_dict = dict()
nsockets_dict = dict()
ndomains_dict = dict()
ncores_dict = dict()
nthreads_dict = dict()

os_dict = dict()

# Memory
memory_dict = dict()

# Load
cpu_loads = dict()
memory_loads = dict()

# Loop over the remotes
for host_id in host_ids:

    # Setup the remote (login)
    remote = Remote()
    try:
        remote.setup(host_id)
    except HostDownException:
        log.warning("Remote host '" + host_id + "' is down: skipping")
        down_hosts.append(host_id)
        continue

    # Multinode
    if remote.is_multinode:

        # Get node status
        status = remote.get_node_status()

        #print(status)

        from lxml import etree
        result = etree.tostring(status, pretty_print=True, method="html")

        # CPU model
        cpu_model = None

        nsockets_list = []
        ndomains_list = []
        ncores_list = []
        nthreads_list = []
        free_list = []
        memory_list = []
        architecture_list = []
        os_list = []
        free_list = []

        for element in status.getchildren():

            #print(element.tag)

            if element.tag != "Node": continue # doesn't happen

            nnodes += 1

            # Get child element of the given element
            name = element.xpath("name")[0].text

            # Get ...
            state = element.xpath("state")[0].text
            power_state = element.xpath("power_state")[0].text
            np = element.xpath("np")[0].text
            ntype = element.xpath("ntype")[0].text
            status = element.xpath("status")[0].text
            total_sockets = int(element.xpath("total_sockets")[0].text)
            total_numa_nodes = int(element.xpath("total_numa_nodes")[0].text)
            total_cores = int(element.xpath("total_cores")[0].text)
            total_threads = int(element.xpath("total_threads")[0].text)

            nsockets_list.append(total_sockets)
            ndomains_list.append(total_numa_nodes)
            ncores_list.append(total_cores)
            nthreads_list.append(total_threads)

            # Get ...
            status_state = status.split("state=")[1].split(",")[0]

            size = status.split("size=")[1].split(",")[0]
            netload = status.split("netload=")[1].split(",")[0]
            physmem = float(status.split("physmem=")[1].split("kb")[0]) * 1e-6 * Unit("GB")
            availmem = float(status.split("availmem=")[1].split("kb")[0]) * 1e-6 * Unit("GB")
            totmem = float(status.split("totmem=")[1].split("kb")[0]) * 1e-6 * Unit("GB")
            idletime = int(status.split("idletime=")[1].split(",")[0]) * Unit("GB")
            opsys = status.split("opsys=")[1]

            architecture = status.split("uname=")[1].split(" #")[0].split(" ")[2]

            architecture_list.append(architecture)

            if state == "free": free_list.append(True)
            else: free_list.append(False)

            memory_list.append(totmem)

            os_list.append(opsys)

        architecture = architecture_list[0]
        os = os_list[0]
        nnodes = len(architecture_list)
        nsockets = nsockets_list[0]
        ndomains = ndomains_list[0]
        ncores = ncores_list[0]
        nthreads = nthreads_list[0]
        memory = memory_list[0]

        nfreenodes = sum(free_list)

        # CPU load and memory load
        cpu_load = None
        memory_load = None

    else:

        # Architecture
        architecture = remote.architecture

        # OS
        os = remote.operating_system

        # CPU model
        cpu_model = remote.cpu_model

        # Nnodes = 1
        nnodes = 1
        nsockets = remote.sockets_per_node
        ndomains = remote.numa_domains
        ncores = remote.cores_per_socket
        nthreads = remote.threads_per_core

        # Virtual memory per node
        memory = remote.virtual_memory_per_node

        # Get CPU load
        cpu_load = remote.cpu_load * 100.

        # Memory load
        memory_load = remote.memory_load * 100.

        # Free nodes
        nfreenodes = nnodes * (1. - cpu_load/100.)

    # Set properties
    architecture_dict[host_id] = architecture
    os_dict[host_id] = os
    models_dict[host_id] = cpu_model
    nnodes_dict[host_id] = nnodes
    nfreenodes_dict[host_id] = nfreenodes
    nsockets_dict[host_id] = nsockets
    ndomains_dict[host_id] = ndomains
    ncores_dict[host_id] = ncores
    nthreads_dict[host_id] = nthreads
    memory_dict[host_id] = memory
    cpu_loads[host_id] = cpu_load
    memory_loads[host_id] = memory_load

# Loop over hosts, report
for host_id in host_ids:

    # Skip down hosts
    if host_id in down_hosts: continue

    # Report
    print("")
    print(fmt.bold + fmt.green + host_id + fmt.reset + ":")
    print("")
    print(" - Architecture:", architecture_dict[host_id])
    print(" - Operating system:", os_dict[host_id])
    print(" - CPU model:", models_dict[host_id])
    print(" - Number of nodes:", nnodes_dict[host_id])
    print(" - Number of sockets per node:", nsockets_dict[host_id])
    print(" - Number of NUMA domains per node:", ndomains_dict[host_id])
    print(" - Number of cores per socket:", ncores_dict[host_id])
    print(" - Number of threads per core:", nthreads_dict[host_id])
    print(" - Virtual memory per node:", memory_dict[host_id])
    print(" - CPU load:", cpu_loads[host_id])
    print(" - Memory load:", memory_loads[host_id])
    print(" - Number of free nodes:", nfreenodes_dict[host_id])
    print("")

# -----------------------------------------------------------------
