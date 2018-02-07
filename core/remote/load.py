#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.load Contains the get_status function.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.map import Map
from ..units.parsing import parse_unit as u

# -----------------------------------------------------------------

def get_status(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    # Multinode
    if remote.is_multinode: return get_status_multinode(remote)

    # Singlenode
    else: return get_status_singlenode(remote)

# -----------------------------------------------------------------

def get_status_singlenode(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    # Architecture
    architecture = remote.architecture

    # OS
    os = remote.operating_system_short

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
    nfreenodes = nnodes * (1. - cpu_load / 100.)

    # Create mapping for the properties
    status = Map()

    # Set properties
    status.architecture = architecture
    status.os = os
    status.cpu_model = cpu_model
    status.nnodes = nnodes
    status.nfreenodes = nfreenodes
    status.nsockets = nsockets
    status.ndomains = ndomains
    status.ncores = ncores
    status.nthreads = nthreads
    status.memory = memory
    status.cpu_load = cpu_load
    status.memory_load = memory_load

    # Return
    return status

# -----------------------------------------------------------------

def get_status_multinode(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    # Get node status
    status = remote.get_node_status()

    # print(status)

    from lxml import etree
    result = etree.tostring(status, pretty_print=True, method="html")

    # CPU model
    cpu_model = None

    nsockets_list = []
    ndomains_list = []
    ncores_list = []
    nthreads_list = []
    memory_list = []
    architecture_list = []
    os_list = []
    free_list = []

    for element in status.getchildren():

        # print(element.tag)

        if element.tag != "Node": continue  # doesn't happen

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
        nthreads_list.append(total_threads / total_cores) # threads per core

        # Get ...
        status_state = status.split("state=")[1].split(",")[0]

        size = status.split("size=")[1].split(",")[0]
        netload = status.split("netload=")[1].split(",")[0]
        physmem = float(status.split("physmem=")[1].split("kb")[0]) * 1e-6 * u("GB")
        availmem = float(status.split("availmem=")[1].split("kb")[0]) * 1e-6 * u("GB")
        totmem = float(status.split("totmem=")[1].split("kb")[0]) * 1e-6 * u("GB")
        idletime = int(status.split("idletime=")[1].split(",")[0]) * u("GB")
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

    # Create mapping for the properties
    status = Map()

    # Set properties
    status.architecture = architecture
    status.os = os
    status.cpu_model = cpu_model
    status.nnodes = nnodes
    status.nfreenodes = nfreenodes
    status.nsockets = nsockets
    status.ndomains = ndomains
    status.ncores = ncores
    status.nthreads = nthreads
    status.memory = memory
    status.cpu_load = cpu_load
    status.memory_load = memory_load

    # Return
    return status

# -----------------------------------------------------------------
