#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.load Contains the LoadChecker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import formatting as fmt
from ..tools.logging import log
from .configurable import RemotesConfigurable
from ..basics.unit import parse_unit as u

# -----------------------------------------------------------------

class LoadChecker(RemotesConfigurable):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(LoadChecker, self).__init__(config, interactive)

        # Architecture
        self.architecture_dict = dict()
        self.models_dict = dict()
        self.nnodes_dict = dict()
        self.nfreenodes_dict = dict()
        self.nsockets_dict = dict()
        self.ndomains_dict = dict()
        self.ncores_dict = dict()
        self.nthreads_dict = dict()

        # OS
        self.os_dict = dict()

        # Memory
        self.memory_dict = dict()

        # Load
        self.cpu_loads = dict()
        self.memory_loads = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Check load
        self.check()

        # 3. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(LoadChecker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Checking on each remote host ...")

        # Loop over the remotes
        for remote in self.remotes:

            # Get host ID
            host_id = remote.host_id

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

            else:

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
                nfreenodes = nnodes * (1. - cpu_load/100.)

            # Set properties
            self.architecture_dict[host_id] = architecture
            self.os_dict[host_id] = os
            self.models_dict[host_id] = cpu_model
            self.nnodes_dict[host_id] = nnodes
            self.nfreenodes_dict[host_id] = nfreenodes
            self.nsockets_dict[host_id] = nsockets
            self.ndomains_dict[host_id] = ndomains
            self.ncores_dict[host_id] = ncores
            self.nthreads_dict[host_id] = nthreads
            self.memory_dict[host_id] = memory
            self.cpu_loads[host_id] = cpu_load
            self.memory_loads[host_id] = memory_load

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Loop over hosts, report
        for host_id in self.host_ids:

            # Report
            print("")
            print(fmt.bold + fmt.green + host_id + fmt.reset + ":")
            print("")
            print(" - Architecture:", self.architecture_dict[host_id])
            print(" - Operating system:", self.os_dict[host_id])
            print(" - CPU model:", self.models_dict[host_id])
            print(" - Number of nodes:", self.nnodes_dict[host_id])
            print(" - Number of sockets per node:", self.nsockets_dict[host_id])
            print(" - Number of NUMA domains per node:", self.ndomains_dict[host_id])
            print(" - Number of cores per socket:", self.ncores_dict[host_id])
            print(" - Number of threads per core:", self.nthreads_dict[host_id])
            print(" - Virtual memory per node:", self.memory_dict[host_id])
            print(" - CPU load:", self.cpu_loads[host_id])
            print(" - Memory load:", self.memory_loads[host_id])
            print(" - Number of free nodes:", self.nfreenodes_dict[host_id])
            print("")

# -----------------------------------------------------------------
