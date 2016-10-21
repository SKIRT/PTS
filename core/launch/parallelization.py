#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.parallelization Contains the ParallelizationTool class and the Parallelization class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import random

# Import the relevant PTS classes and modules
from ..tools.logging import log

# -----------------------------------------------------------------

class Parallelization(object):

    """
    This class ...
    """

    def __init__(self, cores, threads_per_core, processes, data_parallel=False):

        """
        This function ...
        :param cores:
        :param threads_per_core:
        :param processes:
        :param data_parallel:
        """

        # The number of cores
        self.cores = cores

        # The number of threads used per core
        self.threads_per_core = threads_per_core

        # The number of processes
        self.processes = processes

        # Data parallelization mode
        self.data_parallel = data_parallel

    # -----------------------------------------------------------------

    @property
    def threads(self):

        """
        This function ...
        :return:
        """

        threads = self.cores_per_process * self.threads_per_core
        assert int(threads) == threads
        return int(threads)

    # -----------------------------------------------------------------

    @property
    def cores_per_process(self):

        """
        This function ...
        :return:
        """

        corespp = self.cores / self.processes
        assert int(corespp) == corespp
        return int(corespp)

    # -----------------------------------------------------------------

    @classmethod
    def for_host(cls, host, nnodes, processes_per_socket=1, data_parallel=False):

        """
        This function ...
        :param host:
        :param nnodes:
        :param processes_per_socket:
        :param data_parallel:
        :return:
        """

        # Debugging
        log.debug("Determining the parallelization scheme for host " + host.id + " ...")

        # Get the number of cores per node for this host
        cores_per_node = host.cluster.cores_per_socket * host.cluster.sockets_per_node

        # Determine the number of cores corresponding to the number of requested cores
        cores = cores_per_node * nnodes

        # Determine the total number of processes based on the number of process per socket (assume there is enough memory)
        processes_per_node = processes_per_socket * host.cluster.sockets_per_node
        processes = processes_per_node * nnodes

        # Determine the number of threads per core
        if host.use_hyperthreading: threads_per_core = host.cluster.threads_per_core
        else: threads_per_core = 1

        # Create a Parallelization instance
        parallelization = cls(cores, threads_per_core, processes, data_parallel)

        # Return the parallelization scheme
        return parallelization

    # -----------------------------------------------------------------

    @classmethod
    def from_processes_and_threads(cls, processes, threads, threads_per_core=1, data_parallel=False):

        """
        This function ...
        :param processes:
        :param threads:
        :param threads_per_core:
        :param data_parallel:
        :return:
        """

        # Determine the number of required cores
        cores_per_process = threads / threads_per_core
        cores = cores_per_process * processes

        # Create a new class instance and return it
        return cls(cores, threads_per_core, processes, data_parallel)

    # -----------------------------------------------------------------

    @classmethod
    def from_mode(cls, mode, cores, threads_per_core, threads_per_process=None, data_parallel=False):

        """
        This function calculates the number of processes and the number of threads (per process) for
        a certain number of cores, depending on the mode of parallelization (pure 'mpi', pure 'threads' or 'hybrid').
        In other words, this function determines the 'mapping' from a number of cores to an appropriate
        set of threads and processes.
        :param mode:
        :param cores:
        :param threads_per_core:
        :param threads_per_process:
        :param data_parallel:
        :return:
        """

        # Set default values for the number of threads and processes
        processes = 1
        used_threads_per_core = 1

        # In mpi mode, each processor runs a different process
        if mode == "mpi": processes = cores

        # In threads mode, each processor runs a seperate thread within the same process
        if mode == "threads": used_threads_per_core = threads_per_core

        # In hybrid mode, the number of processes depends on how many threads are requested per process
        # and the current number of processors
        if mode == "hybrid":

            # Determine the number of processes
            cores_per_process = threads_per_process / threads_per_core
            processes = cores // cores_per_process

            # Hyperthreading
            used_threads_per_core = threads_per_core

        # Create a new class instance and return it
        return cls(cores, used_threads_per_core, processes, data_parallel)

    # -----------------------------------------------------------------

    @classmethod
    def from_free_cores(cls, free_cores, cores_per_process, threads_per_core, data_parallel=False):

        """
        This function ...
        :param free_cores:
        :param cores_per_process:
        :param threads_per_core:
        :param data_parallel:
        :return:
        """

        # Using the number of cores per process, determine the number of processes
        processes = int(free_cores / cores_per_process)

        # The actual number of cores
        cores = processes * cores_per_process

        # Create a new class instance and return it
        return cls(cores, threads_per_core, processes, data_parallel)

    # -----------------------------------------------------------------

    def get_requirements(self, cores_per_node):

        """
        This function ...
        :param cores_per_node:
        :return:
        """

        # Calculate the necessary amount of nodes
        nodes = self.cores // cores_per_node + (self.cores % cores_per_node > 0)

        # Determine the number of processors per node
        ppn = self.cores if nodes == 1 else cores_per_node

        # Return the required number of nodes and cores per node
        return nodes, ppn

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self.cores == other.cores and self.threads_per_core == other.threads_per_core and self.processes == other.processes

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        """

        return self.__class__.__name__ + " scheme with " + str(self.processes) + " processes and " \
               + str(self.threads_per_core) + " threads per core on a total of " + str(self.cores) + " cores " \
               + "(" + str(self.threads) + " threads per process)"

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        """

        return '<' + self.__class__.__name__ + " cores: " + str(self.cores) + ", threads per core: " \
               + str(self.threads_per_core) + ", processes: " + str(self.processes) + ">"

# -----------------------------------------------------------------
