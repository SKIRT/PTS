#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.parallelization Contains the Parallelization class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..tools import monitoring

# -----------------------------------------------------------------

def represent_parallelization(parallelization):

    """
    This function ...
    :param parallelization:
    :return:
    """

    return str(parallelization.cores) + ":" + str(parallelization.processes) + ":" + str(parallelization.threads_per_core)

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
    def ncores(self):

        """
        This function ...
        :return:
        """

        return self.cores

    # -----------------------------------------------------------------

    @ncores.setter
    def ncores(self, value):

        """
        This fucntion ...
        :return:
        """

        self.cores = value

    # -----------------------------------------------------------------

    @property
    def nprocesses(self):

        """
        This function ...
        :return:
        """

        return self.processes

    # -----------------------------------------------------------------

    @nprocesses.setter
    def nprocesses(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.processes = value

    # -----------------------------------------------------------------

    @property
    def threads(self):

        """
        This function ...
        :return:
        """

        #print(self.cores_per_process, self.threads_per_core)
        threads = self.cores_per_process * self.threads_per_core
        assert int(threads) == threads
        return int(threads)

    # -----------------------------------------------------------------

    @property
    def nthreads(self):

        """
        This function ...
        :return:
        """

        return self.threads

    # -----------------------------------------------------------------

    @nthreads.setter
    def nthreads(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Scale the number of cores per process so to keep threads_per_core constant
        factor = value / self.nthreads
        new_ncores = self.ncores * factor
        self.ncores = new_ncores

    # -----------------------------------------------------------------

    @property
    def nthreads_per_core(self):

        """
        This function ...
        :return:
        """

        return self.threads_per_core

    # -----------------------------------------------------------------

    @nthreads_per_core.setter
    def nthreads_per_core(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.threads_per_core = value

    # -----------------------------------------------------------------

    @property
    def cores_per_process(self):

        """
        This function ...
        :return:
        """

        from ..tools import numbers
        #print(self.cores, self.processes)
        corespp = self.cores / self.processes
        if not numbers.is_integer(corespp): raise RuntimeError("Invalid parallelization")
        return int(corespp)

    # -----------------------------------------------------------------

    @property
    def ncores_per_process(self):

        """
        This function ...
        :return:
        """

        return self.cores_per_process

    # -----------------------------------------------------------------

    @classmethod
    def for_local(cls, nprocesses=1, data_parallel=False):

        """
        This function ...
        :param nprocesses:
        :param data_parallel:
        :return:
        """

        # Determine the number of threads
        free_cpus = monitoring.free_cpus()
        ncores = math.floor(free_cpus)

        # Determine nthreads
        if nprocesses == 1: threads = ncores
        else: threads = 1

        # Create the parallelization
        parallelization = cls(threads, 1, nprocesses, data_parallel=data_parallel)

        # Return the parallelization scheme
        return parallelization

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

        from ..tools import numbers

        # Determine the number of required cores
        cores_per_process = threads / threads_per_core
        if not numbers.is_integer(cores_per_process): raise ValueError("Inconsistent values")
        cores_per_process = int(cores_per_process)
        cores = cores_per_process * processes

        # Create a new class instance and return it
        # cores, threads_per_core, processes, data_parallel
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

            # Check
            if cores_per_process != int(cores_per_process):
                raise ValueError("Cannot set parallelization: number of threads per process (" + str(threads_per_process) + ") must be a multiple of the number of threads per core (" + str(threads_per_core) + ")")

            # Set number of processes
            cores_per_process = int(cores_per_process)
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

    def as_tuple(self):

        """
        This function ...
        :return:
        """

        return (self.cores, self.processes, self.threads_per_core, self.data_parallel)

    # -----------------------------------------------------------------

    def as_string(self):

        """
        This function returns a parsable string
        :return:
        """

        return ":".join(str(prop) for prop in self.as_tuple())

    # -----------------------------------------------------------------

    @classmethod
    def from_string(cls, string, default_nthreads_per_core=1):

        """
        This function ...
        :param string:
        :param default_nthreads_per_core:
        :return:
        """

        from ..tools.parsing import integer, boolean

        # Parse
        splitted = string.split(":")
        if len(splitted) < 2: raise ValueError("Invalid input: must be 'ncores:nprocesses[:nthreads_per_core][:data_parallel?]'")
        if len(splitted) > 4: raise ValueError("Invalid input: must be 'ncores:nprocesses[:nthreads_per_core][:data_parallel?]'")

        ncores = integer(splitted[0])
        nprocesses = integer(splitted[1])

        # Get number of threads per core (hyperthreading)
        if len(splitted) > 2: nthreads_per_core = integer(splitted[2])
        else: nthreads_per_core = default_nthreads_per_core

        # Get data_parallel flag
        if len(splitted) > 3: data_parallel = boolean(splitted[3])
        else: data_parallel = False

        # Create and return the parallelization
        return Parallelization(ncores, nthreads_per_core, nprocesses, data_parallel=data_parallel)

    # -----------------------------------------------------------------

    def __hash__(self):

        """
        This function ...
        :return:
        """

        return hash(self.as_tuple())

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

        if self.processes > 1: mode_string = " in task+data parallel mode" if self.data_parallel else " in task parallel mode"
        else: mode_string = ""

        return "parallelization scheme with " + str(self.processes) + " processes and " \
               + str(self.threads_per_core) + " threads per core on a total of " + str(self.cores) + " cores " \
               + "(" + str(self.threads) + " threads per process)" + mode_string

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        """

        return '<' + self.__class__.__name__ + " cores: " + str(self.cores) + ", threads per core: " \
               + str(self.threads_per_core) + ", processes: " + str(self.processes) + ", data_parallel: " + str(self.data_parallel) + ">"

# -----------------------------------------------------------------

def get_possible_nprocesses_in_memory(free_memory, serial_memory, parallel_memory, data_parallel=False):

    """
    This function ...
    :param free_memory:
    :param serial_memory:
    :param parallel_memory:
    :param data_parallel:
    :return:
    """

    # Calculate the total memory of one process without data parallelization
    total_memory = serial_memory + parallel_memory

    # Calculate the maximum number of processes based on the memory requirements
    processes = int(free_memory / total_memory)

    # If there is too little free memory for the simulation, the number of processes will be smaller than one
    if processes < 1:

        # Exit with an error
        raise RuntimeError("Not enough memory available to run this simulation locally: free memory = " + str(
            free_memory) + ", required memory = " + str(total_memory))

    # Otherwise, return the number of processes
    return processes

# -----------------------------------------------------------------
