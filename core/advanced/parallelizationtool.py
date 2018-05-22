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
import math

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from .memoryestimator import estimate_memory
from ..simulation.skifile import SkiFile
from ..simulation.parallelization import Parallelization
from ..tools.utils import lazyproperty
from ..tools import filesystem as fs
from ..tools import numbers

# -----------------------------------------------------------------

def factors(n, reversed=False):

    """
    This function ...
    :param n:
    :param reversed:
    :return:
    """

    return list(sorted(set(reduce(list.__add__, ([i, n // i] for i in range(1, int(n ** 0.5) + 1) if n % i == 0))), reverse=reversed))

# -----------------------------------------------------------------

# Parameters important for determining parallelization scheme:

# input", "directory_path", "path to the input directory")

# nnodes", "integer", "number of nodes")
# nsockets", "integer", "number of sockets per node")
# ncores", "integer", "number of cores per socket")
# memory", "real", "available virtual memory per node")

# mpi", "mpi available", True)
# hyperthreading", "use hyperthreading", False)
# threads_per_core", "integer", "number of hyperthreads per core")

# Add optional
# ncells", "integer", "number of dust cells (relevant if ski file uses a tree dust grid)")

## MEMORYTOOL:

#

# -----------------------------------------------------------------

def determine_parallelization(ski_path, input_path, memory, nnodes, nsockets, ncores, host_memory, mpi, hyperthreading,
                              threads_per_core, ncells=None, nwavelengths=None, data_parallel=None):

    """
    This function ...
    :param ski_path:
    :param input_path:
    :param memory:
    :param nnodes:
    :param nsockets:
    :param ncores:
    :param host_memory:
    :param mpi:
    :param hyperthreading:
    :param threads_per_core:
    :param ncells:
    :param nwavelengths:
    :param data_parallel:
    :return:
    """

    # Debug
    prefix = fs.strip_extension(fs.name(ski_path))
    log.debug("Determining parallelization scheme for the '" + prefix + "' simulation with the following parameters:")
    log.debug("")
    if memory is not None: log.debug(" - memory requirement: " + str(memory))
    log.debug(" - nnodes: " + str(nnodes))
    log.debug(" - nsockets: " + str(nsockets))
    log.debug(" - ncores: " + str(ncores))
    log.debug(" - host memory: " + str(host_memory))
    log.debug(" - mpi: " + str(mpi))
    log.debug(" - hyperthreading: " + str(hyperthreading))
    log.debug(" - threads per core: " + str(threads_per_core))
    if ncells is not None: log.debug(" - ncells: " + str(ncells))
    if nwavelengths is not None: log.debug(" - nwavelengths: " + str(nwavelengths))
    log.debug("")

    # Create the parallelization tool
    tool = ParallelizationTool()

    # Set configuration options
    tool.config.ski = ski_path
    tool.config.input = input_path

    # Set host properties
    tool.config.nnodes = nnodes
    tool.config.nsockets = nsockets
    tool.config.ncores = ncores
    tool.config.memory = host_memory

    # MPI available and used
    tool.config.mpi = mpi
    tool.config.hyperthreading = hyperthreading
    tool.config.threads_per_core = threads_per_core

    # Number of dust cells
    tool.config.ncells = ncells  # number of dust cells (relevant if ski file uses a tree dust grid)

    # Number of wavelengths
    tool.config.nwavelengths = nwavelengths

    # Data-parallelization
    tool.config.data_parallel = data_parallel

    # Don't show the parallelization
    tool.config.show = False

    # Run the parallelization tool (passing the memory requirement of the simulation as an argument)
    tool.run(memory=memory)

    # Get the parallelization scheme
    return tool.parallelization

# -----------------------------------------------------------------

class ParallelizationTool(Configurable):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ParallelizationTool, self).__init__(*args, **kwargs)

        # The ski file
        self.ski = None

        # The dimension of the dust lib
        self.dustlib_dimension = None

        # The parallelization object
        self.parallelization = None

        # The memory requirement for the simulation
        self.memory = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 2. Get additional properties from the ski file
        self.get_properties()

        # 3. Set the parallelization scheme
        self.set_parallelization()

        # 4. Show the parallelization scheme
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ParallelizationTool, self).setup(**kwargs)

        # Open the ski file
        self.ski = self.config.ski if isinstance(self.config.ski, SkiFile) else SkiFile(self.config.ski)

        # Set the memory requirement
        self.memory = kwargs.pop("memory", None)

    # -----------------------------------------------------------------

    def get_properties(self):

        """
        This function ...
        :return:
        """

        # Get the dustlib dimension
        self.dustlib_dimension = self.ski.dustlib_dimension()

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        if self.config.nwavelengths is not None: return self.config.nwavelengths
        else: return self.ski.nwavelengthsfile(self.config.input) if self.ski.wavelengthsfile() else self.ski.nwavelengths()

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # MPI cannot be used
        if not self.config.mpi: self.set_parallelization_singleprocessing()

        # MPI can be used
        else: self.set_parallelization_multiprocessing()

    # -----------------------------------------------------------------

    def set_parallelization_singleprocessing(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Debugging
        log.debug("Setting parallelization scheme for singleprocessing ...")

        # Determine the number of cores per node
        cores_per_node = self.config.nsockets * self.config.ncores

        # Determine optimal number of cores used for threading (not more than 12)
        cores = min(cores_per_node, 12)

        # if self.config.hyperthreading: threads = cores * self.config.threads_per_core
        # else: threads = cores

        # Determine number of threads per core
        threads_per_core = self.config.threads_per_core if self.config.hyperthreading else 1

        # Create the parallelization object
        self.parallelization = Parallelization.from_mode("threads", cores, threads_per_core)

    # -----------------------------------------------------------------

    @property
    def serial_memory(self):

        """
        This function ...
        :return:
        """

        return self.memory.serial

    # -----------------------------------------------------------------

    @property
    def parallel_memory(self):

        """
        This function ...
        :return:
        """

        return self.memory.parallel

    # -----------------------------------------------------------------

    @property
    def total_memory(self):

        """
        This function ...
        :return:
        """

        # the total memory of one process without data parallelization
        return self.serial_memory + self.parallel_memory

    # -----------------------------------------------------------------

    @lazyproperty
    def nthreads_per_core(self):

        """
        This function ...
        :return:
        """

        # Determine number of threads per core
        return self.config.threads_per_core if self.config.hyperthreading else 1

    # -----------------------------------------------------------------

    @lazyproperty
    def max_ncores_per_process(self):

        """
        This function ...
        :return:
        """

        return int(math.ceil(self.config.max_nthreads / self.nthreads_per_core)) # round up because when hyperthreading is envolved, we can as well allow a slightly higher effective nthreads per process

    # -----------------------------------------------------------------

    @lazyproperty
    def min_nprocesses_per_socket(self):

        """
        This function ...
        :return:
        """

        return numbers.as_integer_check(self.config.ncores / self.max_ncores_per_process)

    # -----------------------------------------------------------------

    @lazyproperty
    def ideal_ncores_per_process(self):

        """
        This function ...
        :return:
        """

        # Lesser cores per socket than the maximum number of cores per process
        if self.config.ncores <= self.max_ncores_per_process: return self.config.ncores
        else:

            # The number of cores per process must be a divisor of the number of cores per socket (so that the number of processes is integer)
            divisors = factors(self.config.ncores)
            allowed_divisors = [divisor for divisor in divisors if divisor <= self.max_ncores_per_process]

            # Pick the largest divisor that is still smaller than the maximum number of cores per process
            return max(allowed_divisors)

    # -----------------------------------------------------------------

    @lazyproperty
    def ideal_nprocesses_per_socket(self):

        """
        This function ...
        :return:
        """

        return numbers.as_integer_check(self.config.ncores / self.ideal_ncores_per_process)
        #if not numbers.is_integer(result): raise RuntimeError("Something went wrong")
        #return int(result)

    # -----------------------------------------------------------------

    @lazyproperty
    def ideal_nprocesses_per_node(self):

        """
        This function ...
        :return:
        """

        return self.ideal_nprocesses_per_socket * self.config.nsockets

    # -----------------------------------------------------------------

    def set_parallelization_multiprocessing(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting parallelization scheme for multiprocessing ...")

        # If memory as not been set, estimate it
        if self.memory is None: self.memory = estimate_memory(self.ski, input_path=self.config.input, ncells=self.config.ncells)

        # The total memory consumption is larger than the node memory
        if self.total_memory > self.config.memory:

            # Check data parallel flag
            if self.config.data_parallel is False: raise ValueError("Data-parallelization is disabled but memory usage is larger than node memory")

            # Check number of nodes
            if self.config.nnodes == 1: raise ValueError("Simulation cannot be run: decrease resolution, use system with more memory per node, or use more nodes")

            # Set parallelization
            self.set_parallelization_high_memory() # Ms > Mn

        # If more than one simulation memory fits on a node
        else: self.set_parallelization_low_memory()

    # -----------------------------------------------------------------

    def set_parallelization_high_memory(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting parallelization for simulation with high memory requirement ...")

        # Determine the amount of required memory per process
        memory_per_process = self.serial_memory + self.parallel_memory / self.config.nnodes

        # Mn x Nn < Ms?
        # if self.config.memory * self.config.nnodes < memory_simulation_nnode_processes:
        if memory_per_process > self.config.memory:  # is the memory used per process larger than the memory per node?
            raise ValueError("Simulation cannot be run: decrease resolution, use system with more memory per node, or use more nodes")

        # Try setting a parallelization with the ideal number of processes per socket
        success = self.try_parallelization_pps(self.ideal_nprocesses_per_socket)
        if success: return

        # Try other nprocesses per socket
        max_nprocesses_per_socket = self.ideal_nprocesses_per_socket
        min_nprocesses_per_socket = self.min_nprocesses_per_socket

        # Factors of the number of cores per socket, from highest to lowest
        divisors = factors(self.config.ncores, reversed=True)

        # Loop over the divisors
        for pps in divisors:

            # Check
            if pps >= max_nprocesses_per_socket: continue
            if pps < min_nprocesses_per_socket:
                log.warning("Could not find an efficient parallelization scheme for data parallelization")
                self.set_parallelization_pps(min_nprocesses_per_socket)

            # Try
            success = self.try_parallelization_pps(pps)
            if success: break

    # -----------------------------------------------------------------

    @property
    def total_ncores(self):

        """
        This function ...
        :return:
        """

        return self.config.nnodes * self.config.nsockets * self.config.ncores

    # -----------------------------------------------------------------

    def try_parallelization_pps(self, processes_per_socket):

        """
        This function ...
        :param processes_per_socket:
        :return:
        """

        # Debuggging
        log.debug("Trying to set paralellelization with number of processes per socket of " + str(processes_per_socket) + " ...")

        # Determine number of cores per process
        cores_per_process = numbers.as_integer_check(self.config.ncores / processes_per_socket)

        # Determine number of threads per process
        threads_per_process = cores_per_process * self.nthreads_per_core

        # Determine total number of processes
        nprocesses = processes_per_socket * self.config.ncores * self.config.nnodes

        # Nlambda >= 10 x Np: can we use data parallelization?
        # -> nprocesses cannot be too high (compared to nwavelengths) for load balancing
        if self.nwavelengths >= self.config.min_nwavelengths_per_process * nprocesses:

             # data parallelization
             # Create the parallelization object
             self.parallelization = Parallelization.from_mode("hybrid", self.total_ncores, self.nthreads_per_core,
                                                              threads_per_process=threads_per_process,
                                                              data_parallel=True)

             # Success
             return True

        # Try again with less processes?
        else: return False

    # -----------------------------------------------------------------

    def set_parallelization_pps(self, processes_per_socket):

        """
        This function ...
        :param processes_per_socket:
        :return:
        """

        # Debugging
        log.debug("Setting parallelization scheme with " + str(processes_per_socket) + " processes per socket ...")

        # Determine number of cores per process
        cores_per_process = numbers.as_integer_check(self.config.ncores / processes_per_socket)

        # Determine number of threads per process
        threads_per_process = cores_per_process * self.nthreads_per_core

        # Determine total number of processes
        nprocesses = processes_per_socket * self.config.ncores * self.config.nnodes

        # Show the number of wavelengths per process
        nwavelengths_per_process = self.nwavelengths / nprocesses
        log.debug("The number of wavelengths per process is " + str(nwavelengths_per_process))

        # data parallelization
        # Create the parallelization object
        self.parallelization = Parallelization.from_mode("hybrid", self.total_ncores, self.nthreads_per_core,
                                                         threads_per_process=threads_per_process,
                                                         data_parallel=True)

    # -----------------------------------------------------------------

    def set_parallelization_low_memory(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting parallelization for simulation with low memory requirement ...")

        # Np = min(Mn / Ms, Nppn)
        ppn = self.config.nsockets * self.config.ncores

        # Determine the optimal number of processes per node, taking into account memory considerations and not having too many processes
        nprocesses_per_node = int(math.floor(min(self.config.memory / self.total_memory, ppn)))
        nprocesses_per_node = min(self.ideal_nprocesses_per_node, nprocesses_per_node)

        # Determine the total number of processes
        nprocesses = nprocesses_per_node * self.config.nnodes

        # nthreads = ppn / nprocesses_per_node
        ncores_per_process = ppn / nprocesses_per_node
        # ncores_per_process = int(math.ceil(ncores_per_process)) -> this leads to total_ncores that is larger than self.config.nnodes * self.config.nsockets * self.config.ncores
        ncores_per_process = int(ncores_per_process)

        # Threads per process
        threads_per_process = self.nthreads_per_core * ncores_per_process

        # Total number of cores
        total_ncores = nprocesses * ncores_per_process

        # Data-parallelization flag is set
        if self.config.data_parallel is not None: data_parallel = self.config.data_parallel

        # Nlambda >= 10 * Np?
        # nwavelengths = self.ski.nwavelengthsfile(self.config.input) if self.ski.wavelengthsfile() else self.ski.nwavelengths()
        elif self.nwavelengths >= self.config.min_nwavelengths_per_process * nprocesses and self.dustlib_dimension == 3: data_parallel = True

        # No data-parallelization should be used
        else: data_parallel = False

        # data parallelization
        # Create the parallelization object
        self.parallelization = Parallelization.from_mode("hybrid", total_ncores, self.nthreads_per_core,
                                                         threads_per_process=threads_per_process,
                                                         data_parallel=data_parallel)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        if self.parallelization is not None:

            log.info("The paralleliation scheme is:")
            print(self.parallelization)

# -----------------------------------------------------------------
