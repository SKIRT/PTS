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
import random

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from .memoryestimator import MemoryEstimator, estimate_memory
from ..simulation.skifile import SkiFile
from ..simulation.parallelization import Parallelization
from ..tools.utils import lazyproperty
from ..tools import filesystem as fs
from ..tools import sequences

# -----------------------------------------------------------------

def factors(n):

    """
    This function ...
    :param n:
    :return:
    """

    return list(sorted(set(reduce(list.__add__, ([i, n // i] for i in range(1, int(n ** 0.5) + 1) if n % i == 0)))))

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

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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

    def set_parallelization_multiprocessing(self):

        """
        This function ...
        :return:
        """

        # If memory as not been set, estimate it
        if self.memory is None: self.memory = estimate_memory(self.ski, input_path=self.config.input, ncells=self.config.ncells)

        # Get serial and parallel parts of the memory requirement
        serial_memory = self.memory.serial
        parallel_memory = self.memory.parallel

        # Calculate the total memory of one process without data parallelization
        total_memory = serial_memory + parallel_memory

        # The total memory consumption is larger than the node memory
        if total_memory > self.config.memory:

            # Check data parallel flag
            if self.config.data_parallel is False: raise ValueError("Data-parallelization is disabled but memory usage is larger than node memory")

            # Check number of nodes
            if self.config.nnodes == 1: raise ValueError("Simulation cannot be run: decrease resolution, use system with more memory per node, or use more nodes")

            # Set parallelization
            self.set_parallelization_high_memory(serial_memory, parallel_memory) # Ms > Mn

        # If more than one simulation memory fits on a node
        else: self.set_parallelization_low_memory(total_memory)

    # -----------------------------------------------------------------

    def set_parallelization_high_memory(self, serial_memory, parallel_memory):

        """
        This function ...
        :param serial_memory:
        :param parallel_memory:
        :return:
        """

        # TODO: use max_nthreads
        # TODO: don't use hyperthreading for SKIRT by default
        # TODO: use newer cpus-per-proc syntax

        # Debugging
        log.debug("Setting parallelization for simulation with high memory requirement ...")

        # Determine the amount of required memory per process
        memory_per_process = serial_memory + parallel_memory / self.config.nnodes

        # Mn x Nn < Ms?
        # if self.config.memory * self.config.nnodes < memory_simulation_nnode_processes:
        if memory_per_process > self.config.memory:  # is the memory used per process larger than the memory per node?
            raise ValueError("Simulation cannot be run: decrease resolution, use system with more memory per node, or use more nodes")

        # Arbitrarily pick a divisor (DIV) of Npps
        divisors = factors(self.config.ncores)

        # Loop over divisors, try to set parallelization
        for divisor in sequences.iterate_from_middle(divisors):

            # Try to set parallelization
            if self.try_parallelization_divisor(divisor): break

    # -----------------------------------------------------------------

    def try_parallelization_divisor(self, divisor):

        """
        This function ...
        :param divisor:
        :return:
        """

        # Debuggging
        log.debug("Trying to set paralellelization with divisor " + str(divisor) + " ...")

        # First implementation was random
        #divisor = random.choice(divisors)

        # Nt = min(Npps, DIV)
        nthreads = min(self.config.ncores, divisor)

        # Np = Nn x Nppn / Nt
        ppn = self.config.nsockets * self.config.ncores
        nprocesses = self.config.nnodes * ppn / nthreads

        total_ncores = self.config.nnodes * self.config.nsockets * self.config.ncores

        # Nlambda >= 10 x Np
        if self.nwavelengths >= 10 * nprocesses:

            # Determine number of threads per core
            threads_per_core = self.config.threads_per_core if self.config.hyperthreading else 1

            # data parallelization
            # Create the parallelization object
            self.parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core,
                                                             threads_per_process=nthreads,
                                                             data_parallel=True)

            # Success
            return True

        # Try again from picking divisor, but now a larger one (less processes)
        else: return False

    # -----------------------------------------------------------------

    def set_parallelization_low_memory(self, total_memory):

        """
        This function ...
        :param total_memory:
        :return:
        """

        # Debugging
        log.debug("Setting parallelization for simulation with low memory requirement ...")

        # Np = min(Mn / Ms, Nppn)
        ppn = self.config.nsockets * self.config.ncores
        nprocesses_per_node = int(math.floor(min(self.config.memory / total_memory, ppn)))

        nprocesses = nprocesses_per_node * self.config.nnodes

        # nthreads = ppn / nprocesses_per_node
        ncores_per_process = ppn / nprocesses_per_node
        # ncores_per_process = int(math.ceil(ncores_per_process)) -> this leads to total_ncores that is larger than self.config.nnodes * self.config.nsockets * self.config.ncores
        ncores_per_process = int(ncores_per_process)

        # Determine number of threads per core
        threads_per_core = self.config.threads_per_core if self.config.hyperthreading else 1

        # Threads per process
        threads_per_process = threads_per_core * ncores_per_process

        # Total number of cores
        total_ncores = nprocesses * ncores_per_process

        # Data-parallelization flag is set
        if self.config.data_parallel is not None: data_parallel = self.config.data_parallel

        # Nlambda >= 10 * Np?
        # nwavelengths = self.ski.nwavelengthsfile(self.config.input) if self.ski.wavelengthsfile() else self.ski.nwavelengths()
        elif self.nwavelengths >= 10 * nprocesses and self.dustlib_dimension == 3: data_parallel = True

        # No data-parallelization should be used
        else: data_parallel = False

        # data parallelization
        # Create the parallelization object
        self.parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core,
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
