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
from ..tools.logging import log
from ..basics.configurable import Configurable
from .memoryestimator import MemoryEstimator
from ..simulation.skifile import SkiFile
from ..simulation.parallelization import Parallelization

# -----------------------------------------------------------------

def factors(n):

    """
    This function ...
    :param n:
    :return:
    """

    return set(reduce(list.__add__, ([i, n // i] for i in range(1, int(n ** 0.5) + 1) if n % i == 0)))

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

        # The memory estimator
        self.estimator = MemoryEstimator()

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

        # Get additional properties from the ski file
        self.get_properties()

        # 2. Set the parallelization scheme
        self.set_parallelization()

        # 3. Show the parallelization scheme
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

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # If MPI cannot be used
        if not self.config.mpi:

            # Determine the number of cores per node
            cores_per_node = self.config.nsockets * self.config.ncores

            # Determine optimal number of cores used for threading (not more than 12)
            cores = min(cores_per_node, 12)

            #if self.config.hyperthreading: threads = cores * self.config.threads_per_core
            #else: threads = cores

            # Determine number of threads per core
            threads_per_core = self.config.threads_per_core if self.config.hyperthreading else 1

            # Create the parallelization object
            self.parallelization = Parallelization.from_mode("threads", cores, threads_per_core)

        # If MPI can be used
        else:

            # If memory as not been set
            if self.memory is None:

                # Configure the memory estimator
                self.estimator.config.ski = self.ski
                self.estimator.config.input = self.config.input
                self.estimator.config.ncells = self.config.ncells

                # Don't show the memory
                self.estimator.config.show = False

                # Estimate the memory
                self.estimator.run()

                # Get the serial and parallel parts of the simulation's memory
                #serial_memory = self.estimator.serial_memory
                #parallel_memory = self.estimator.parallel_memory

                # Set the memory requirement
                self.memory = self.estimator.memory

            # Get serial and parallel parts of the memory requirement
            serial_memory = self.memory.serial
            parallel_memory = self.memory.parallel

            # Calculate the total memory of one process without data parallelization
            total_memory = serial_memory + parallel_memory

            # The total memory consumption is larger than the node memory
            if total_memory > self.config.memory: # Ms > Mn

                # If multiple nodes can be used
                if self.config.nnodes > 1:

                    # Total memory consumed by the simulation for all processes (=nnodes) combined
                    #memory_simulation_nnode_processes = serial_memory * self.config.nnodes + parallel_memory

                    memory_per_process = serial_memory + parallel_memory / self.config.nnodes

                    # Mn x Nn < Ms?
                    #if self.config.memory * self.config.nnodes < memory_simulation_nnode_processes:
                    if memory_per_process > self.config.memory: # is the memory used per process larger than the memory per node?
                        raise ValueError("Simulation cannot be run: decrease resolution, use system with more memory per node, or use more nodes")

                    else:

                        # Arbitrarily pick a divisor (DIV) of Npps between 4 and 10
                        divisors = factors(self.config.ncores)
                        divisor = random.choice(divisors)

                        # Nt = min(Npps, DIV)
                        nthreads = min(self.config.ncores, divisor)

                        # Np = Nn x Nppn / Nt
                        ppn = self.config.nsockets * self.config.ncores
                        nprocesses = self.config.nnodes * ppn / nthreads

                        total_ncores = self.config.nnodes * self.config.nsockets * self.config.ncores

                        # Nlambda >= 10 x Np
                        nwavelengths = self.ski.nwavelengthsfile(self.config.input) if self.ski.wavelengthsfile() else self.ski.nwavelengths()
                        if nwavelengths >= 10 * nprocesses:

                            # Determine number of threads per core
                            threads_per_core = self.config.threads_per_core if self.config.hyperthreading else 1

                            # data parallelization
                            # Create the parallelization object
                            self.parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core, threads_per_process=nthreads, data_parallel=True)

                        # try again from picking divisor, but now a larger one (less processes)
                        else: pass

                # If only one node can be used
                else: raise ValueError("Simulation cannot be run: decrease resolution, use system with more memory per node, or use more nodes")

            # If more than one simulation memory fits on a node
            else:

                #Np = min(Mn / Ms, Nppn)
                ppn = self.config.nsockets * self.config.ncores
                nprocesses_per_node = int(math.floor(min(self.config.memory / total_memory, ppn)))

                nprocesses = nprocesses_per_node * self.config.nnodes

                #nthreads = ppn / nprocesses_per_node
                ncores_per_process = ppn / nprocesses_per_node

                # Determine number of threads per core
                threads_per_core = self.config.threads_per_core if self.config.hyperthreading else 1

                threads_per_process = threads_per_core * ncores_per_process

                total_ncores = self.config.nnodes * self.config.nsockets * self.config.ncores

                # Nlambda >= 10 * Np?
                nwavelengths = self.ski.nwavelengthsfile(self.config.input) if self.ski.wavelengthsfile() else self.ski.nwavelengths()
                if nwavelengths >= 10 * nprocesses and self.dustlib_dimension == 3:

                    # data parallelization
                    # Create the parallelization object
                    self.parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core, threads_per_process=threads_per_process, data_parallel=True)

                else:

                    # task parallelization
                    self.parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core, threads_per_process=threads_per_process, data_parallel=False)

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
