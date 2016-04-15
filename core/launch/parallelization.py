#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.parallelization Contains the Parallelization class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class Parallelization(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        # The number of cores
        self.cores = None

        # The number of threads used per core
        self.threads_per_core = None

        # The number of processes
        self.processes = None

    # -----------------------------------------------------------------

    @property
    def threads(self):

        """
        This function ...
        :return:
        """

        return self.cores_per_process * self.threads_per_core

    # -----------------------------------------------------------------

    @property
    def cores_per_process(self):

        """
        This function ...
        :return:
        """

        return self.cores / self.processes

    # -----------------------------------------------------------------

    @classmethod
    def from_processes_and_threads(cls, processes, threads, threads_per_core=1):

        """
        This function ...
        :param processes:
        :param threads:
        :param threads_per_core:
        :return:
        """

        # Create a new class instance
        parallelization = cls()

        # Determine the number of required cores
        cores_per_process = threads / threads_per_core
        cores = cores_per_process * processes

        # Set the properties
        parallelization.cores = cores
        parallelization.threads_per_core = threads_per_core
        parallelization.processes = processes

        # Return the new instance
        return parallelization

    # -----------------------------------------------------------------

    @classmethod
    def from_mode(cls, mode, cores, threads_per_core, threads_per_process=None):

        """
        This function ...
        :param mode:
        :param cores:
        :param threads_per_core:
        :param threads_per_process:
        :return:
        """

        # Create a new class instance
        parallelization = cls()

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

        # Set the properties
        parallelization.cores = cores
        parallelization.threads_per_core = used_threads_per_core
        parallelization.processes = processes

        # Return the new instance
        return parallelization

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self.cores == other.cores and self.threads_per_core == other.threads_per_core and self.processes == other.processes

# -----------------------------------------------------------------
