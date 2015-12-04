#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module is used for managing a jobscript for executing SKIRT in a remote scheduling environment
"""

# -----------------------------------------------------------------

# Import standard modules
import os

# Import the relevant PTS classes and modules

# -----------------------------------------------------------------

class JobScript(object):

    """
    An instance of the JobScript class manages a job script for executing SKIRT in a remote scheduling environment
    """

    def __init__(self, path, arguments, cluster, skirt_path, mpi_command, modules, walltime, nodes, ppn, name=None, mail=False, full_node=False):

        """
        The constructor takes the following arguments:
        :param path: the path of the job script to be created
        :param arguments: the SKIRT command-line arguments
        :param cluster: the cluster configuration settings
        :param skirt_path: the path to the SKIRT executable
        :param mpi_command: the command to be used for starting the multiprocessing environment
        :param walltime: the estimated walltime
        :param nodes: the number of nodes to be used for the simulation
        :param ppn: the requested number of processors per node
        :param full_node: this flag can be set to True if one wants to request at least one full node, irrespective of the
            number of processors (ppn) that is needed. When threadspp > 1, this is the default behaviour.
            If a job is launched with pure mpi (threadspp = 1) where the number of processes is less than the
            number of cpu's on a node, these parallel processes could get scattered amongst different nodes,
            potentially increasing communication time and being affected by interference of other programs (from
            other HPC users). Do not set this flag if you don't care about the reproducibility of your simulation
            in terms of computation time.
        :return:
        """

        # Set the file path
        self.path = path

        # Open the job script file
        self.script = open(path, 'w')

        # Write a general header to the job script
        self.script.write("#!/bin/sh\n")
        self.script.write("# Batch script for running SKIRT on a remote system\n")
        self.script.write("#\n")

        # Determine the walltime in "hours, minutes and seconds" format
        minutes, seconds = divmod(walltime, 60)
        hours, minutes = divmod(minutes, 60)

        # Determine an appropriate name for this job
        if name is None:
            prefix = os.path.basename(arguments.ski_pattern)
            name = prefix + "_" + str(nodes) + "_" + str(ppn)

        # Check whether we are dealing with multithreading. If so, we calculate the number of processes per
        # node and the requested number of processors per node is set to the maximum (for performance reasons).
        hybrid_processes = 1
        if arguments.parallel.threads > 1:

            # The number of processes per node = [processors per node] / [threads (processors) per process]
            hybrid_processes = ppn / arguments.parallel.threads

            # For hybrid (or threads) mode we always request the full node.
            # Therefore, we determine the number of cores on the node.
            ppn = cluster.cores

        # In MPI mode, we also request a full node for processors < cpu count of a node, if specified by the fullnode flag
        elif full_node:

            # Set the number of processes per node
            hybrid_processes = ppn

            # Set the requested number of processors on the node to the maximum (a full node)
            ppn = cluster.cores

        # Set the environment variables
        self.script.write("#PBS -N " + name + "\n")
        self.script.write("#PBS -o output_" + name + ".txt\n")
        self.script.write("#PBS -e error_" + name + ".txt\n")
        self.script.write("#PBS -l walltime=%d:%02d:%02d\n" % (hours, minutes, seconds))
        self.script.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n")
        if mail: self.script.write("#PBS -m bae\n")
        self.script.write("#\n")
        self.script.write("\n")

        # Load cluster modules
        if modules:
            self.script.write("# Load the necessary modules\n")
            #self.script.write("module load jobs\n")
            for module_name in modules:
                self.script.write("module load " + module_name)

        self.script.write("\n")

        # Run the simulation
        self.script.write("# Run the simulation\n")

        # Add the appropriate syntax for hybrid / multithreaded runs
        if arguments.parallel.threads > 1 or full_node: mpi_command += " --hybrid " + str(hybrid_processes) + " "

        # Write the command string to the job script
        self.script.write(arguments.to_command(skirt_path, mpi_command, scheduler=True, to_string=True) + "\n")

        # Close the script file
        self.close()

    # -----------------------------------------------------------------

    def close(self):

        """
        This function ...
        :return:
        """

        # Close the file
        self.script.close()

    # -----------------------------------------------------------------

    def remove(self):

        """
        This function removes the jobscript
        :return:
        """

        # Remove the script file from disk
        if os.path.isfile(self.path): os.remove(self.path)

# -----------------------------------------------------------------
