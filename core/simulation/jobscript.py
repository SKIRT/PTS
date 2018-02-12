#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.jobscript Contains the JobScript class, used for managing a jobscript
#  for executing SKIRT in a remote scheduling environment.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from ..remote.jobscript import JobScript as _JobScript
from ..tools import filesystem as fs
from ..tools import numbers

# -----------------------------------------------------------------

class JobScript(object):

    """
    An instance of the JobScript class manages a job script for executing SKIRT in a remote scheduling environment
    """

    def __init__(self, path, arguments, cluster, skirt_path, mpi_command, modules, walltime, nodes, ppn, name=None,
                 mail=False, full_node=False, bind_to_cores=False):

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
        :param bind_to_cores: force process binding to cores
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

            # The number of processes per node = [processors per node] //(integer division) [threads (processors) per process]
            hybrid_processes = ppn // arguments.parallel.threads

            # For hybrid (or threads) mode we always request the full node.
            # Therefore, we determine the number of cores on the node.
            ppn = cluster.cores_per_socket * cluster.sockets_per_node

        # In MPI mode, we also request a full node for processors < cpu count of a node, if specified by the fullnode flag
        elif full_node:

            # Set the number of processes per node
            hybrid_processes = ppn

            # Set the requested number of processors on the node to the maximum (a full node)
            ppn = cluster.cores_per_socket * cluster.sockets_per_node

        # Determine the paths to the output and error files
        output_file_path = os.path.join(arguments.output_path, "output_" + name + ".txt")
        error_file_path = os.path.join(arguments.output_path, "error_" + name + ".txt")

        # Set the environment variables
        self.script.write("#PBS -N " + name + "\n")
        self.script.write("#PBS -o " + output_file_path + "\n")
        self.script.write("#PBS -e " + error_file_path + "\n")
        self.script.write("#PBS -l walltime=%d:%02d:%02d\n" % (hours, minutes, seconds))
        self.script.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n")
        if mail: self.script.write("#PBS -m bae\n")
        self.script.write("#\n")
        self.script.write("\n")

        # Load cluster modules if specified
        if modules:

            self.script.write("# Load the necessary modules\n")
            for module_name in modules: self.script.write("module load " + module_name + "\n")

        # Add whiteline
        self.script.write("\n")

        # Run the simulation
        self.script.write("# Run the simulation\n")

        # Add the appropriate syntax for hybrid / multithreaded runs
        if arguments.parallel.threads > 1 or full_node: mpi_command += " --hybrid " + str(hybrid_processes)

        # Write the command string to the job script
        command = arguments.to_command(scheduler=True, skirt_path=skirt_path, mpirun_path=mpi_command, bind_to_cores=bind_to_cores, to_string=True)
        self.script.write(command + "\n")

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

class MultiJobScript(object):

    """
    This class ...
    """

    def __init__(self, path, arguments, cluster, skirt_path, mpi_command, modules, walltime, nodes, ppn, name=None,
                 mail=False, full_node=False, bind_to_cores=False):

        """
        The constructor ...
        :param path:
        :param arguments: IS A LIST OF SKIRTARGUMENT OBJECTS
        :param cluster:
        :param skirt_path:
        :param mpi_command:
        :param modules:
        :param walltime:
        :param nodes:
        :param ppn:
        :param name:
        :param mail:
        :param full_node:
        :param bind_to_cores:
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

            # The number of processes per node = [processors per node] //(integer division) [threads (processors) per process]
            hybrid_processes = ppn // arguments.parallel.threads

            # For hybrid (or threads) mode we always request the full node.
            # Therefore, we determine the number of cores on the node.
            ppn = cluster.cores_per_socket * cluster.sockets_per_node

        # In MPI mode, we also request a full node for processors < cpu count of a node, if specified by the fullnode flag
        elif full_node:

            # Set the number of processes per node
            hybrid_processes = ppn

            # Set the requested number of processors on the node to the maximum (a full node)
            ppn = cluster.cores_per_socket * cluster.sockets_per_node

        # Determine the paths to the output and error files
        output_file_path = os.path.join(arguments.output_path, "output_" + name + ".txt")
        error_file_path = os.path.join(arguments.output_path, "error_" + name + ".txt")

        # Set the environment variables
        self.script.write("#PBS -N " + name + "\n")
        self.script.write("#PBS -o " + output_file_path + "\n")
        self.script.write("#PBS -e " + error_file_path + "\n")
        self.script.write("#PBS -l walltime=%d:%02d:%02d\n" % (hours, minutes, seconds))
        self.script.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n")
        if mail: self.script.write("#PBS -m bae\n")
        self.script.write("#\n")
        self.script.write("\n")

        # Load cluster modules if specified
        if modules:

            self.script.write("# Load the necessary modules\n")
            for module_name in modules: self.script.write("module load " + module_name + "\n")

        # Add whiteline
        self.script.write("\n")

        # Run the simulation
        self.script.write("# Run the simulation\n")

        # Add the appropriate syntax for hybrid / multithreaded runs
        if arguments.parallel.threads > 1 or full_node: mpi_command += " --hybrid " + str(hybrid_processes)

        # Write the command string to the job script
        command = arguments.to_command(scheduler=True, skirt_path=skirt_path, mpirun_path=mpi_command, bind_to_cores=bind_to_cores, to_string=True)
        self.script.write(command + "\n")

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

class SKIRTJobScript(_JobScript):

    """
    This class
    """

    def __init__(self, name, arguments, cluster, skirt_path, mpi_command, walltime, modules, mail=False,
                 bind_to_cores=False, extra_header_lines=None, hyperthreading=None):

        """
        The constructor ...
        :param name:
        :param arguments:
        :param cluster:
        :param skirt_path:
        :param mpi_command:
        :param walltime:
        :param mail:
        :param bind_to_cores:
        :param extra_header_lines:
        :param hyperthreading:
        """

        # Determine the paths to the output and error files
        output_file_path = fs.join(arguments.output_path, "output_" + name + ".txt")
        error_file_path = fs.join(arguments.output_path, "error_" + name + ".txt")

        cores_per_node = cluster.cores_per_socket * cluster.sockets_per_node

        # Set number of threads per core
        #if hyperthreading is None: hyperthreading = cluster.threads_per_core > 1

        # Set number of threads per core based on hyperthreading flag
        #if hyperthreading: threads_per_core = cluster.threads_per_core
        #else: threads_per_core = 1

        threads_per_core = arguments.parallel.threads_per_core

        # Determine number of threads per node
        threads_per_node = threads_per_core * cores_per_node

        # Determine number of processors per process
        processors_per_process = arguments.parallel.threads / threads_per_core
        assert int(processors_per_process) == processors_per_process
        processors_per_process = int(processors_per_process)

        # HYBRID
        if arguments.parallel.processes > 1 and arguments.parallel.threads > 1:

            processes_per_node = threads_per_node // arguments.parallel.threads

            if processes_per_node == 0: raise RuntimeError("Impossible parallelization scheme: too many threads to fit on a node")

            processors_per_node = processes_per_node * processors_per_process
            processors = processors_per_process * arguments.parallel.processes
            nodes, ppn = get_requirements(processors, cores_per_node, full_node=True)

        # Multiprocessing
        elif arguments.parallel.processes > 1 and arguments.parallel.threads == 1:

            if arguments.parallel.processes > cores_per_node:

                # Determine necessary number of nodes
                nodes, ppn = get_requirements(arguments.parallel.processes, cores_per_node, full_node=True)

                # Divide the processes over the nodes
                processes_per_node = arguments.parallel.processes / nodes
                if not numbers.is_integer(processes_per_node): raise RuntimeError("Impossible parallelization scheme: " + str(arguments.parallel.processes) + " cannot be distributed over " + str(nodes) + " nodes")
                processes_per_node = int(processes_per_node)

            else:

                processes_per_node = cores_per_node
                nodes = 1
                ppn = cores_per_node

        # Multithreading
        elif arguments.parallel.threads > 1 and arguments.parallel.processes == 1:

            if arguments.parallel.threads > threads_per_node: raise RuntimeError("Impossible parallelization schem: too many threads to fit on a node")
            else:

                processes_per_node = 1
                nodes = 1
                ppn = cores_per_node

        # Singlethreading (nthreads = nprocesses = 1)
        else:

            processes_per_node = 1
            nodes = 1
            ppn = cores_per_node

        # Call the constructor of the base class
        super(SKIRTJobScript, self).__init__(name, walltime, nodes, ppn, output_file_path, error_file_path, mail, extra_header_lines=extra_header_lines)

        # Add the appropriate syntax for hybrid / multithreaded runs
        mpi_command += " --hybrid " + str(processes_per_node)

        # Write the command string to the job script
        command = arguments.to_command(scheduler=True, skirt_path=skirt_path, mpirun_path=mpi_command, bind_to_cores=bind_to_cores, to_string=True, report_bindings=False)

        # Add the SKIRT command
        self.add_command(command, "Launch SKIRT")

        # Add modules to load
        for module in modules: self.import_module(module)

# -----------------------------------------------------------------

def get_requirements(processors, cores_per_node, full_node=False):

    """
    This function calculates the required amount of nodes and processors per node, given a certain number of
    processors.
    :param processors:
    :param cores_per_node:
    :param full_node:
    :return:
    """

    # Calculate the necessary amount of nodes
    nodes = processors // cores_per_node + (processors % cores_per_node > 0)

    # Determine the number of processors per node
    ppn = processors if nodes == 1 else cores_per_node

    # Always use full nodes if requested
    if full_node: ppn = cores_per_node

    # Return the number of nodes and processors per node
    return nodes, ppn

# -----------------------------------------------------------------
