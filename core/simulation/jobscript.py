#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.jobscript Managing a jobscript for SKIRT
#
# An instance of the JobScript class in this module manages a job script for the STEVIN infrastructure

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import subprocess
import datetime
from distutils.spawn import find_executable

# Import the relevant PTS classes and modules

# -----------------------------------------------------------------

# Define the number of physical cores on the different clusters
cores = {'delcatty': 16, 'gulpin': 32}

# Determine the full path to the SKIRT executable
skirtpath = find_executable("skirt")

# Determine the path to the SKIRT directory
skirtdir = os.path.dirname(os.path.dirname(os.path.dirname(skirtpath)))

# Determine the path to the SKIRT/run directory
rundir = os.path.join(skirtdir, "run")

# -----------------------------------------------------------------

## An instance of the JobScript class in this module manages a job script for the STEVIN infrastructure
#
class JobScript(object):
    
    ## The constructor of the JobScript class takes the following arguments:
    #
    #  - path: the path of the job script to be created
    #  - skifilepath: the path to the ski file that should be run
    #  - nodes: the number of nodes to use
    #  - ppn: the number of desired processors per node
    #  - threadspp: the number of threads (per process), passed directly to the SKIRT executable
    #  - inputpath: the path of the directory to contain the input for the simulation
    #  - outputpath: the path of the directory to contain the output of the simulation
    #  - walltime: an (over)estimate of the required time to complete the simulation
    #  - mail: this flag indicates whether the user wants to receive e-mails when the job is started, completed or aborted.
    #  - verbose: this flag can be used to turn SKIRT's verbose logging mode on or off
    #  - fullnode: this flag can be set to True if one wants to request at least one full node, irrespective of the
    #              number of processors (ppn) that is needed. When threadspp > 1, this is the default behaviour.
    #              If a job is launched with pure mpi (threadspp = 1) where the number of processes is less than the
    #              number of cpu's on a node, these parallel processes could get scattered amongst different nodes,
    #              potentially increasing communication time and being affected by interference of other programs (from
    #              other HPC users). Do not set this flag if you don't care about the reproducibility of your simulation
    #              in terms of computation time.
    #
    def __init__(self, path, parameters, walltime):

        """
        The constructor ...
        :return:
        """

        # Set the file path
        self._path = path

        # Open the job script file
        self._script = open(path, 'w')

        # Write a general header to the job script
        self._script.write("#!/bin/sh\n")
        self._script.write("# Batch script for running SKIRT on a remote system\n")
        self._script.write("#\n")

        # Determine the walltime in "hours, minutes and seconds" format
        m, s = divmod(walltime, 60)
        h, m = divmod(m, 60)

        # Determine an appropriate name for this job
        name = self._skifilename + "_" + str(nodes) + "_" + str(ppn)

        # Set the environment variables
        self._script.write("#PBS -N " + name + "\n")
        self._script.write("#PBS -o output_" + name + ".txt\n")
        self._script.write("#PBS -e error_" + name + ".txt\n")
        self._script.write("#PBS -l walltime=%d:%02d:%02d\n" % (h, m, s))
        self._script.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n")
        #if mail: self._script.write("#PBS -m bae\n")
        self._script.write("#\n")
        self._script.write("\n")

        # Load cluster modules
        self._script.write("# Load the necessary modules\n")
        self._script.write("module load jobs\n")
        self._script.write("module load lxml/3.4.2-intel-2015a-Python-2.7.9\n")
        self._script.write("\n")

        # Run the simulation
        self._script.write("# Run the simulation\n")
        #self._script.write("cd " + self._directorypath + "\n")

        # Write the command string to the job script
        self._script.write(parameters.to_command() + "\n")

    # -----------------------------------------------------------------

    def old__init__(self, path, skifilepath, cluster, nodes, ppn, threadspp, inputpath, outputpath, walltime, mail=False, brief=False, verbose=False, fullnode=False):

        # Set the name of the cluster
        self._clustername = cluster

        # Set the file path
        self._path = path

        # Open the job script file
        self._script = open(path, 'w')

        # The name of the ski file
        self._skifilename = os.path.splitext(os.path.basename(skifilepath))[0]

        # The path of the directory containing the ski file
        self._directorypath = os.path.dirname(skifilepath)

        # Set the input and output paths
        self._inputpath = inputpath
        self._outputpath = outputpath

        # Determine the walltime in "hours, minutes and seconds" format
        m, s = divmod(walltime, 60)
        h, m = divmod(m, 60)

        # Determine an appropriate name for this job
        name = self._skifilename + "_" + str(nodes) + "_" + str(ppn)

        # Check whether we are dealing with multithreading. If so, we calculate the number of processes per
        # node and the requested number of processors per node is set to the maximum (for performance reasons).
        hybrid_processes = 1
        if threadspp > 1:

            # The number of processes per node = [processors per node] / [threads (processors) per process]
            hybrid_processes = ppn / threadspp

            # For hybrid (or threads) mode we always request the full node.
            # Therefore, we determine the number of cores on the node.
            ppn = cores[self._clustername]

        # In MPI mode, we also request a full node for processors < cpu count of a node, if specified by the fullnode flag
        elif fullnode:

            # Set the number of processes per node
            hybrid_processes = ppn

            # Set the requested number of processors on the node to the maximum (a full node)
            ppn = cores[self._clustername]

        # Write a general header to the job script
        self._script.write("#!/bin/sh\n")
        self._script.write("# Batch script for running SKIRT on the UGent HPC infrastructure\n")
        self._script.write("#\n")

        # Set the environment variables
        self._script.write("#PBS -N " + name + "\n")
        self._script.write("#PBS -o output_" + name + ".txt\n")
        self._script.write("#PBS -e error_" + name + ".txt\n")
        self._script.write("#PBS -l walltime=%d:%02d:%02d\n" % (h, m, s))
        self._script.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n")
        if mail: self._script.write("#PBS -m bae\n")
        self._script.write("#\n")
        self._script.write("\n")

        # Load cluster modules
        self._script.write("# Load the necessary modules\n")
        self._script.write("module load jobs\n")
        self._script.write("module load lxml/3.4.2-intel-2015a-Python-2.7.9\n")
        self._script.write("\n")

        # Run the simulation
        self._script.write("# Run the simulation\n")
        self._script.write("cd " + self._directorypath + "\n")

        # Construct a string that represents the SKIRT execution command
        commandstring = "mympirun "

        # Add the appropriate syntax for hybrid / multithreaded runs
        if threadspp > 1 or fullnode: commandstring += "--hybrid " + str(hybrid_processes) + " "

        # Add the number of threads per process and the SKIRT output path to the command string
        commandstring += "skirt -t " + str(threadspp) + " -i " + inputpath + " -o " + outputpath + " "

        # If brief console logging is desired, we pass the '-b' flag to SKIRT
        if brief: commandstring += "-b "

        # If verbose mode is desired, we pass the '-v' flag to SKIRT
        if verbose: commandstring += "-v "

        # Finally, give the name of the ski file to SKIRT
        commandstring += self._skifilename + ".ski"

        # Write the command string to the job script
        self._script.write(commandstring + "\n")

    # -----------------------------------------------------------------

    ## Add an additional command to the job script, optionally preceeded by a comment line
    def addcommand(self, command, comment=""):

        # Add a white line
        self._script.write("\n")

        # Write a comment line preceeding the actual command
        if comment: self._script.write("# " + comment + "\n")

        # Add the command to the job script
        self._script.write(command + "\n")

    # -----------------------------------------------------------------

    ## Submit the script on the cluster
    def submit(self):

        # First, close the file
        self._script.close()

        # Then, swap to the desired cluster and launch the job script
        output = subprocess.check_output("module swap cluster/" + self._clustername + "; qsub " + self._path, shell=True, stderr=subprocess.STDOUT)

        # Read the job ID from the qsub output
        jobid = int(output.split(".")[0])

        # Create the job file
        jobfilepath = os.path.join(rundir, str(jobid) + ".txt")
        jobfile = open(jobfilepath, 'w')

        # Add the information
        jobfile.write("ski file name: " + self._skifilename + "\n")
        jobfile.write("ski file directory: " + self._directorypath + "\n")
        jobfile.write("input directory: " + self._inputpath + "\n")
        jobfile.write("output directory: " + self._outputpath + "\n")
        jobfile.write("submitted at: " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S.%f") + "\n")

        # Close the file
        jobfile.close()

    # -----------------------------------------------------------------

    ## This function removes the job script
    def remove(self):

        # Remove the script file from disk
        os.remove(self._path)

# -----------------------------------------------------------------
