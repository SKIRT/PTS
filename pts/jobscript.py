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
import multiprocessing

# -----------------------------------------------------------------

## An instance of the JobScript class in this module manages a job script for the STEVIN infrastructure
#
class JobScript:
    
    ## The constructor of the JobScript class takes the following arguments:
    #
    #  - path: the path of the job script to be created
    #  - skifilepath: the path to the ski file that should be run
    #  - nodes: the number of nodes to use
    #  - ppn: the number of desired processors per node
    #  - threadspp: the number of threads (per process), passed directly to the SKIRT executable
    #  - outputpath: the path of the directory to contain the output of the simulation
    #  - walltime: an (over)estimate of the required time to complete the simulation
    #
    def __init__(self, path, skifilepath, nodes, ppn, threadspp, outputpath, walltime, mail=False):

        # Save the file path
        self._path = path

        # Open the job script file
        self._script = open(path, 'w')

        # The name of the ski file
        skifilename = os.path.splitext(os.path.basename(skifilepath))[0]

        # The path of the directory containing the ski file
        directorypath = os.path.dirname(skifilepath)

        # Determine the walltime in "hours, minutes and seconds" format
        m, s = divmod(walltime, 60)
        h, m = divmod(m, 60)

        # Check whether we are dealing with multithreading. If so, we calculate the number of processes per
        # node and the requested number of processors per node is set to the maximum (for performance reasons).
        hybrid_processes = 1
        if threadspp > 1:

            # The number of processes per node = [processors per node] / [threads (processors) per process]
            hybrid_processes = ppn / threadspp

            # For hybrid (or threads) mode we always request the full node.
            # Therefore, we determine the number of cores on the node.
            ppn = multiprocessing.cpu_count()

        # Write a general header to the job script
        self._script.write("#!/bin/sh\n")
        self._script.write("# Batch script for running SKIRT on the UGent HPC infrastructure\n")
        self._script.write("#\n")

        # Set the environment variables
        self._script.write("#PBS -N " + skifilename + "_" + str(nodes) + "_" + str(ppn) + "\n")
        self._script.write("#PBS -o output_" + skifilename + "_" + str(nodes) + "_" + str(ppn) + ".txt\n")
        self._script.write("#PBS -e error_" + skifilename + "_" + str(nodes) + "_" + str(ppn) + ".txt\n")
        self._script.write("#PBS -l walltime=%d:%02d:%02d\n" % (h, m, s))
        self._script.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n")
        if mail:
            self._script.write("#PBS -m bae\n")
        self._script.write("#\n")
        self._script.write("\n")

        # Load cluster modules
        self._script.write("# Load the necessary modules\n")
        self._script.write("module load jobs\n")
        self._script.write("module load GCC/4.8.3\n")
        #self._script.write("module load ictce/7.1.2\n")
        self._script.write("module load Python/2.7.3-ictce-4.0.6\n")
        self._script.write("\n")

        # Run the simulation
        self._script.write("# Run the simulation\n")
        self._script.write("cd " + directorypath + "\n")

        hybridoptions = ""
        if threadspp > 1:

            hybridoptions = "--hybrid " + str(hybrid_processes) + " "

        self._script.write("mympirun " + hybridoptions + "skirt -t " + str(threadspp) + " -o " + outputpath + " " + skifilename + ".ski\n")

    ## Add an additional command to the job script, optionally preceeded by a comment line
    def addcommand(self, command, comment=""):

        # Add a white line
        self._script.write("\n")

        # Write a comment line preceeding the actual command
        if comment:
            self._script.write("# " + comment + "\n")

        # Add the command to the job script
        self._script.write(command + "\n")

    ## Submit the script on the cluster
    def submit(self):

        # First, close the file
        self._script.close()

        # Then, launch the job script
        FNULL = open(os.devnull, 'w')   # We ignore the output of the qsub command
        subprocess.call(("qsub",), stdin=open(self._path), stdout=FNULL, stderr=subprocess.STDOUT)

    ## Remove the job script
    def remove(self):

        os.remove(self._path)

# -----------------------------------------------------------------
