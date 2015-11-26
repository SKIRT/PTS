#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.submit Submit a job for running a ski file to the HPC infrastructure
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import os.path
import argparse

# Import the relevant PTS classes and modules
from pts.core.basics import Log
from pts.core.simulation import SkiFile
from pts.core.simulation import JobScript

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('skifile', type=str, help='the name of the ski file')
parser.add_argument('inpath', type=str, help='the simulation input path')
parser.add_argument('outpath', type=str, help='the simulation output path')
parser.add_argument('nodes', type=int, help='the number of nodes to use for the simulation')
parser.add_argument('walltime', type=str, help='the expected walltime for this simulation in hh:mm:ss format')
parser.add_argument('--threadspp', nargs='?', type=int, help='the number of parallel threads per process', default=0)
parser.add_argument('--cluster', nargs='?', type=str, help='the cluster to which the job should be submitted', default="delcatty")
parser.add_argument('--verbose', action='store_true', help='add this option to enable verbose logging mode')

# Parse the command line arguments
args = parser.parse_args()
skifile = args.skifile
inpath = args.inpath
outpath = args.outpath
nodes = args.nodes
threadspp = args.threadspp
walltime = args.walltime
cluster = args.cluster
verbose = args.verbose

# -----------------------------------------------------------------

# Define the number of cores per node for this different clusters
cores = {'delcatty': 16, 'gulpin': 32}

# Define the available memory for the different clusters (in GB)
memory = {'delcatty': 64, 'gulpin': 64}

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Calculate the walltime in seconds
hours, minutes, seconds = walltime.split(':')
time = int(hours)*3600 + int(minutes)*60 + int(seconds)

# Determine the full path to the ski file
skifilepath = os.path.abspath(skifile)

# Determine the full path to the input and output directories
inputpath = os.path.abspath(inpath)
outputpath = os.path.abspath(outpath)

# Get the number of processors per node on the specified cluster
ppn = cores[cluster]

# If the number of threads per process is not specified, calculate the optimal value
if not threadspp:

    # Calculate the amount of required memory for this simulation (in gigabytes)
    skifile = SkiFile(skifilepath)
    required = estimate_memory(skifile)

    # Determine the number of processes per node and the number of threads per process
    processes = int(memory[cluster]/required)

    if processes < 1:

        log.error("Not enough memory available to run this simulation")
        exit()

    # Calculate the number of threads per process
    threadspp = int(cores[cluster] / processes)

# Create a job script for this simulation
jobscript = JobScript("job.sh", skifilepath, cluster, nodes, ppn, threadspp, inputpath, outputpath, time, mail=True, brief=True, verbose=verbose)

# Submit the job script to the cluster
jobscript.submit()

# Remove the job script
jobscript.remove()

# -----------------------------------------------------------------
