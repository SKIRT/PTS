#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.scaling Test the scaling of SKIRT on a particular system
#

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import argparse
import multiprocessing

# Import the relevant PTS class
from pts.scalingtest import ScalingTest

# -----------------------------------------------------------------

logical_cores = multiprocessing.cpu_count()

# Create the command-line parser and a set of subparsers
parser = argparse.ArgumentParser()
parser.add_argument('system', type=str, help='a name identifying the system')

subparsers = parser.add_subparsers(dest='mode')

# Create the mpi subparser
mpi = subparsers.add_parser('mpi', help="select the MPI mode for the scaling test")
mpi.add_argument('processes', type=int, help='the maximum number of processes')

# Create the threads subparser
threads = subparsers.add_parser('threads', help="select the multithreading mode for the scaling test")
threads.add_argument('threads', type=int, help='the maximum number of threads')

# Create the hybrid subparser
hybrid = subparsers.add_parser('hybrid', help="select the hybrid mode for the scaling test")
hybrid.add_argument('processes', type=int, help='the maximum number of processes')
hybrid.add_argument('threads', type=int, help='the maximum number of threads per process')

args = parser.parse_args()

# Set the command-line options
system = args.system
mode = args.mode
processes = args.processes if hasattr(args, 'processes') else 1
threads = args.threads if hasattr(args, 'threads') else 1

# -----------------------------------------------------------------

# Set the path for the scaling test
scalingname = "SKIRT-SCALING"
scalingpath = os.path.join(os.getenv("HOME"), scalingname)

# -----------------------------------------------------------------

# Run the test
test = ScalingTest(scalingpath, system, mode)
test.run(processes, threads)

# -----------------------------------------------------------------

# Make the plots
test.plot()

# -----------------------------------------------------------------
