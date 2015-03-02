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

# Import the relevant PTS class
from pts.scalingtest import ScalingTest

# -----------------------------------------------------------------

# Create the command-line parser and a set of subparsers
parser = argparse.ArgumentParser()
parser.add_argument('system', type=str, help='a name identifying the system')
parser.add_argument('simulation', type=str, help='the name of the simulation to use for the test', nargs='?', default="")
parser.add_argument('mode', type=str, help='the mode for the scaling test', choices=['mpi', 'hybrid', 'threads'])
parser.add_argument('maxnodes', type=float, help='the maximum number of nodes', nargs='?', default=1)
parser.add_argument('minnodes', type=float, help='the minimum number of nodes', nargs='?', default=0)
parser.add_argument('threadspp', type=int, help='if hybrid mode is selected, this argument defines the number of threads per process', nargs='?', default=0)
parser.add_argument('--keep', action='store_true')

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
system = args.system
simulation = args.simulation
mode = args.mode
maxnodes = args.maxnodes
minnodes = args.minnodes
threadspp = args.threadspp
keepoutput = args.keep

# -----------------------------------------------------------------

# Set the path for the scaling test
scalingname = "SKIRTscaling"
scalingpath = os.path.join(os.getenv("HOME"), scalingname)

# -----------------------------------------------------------------

# Run the test
test = ScalingTest(scalingpath, simulation, system, mode, threadspp)
test.run(maxnodes, minnodes, keepoutput)

# -----------------------------------------------------------------
