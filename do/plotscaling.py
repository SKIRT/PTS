#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotscaling Make plots for the scaling benchmark output
#

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import argparse

# Import the relevant PTS class
from pts.plotscaling import ScalingPlotter

# -----------------------------------------------------------------

# Create the command-line parser and a set of subparsers
parser = argparse.ArgumentParser()
parser.add_argument('system', nargs='?', type=str, help='the system for which you want to create the plot', default="")

# Parse the command line arguments
args = parser.parse_args()
system = args.system

# -----------------------------------------------------------------

# Get the current working directory
directory = os.getcwd()

# Create a plotting object
plotter = ScalingPlotter(directory, system)

# Plot the runtimes as a function of the number of threads
plotter.plottimes(xlim=(0,40))

# Plot the speedups as a function of the number of threads
plotter.plotspeedups(xlim=(0,40))

# Plot the efficiencies as a function of the number of threads
plotter.ploteffs(xlim=(0,40))

# -----------------------------------------------------------------
