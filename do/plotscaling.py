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

# The choices for the simulation phase
phases = ['setup', 'stellar', 'dustselfabs', 'dustem', 'writing', 'total']

# -----------------------------------------------------------------

# Create the command-line parser and a set of subparsers
parser = argparse.ArgumentParser()
parser.add_argument('phase', type=str, help='the simulation phase for which you want to plot the scaling', choices=phases)
parser.add_argument('system', nargs='?', type=str, help='the system for which you want to plot the scaling', default="")
parser.add_argument('--fit', action='store_true', help='fit a theoretical relation to the speedups and save the parameters')

# Parse the command line arguments
args = parser.parse_args()
phase = args.phase
system = args.system
fit = args.fit

# -----------------------------------------------------------------

# Get the current working directory
directory = os.getcwd()

# Create a plotting object
plotter = ScalingPlotter(directory, phase, system)

# Plot the runtimes as a function of the number of threads
plotter.plottimes()

# Plot the speedups as a function of the number of threads
plotter.plotspeedups(fit=fit)

# Plot the efficiencies as a function of the number of threads
plotter.ploteffs()

# -----------------------------------------------------------------
