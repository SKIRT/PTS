#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotscaling Make plots for the scaling benchmark output
#

# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import os.path
import argparse

# Import the relevant PTS modules
from plotting.scaling import ScalingPlotter

# *****************************************************************

# The choices for the simulation phase
phases = ['setup', 'stellar', 'dustselfabs', 'dustem', 'writing', 'total']

# *****************************************************************

# Create the command-line parser and a set of subparsers
parser = argparse.ArgumentParser()
parser.add_argument('phase', type=str, help='the simulation phase for which you want to plot the scaling', choices=phases)
parser.add_argument('system', nargs='?', type=str, help='the system for which you want to plot the scaling', default="")
parser.add_argument('--plotfit', action='store_true', help='include the fit to the speedups of a theoretical scaling relation in the plot')

# Parse the command line arguments
args = parser.parse_args()
phase = args.phase
system = args.system
plotfit = args.plotfit

# -----------------------------------------------------------------

# Get the current working directory
directory = os.getcwd()

# Create a plotting object
plotter = ScalingPlotter(directory, phase, system)

# Plot the runtimes as a function of the number of threads
plotter.plottimes()

# Plot the speedups as a function of the number of threads
plotter.plotspeedups(plotfit=plotfit)

# Plot the efficiencies as a function of the number of threads
plotter.ploteffs()

# -----------------------------------------------------------------
