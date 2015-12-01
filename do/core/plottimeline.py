#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plottimeline Plot a timeline for the different simulation phases of a SKIRT simulation
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.core.plot.timeline import TimelinePlotter

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('simulations', type=str, help='a string identifying the simulation(s)', nargs='?', default="")
parser.add_argument('--force', action='store_true', help='force replotting already existing plots')

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
simulations = args.simulations
force = args.force

# Get the current working directory
directory = os.getcwd()

# Create a plotting object
plotter = TimelinePlotter(directory, simulations)

# Plot the timelines
plotter.plot(force)

# -----------------------------------------------------------------
