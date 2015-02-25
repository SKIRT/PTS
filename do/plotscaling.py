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
from pts.plotscaling import plotscaling

# -----------------------------------------------------------------

# Create the command-line parser and a set of subparsers
parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='?', type=str, help='the file for which you want to create the plot', default="")

# Parse the command line arguments
args = parser.parse_args()
filename = args.file

# -----------------------------------------------------------------

# Get the current working directory
directory = os.getcwd()

# Make the scaling plots
# TODO: combine the results of different scaling tests of the same simulation
# TODO: choose multiple results files to include in the plots
plotscaling(directory, filename)

# -----------------------------------------------------------------
