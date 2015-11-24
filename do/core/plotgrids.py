#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotgrids Plot dust grids for SKIRT output files.
#
# This script plots dust grids for SKIRT \c prefix_ds_gridxx.dat output files, placing PDF files with the same
# name (but with the ".pdf" extension instead of ".dat") next to the original files.
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys

# Import the relevant PTS modules
from pts.simulation import createsimulations
from plotting.grids import plotgrids

# -----------------------------------------------------------------

print("Starting plotgrids...")

# get the command-line argument specifying the simulation(s)
argument = sys.argv[1] if len(sys.argv) > 1 else ""

# construct the list of simulation objects and make the plots
for simulation in createsimulations(argument):
    plotgrids(simulation)

print("Finished plotgrids")

# -----------------------------------------------------------------
