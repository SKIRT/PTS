#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotpolarization Plot polarization maps for the output of a SKIRT simulation.
#
# This script plots a polarization map based on the SKIRT \c prefix_instr_stokes*.fits output files for a particular
# instrument, placing a PDF file with the name \c prefix_instr_stokes.pdf next to the original set of files.
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.

# -----------------------------------------------------------------

# import standard modules
import sys

# import relevant PTS modules
from pts.skirtsimulation import createsimulations
from pts.plotpolarization import plotpolarization

# -----------------------------------------------------------------

print "Starting plotpolarization..."

# get the command-line argument specifying the simulation(s)
argument = sys.argv[1] if len(sys.argv) > 1 else ""

# construct the list of simulation objects and make the plots
for simulation in createsimulations(argument):
    plotpolarization(simulation)

print "Finished plotpolarization."

# -----------------------------------------------------------------
