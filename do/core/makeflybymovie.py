#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.makeflybymovie Create a flyby movie for a SKIRT simulation.
#
# This script creates a flyby movie for the output of a SKIRT simulation, which is saved next to the simulation output.
# Each instrument in the simulation represents a single frame in the movie, in order of appearance in the ski file.
# Use the \c prepareflybymovie script to equip a ski file with instruments that correspond to a given timelime.
#
# Some parameters, such as the frame rate, are hardcoded in the script and can thus be adjusted by hand.
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
from misc.makeflybymovie import makeflybymovie

# -----------------------------------------------------------------

rate = 15                       # movie frame rate
contrast = False                # if true, a contrast curve is applied to each frame (takes quite a while)
from_percentile = 30            # lower percentile used to clip the luminosity values
to_percentile = 100             # upper percentile used to clip the luminosity values

# -----------------------------------------------------------------

print("Starting makeflybymovie...")

# get the command-line argument specifying the simulation(s)
argument = sys.argv[1] if len(sys.argv) > 1 else ""

# construct the list of simulation objects and make the plots
for simulation in createsimulations(argument):
    makeflybymovie(simulation, rate=rate, contrast=contrast,
                   from_percentile=from_percentile, to_percentile=to_percentile)

print("Finished makeflybymovie.")

# -----------------------------------------------------------------
