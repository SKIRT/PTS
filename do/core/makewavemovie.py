#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.makewavemovie Create a movie that runs through all wavelengths in the SKIRT simulation output.
#
# This script creates a movie for the output of each SKIRT simulation specified through the command line argument
# (see below). The movie combines the SEDs (bottom panel) and the pixel frames (top panel, from left to right)
# for up to three instruments,  running through all wavelengths in the simulation. The movie is placed next to the
# original file(s) with a similar name (omitting the instrument name) but a different extension.
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.
#
# By default both axes of the SED plot and the luminosity of the frames are autoscaled. You can hardcode specific
# ranges in the script.

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.plot.wavemovie import makewavemovie

# -----------------------------------------------------------------

# a value of None means that the axis is autoscaled;
# alternatively specify a range through a tuple with min and max values
xlim = None
ylim = None
#xlim = ( 5e-2, 1e3 )
#ylim = ( 1e-13, 1e-9 )

# the percentile values, in range [0,100], used to clip the luminosity values
# loaded from the fits files; the default values are 30 and 100 respectively
from_percentile = 30
to_percentile = 100

# -----------------------------------------------------------------

print "Starting makewavemovie..."

# get the command-line argument specifying the simulation(s)
argument = sys.argv[1] if len(sys.argv) > 1 else ""

# construct the list of simulation objects and make the movies
for simulation in createsimulations(argument):
    makewavemovie(simulation, xlim=xlim, ylim=ylim, from_percentile=from_percentile, to_percentile=to_percentile)

print "Finished makewavemovie"

# -----------------------------------------------------------------
