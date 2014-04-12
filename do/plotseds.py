#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotseds Plot SEDs for the output of a SKIRT simulation.
#
# This script plots SEDs listed in SKIRT \c prefix_instr_sed.dat output files. All SEDs resulting from
# a particular simulation are plotted on the same axes. The result is saved as a PDF file placed next
# to the original files, and with a similar name constructed as \c prefix_sed.pdf (i.e. leaving out the
# instrument names).
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.
#
# By default both axes are autoscaled; you can hardcode specific ranges in the script.

# -----------------------------------------------------------------

xlim = None                 # default: autoscale x axis
#xlim = ( 5e-2, 1e3 )        # specify range for x axis

ylim = None                 # default: autoscale y axis
#ylim = ( 1e-13, 1e-9 )      # specify range for y axis

# -----------------------------------------------------------------

# import standard modules
import sys
import os.path

# import relevant PTS modules
from pts.visualizer import Visualizer
from pts.skirtsimulation import SkirtSimulation

# -----------------------------------------------------------------

print "Starting plotseds..."

# get the command-line argument specifying the simulation prefix, if any
argument = sys.argv[1] if len(sys.argv) > 1 else ""

# construct the visualizer instance
visualizer = Visualizer(argument, log=True)

# make the plots
visualizer.plotseds(xlim=xlim, ylim=ylim)

print "Finished plotseds."

# -----------------------------------------------------------------
