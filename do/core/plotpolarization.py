#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plotpolarization Plot polarization maps for the output of a SKIRT simulation.
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

# Import standard modules
import sys
import argparse
import time

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.plot.polarization import plotpolarization

# -----------------------------------------------------------------

start = time.time()

#get the command-line argument(s)
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--simulation", metavar='sim', nargs='+', type=str, default="",
            help='Simulation(s) to be plotted. Default:all')
parser.add_argument("-i", "--instrument", metavar='instr', nargs='+', type=str, default='all',
            help='NOT YET IMPLEMENTED: Intstrument(s) to be plotted. Default:all')
parser.add_argument("-w", "--wavelength", default='all',
            help='Wavelenght to be plotted. Default:all')
parser.add_argument("-x", "--binx", type = int, default = 10,
            help = "Binning range in x-direction. Default:10")
parser.add_argument("-y", "--biny", type = int, default = 10,
            help = "Binning range in y-direction. Default:10")
parser.add_argument("-nm", "--noMaps", action="store_true", default = False,
            help="Do not plot the polarization maps, just the additional options (speeds things up).")
parser.add_argument("-pay", "--polAvY", action="store_true", default = False,
            help="Plot the polarization degree integrated over y-direction for all x-pixels")
parser.add_argument("-e", "--export", action="store_true", default = False,
            help="Exports data files in addition to the plots.")
args = parser.parse_args()
binsize = (args.binx, args.biny)
# construct the list of simulation objects and make the plots
for simulation in createsimulations(args.simulation):
    print "Starting plotpolarization for simulation '" + simulation.prefix() + "'",
    plotpolarization(simulation, instrumentList=args.instrument, binsize=binsize,
                    wavelength=args.wavelength, polAvY=args.polAvY, noMaps = args.noMaps,
                    export = args.export)
end = time.time()
print "Finished plotpolarization in {0:0.2f} s".format(end-start)

# -----------------------------------------------------------------
