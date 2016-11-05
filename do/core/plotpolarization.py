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
parser.add_argument("-sim", "--simulation", metavar='sim', nargs='+', type=str, default="",
            help='Simulation(s) to be plotted. Default:all')
parser.add_argument("-i", "--instrument", metavar='instr', nargs='+', type=str, default='all',
            help='Instrument(s) to be plotted. Default:all')
parser.add_argument("-w", "--wavelength", default='all',
            help='Wavelength to be plotted, in micron. Default:all')
parser.add_argument("-x", "--binx", type = int, default = 10,
            help = "Binning range in x-direction. Default:10")
parser.add_argument("-y", "--biny", type = int, default = 10,
            help = "Binning range in y-direction. Default:10")
parser.add_argument("-sc", "--scale", nargs=2, type = float, metavar=('d', 's'), default = [None,None],
            help = "Scale of polarization segments. [degree, length]. Default:[None,None] (automatic)")
parser.add_argument("-vr", "--vertRange", nargs=2, type = float, metavar=('min', 'max'), default = [None,None],
            help = "Range of the background plot. [min, max]. Default:[None, None] (automatic)")
parser.add_argument("-ncb", "--noColBar", action="store_true", default = False,
            help="Plot the colorbar(s) separately. Helpful for combining multiple plots, together with --vertRange.")
parser.add_argument("-pay", "--polAvY", action="store_true", default = False,
            help="Plot the polarization degree integrated over y-direction for all x-pixels")
parser.add_argument("-ppdm", "--plotPolDegMap", action="store_true", default = False,
            help="Plot the polarization degree map as a separate picture")
parser.add_argument("-e", "--export", action="store_true", default = False,
            help="Exports data files in addition to the plots.")
circPol_parser = parser.add_mutually_exclusive_group(required=False)
circPol_parser.add_argument("-pc", '--plotCircular', dest='plotCircular', action='store_true',
            help="Plot the circular polarization maps. Default: automatic")
circPol_parser.add_argument("-npc", '--no-plotCircular', dest='plotCircular', action='store_false',
            help="Do not plot the circular polarization maps. Default: automatic")
parser.set_defaults(plotCircular=None)
parser.add_argument("-npl", "--no-plotLinear", action="store_false", default = True,
            help="Do not plot the linear polarization maps.")
args = parser.parse_args()
binsize = (args.binx, args.biny)
# construct the list of simulation objects and make the plots
for simulation in createsimulations(args.simulation):
    print "Starting plotpolarization for simulation '" + simulation.prefix() + "'",
    plotpolarization(simulation, instrumentList=args.instrument, binsize=binsize,
                    wavelength=args.wavelength, polAvY=args.polAvY, export=args.export,
                    degreeLength=args.scale, vertRange=args.vertRange,
                    noColBar=args.noColBar, plotCircular=args.plotCircular,
                    plotPolDegMap=args.plotPolDegMap, plotLinear=args.no_plotLinear)
end = time.time()
print "Finished plotpolarization in {0:0.2f} s".format(end-start)

# -----------------------------------------------------------------
