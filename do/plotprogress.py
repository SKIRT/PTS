#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotprogress Plot progress for the various phases of a SKIRT simulation.
#
# This script plots the progress in function of time for certain phases of a SKIRT simulation, based on the log
# messages. A seperate PDF plot is created for each of the following phases, if present in the simulation:
# - shooting photons for stellar emission ("prefix_progress_stellar_photons.pdf");
# - calculating dust emission spectra ("prefix_progress_dust_spectra.pdf");
# - shooting photons for dust emission ("prefix_progress_dust_photons.pdf");
#
# The dust self-absorption phase, if present, is ignored in the current implementation of the script.
#
# For multi-process (MPI) simulations with verbose logging (i.e. with a separate log file per process),
# the progress for all processes is displayed on the same plot.
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.
#

# -----------------------------------------------------------------

# Import standard modules
import os
import sys
import argparse

# Import relevant PTS modules
from pts.skirtsimulation import createsimulations
from pts.plotprogress import plotprogress

# -----------------------------------------------------------------

print "Starting plotprogress..."

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('simulations', type=str, help='a string identifying the simulation(s)', nargs='?', default="")
parser.add_argument('plotdir', type=str, help='the path to the directory where the plots should be placed', nargs='?', default="")
parser.add_argument('--scaling', action='store_true', help='use this option to construct the plots from the progress data files resulting from a scaling test')

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
simulations = args.simulations
plotdir = args.plotdir if args.plotdir else os.getcwd()
scaling = args.scaling

# Construct the list of simulation objects and make the plots
for simulation in createsimulations(simulations):

    # Plot the progress for this simulation
    plotprogress(simulation, plotdir)

print "Finished plotprogress."

# -----------------------------------------------------------------
