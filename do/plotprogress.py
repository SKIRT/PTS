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
import argparse

# Import relevant PTS modules
from pts.log import Log
from pts.skirtsimulation import createsimulations
from pts.plotprogress import plotprogress

# -----------------------------------------------------------------

# The choices for the simulation phase
phases = ['stellar', 'dust', 'spectra']

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Inform the user
log.info("Starting plotprogress...")

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('simulations', type=str, help='a string identifying the simulation(s)', nargs='?', default="")
parser.add_argument('phase', type=str, help='the simulation phase for which you want to plot the progress', choices=phases)

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
simulations = args.simulations
phase = args.phase

# Initialize an empty list to contain the paths of progress files
progressfiles = []

# Find .dat files that contain extracted progress information
for root, dirs, filenames in os.walk(os.getcwd()):

    # For each file in the current (sub)directory
    for filename in filenames:

        # Check if this is a progress data file
        if filename.endswith('.dat') and not filename.startswith(".") and "progress" in filename:

            # Add this file to the list of progress files
            progressfiles.append((root, filename))

# If the progressfiles list is not empty
if progressfiles:

    # Determine the full path to the visualization directory
    vispath = os.path.join(os.getcwd(), "vis")

    # For each progress file in the list
    for directory, filename in progressfiles:

        # Determine the full path to the progress file
        progressfilepath = os.path.join(directory, filename)

        # Determine the path to the directory that will contain the plot from this progress file
        plotpath = os.path.join(vispath, os.path.basename(directory))

        # Create this directory, if it didn't already exist
        try: os.makedirs(plotpath)
        except OSError: pass

        # Determine the path to the plot file
        plotfilepath = os.path.join(plotpath, os.path.splitext(filename)[0] + "_" + phase + ".pdf")

        # Plot the progress information in this file
        plotprogress(progressfilepath, plotfilepath, phase)

# If no extracted progress information could be found, create a list of simulations from the current directory or
# a string given as a command-line argument and first extract the progress for these simulations into a temporary
# progress data file
else:

    # Initialize an empty list to contain the paths of progress files
    progressfiles = []

    # Construct the list of simulation objects and make the plots
    for simulation in createsimulations(simulations):

        # Load the extract function to extract progress information from simulation log files
        from do.extractprogress import extract

        # Create a temporaray progress file
        progressfilepath = os.path.join(os.getcwd(), "progress_" + simulation.prefix() + ".dat")
        progressfile = open(progressfilepath, 'w')

        # Write a header to this new file which contains some general info about its contents
        progressfile.write("# Progress results for simulation " + simulation.prefix() + "\n")
        progressfile.write("# Column 1: Process rank\n")
        progressfile.write("# Column 2: Simulation phase (stellar, dust or spectra)\n")
        progressfile.write("# Column 3: Execution time (s)\n")
        progressfile.write("# Column 4: Progress (%)\n")

        # Close the progress file (results will be appended)
        progressfile.close()

        # Determine the name of the ski file for this simulation
        skifilename = simulation.prefix() + ".ski"

        # Determine the output path of this simulation
        outputpath = simulation.outpath()

        # Extract the progress information
        extract(skifilename, outputpath, progressfilepath)

        # Determine the path to the plot file
        plotfilepath = "progress_" + simulation.prefix() + "_" + phase + ".pdf"

        # Plot the progress for this simulation
        plotprogress(progressfilepath, plotfilepath, phase)

# Inform the user
log.info("Finished plotprogress.")

# -----------------------------------------------------------------
