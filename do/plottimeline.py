#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plottimeline Plot a timeline for the different simulation phases of a SKIRT simulation
#

# -----------------------------------------------------------------

# Import standard modules
import os
import argparse

# Import relevant PTS modules
from pts.log import Log
from pts.skirtsimulation import createsimulations
from pts.plottimeline import plottimeline

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Inform the user
log.info("Starting plottimeline...")

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('simulations', type=str, help='a string identifying the simulation(s)', nargs='?', default="")

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
simulations = args.simulations

# Initialize an empty list to contain the paths of timeline files
timelinefiles = []

# Find .dat files that contain extracted timeline information
for root, dirs, filenames in os.walk(os.getcwd()):

    # For each file in the current (sub)directory
    for filename in filenames:

        # Check if this is a timeline data file
        if filename.endswith('.dat') and not filename.startswith(".") and "timeline" in filename:

            # Add this file to the list of progress files
            timelinefiles.append((root, filename))

# If the timelinefiles list is not empty
if timelinefiles:

    # Determine the full path to the visualization directory
    vispath = os.path.join(os.getcwd(), "vis")

    # For each progress file in the list
    for directory, filename in timelinefiles:

        # Determine the full path to the progress file
        timelinefilepath = os.path.join(directory, filename)

        # Determine the path to the directory that will contain the plot from this progress file
        plotpath = os.path.join(vispath, os.path.basename(directory))

        # Create this directory, if it didn't already exist
        try: os.makedirs(plotpath)
        except OSError: pass

        # Determine the path to the plot file
        plotfilepath = os.path.join(plotpath, os.path.splitext(filename)[0] + "_" + ".pdf")

        # Plot the timeline for this simulation
        plottimeline(timelinefilepath, plotfilepath)

# If no extracted timeline information could be found, create a list of simulations from the current directory or
# a string given as a command-line argument and first extract the timeline information for these simulations into
# a temporary timeline data file
else:

    # Initialize an empty list to contain the paths of timeline data files
    progressfiles = []

    # Construct the list of simulation objects and make the plots
    for simulation in createsimulations(simulations):

        # Load the extract function to extract timeline information from simulation log files
        from do.extracttimeline import extract

        # Create a temporaray progress file
        timelinefilepath = os.path.join(os.getcwd(), "timeline_" + simulation.prefix() + ".dat")
        timelinefile = open(timelinefilepath, 'w')

        # Write a header to this new file
        timelinefile.write("# Timeline information for simulation " + simulation.prefix() + "\n")

        # Close the timeline file (results will be appended)
        timelinefile.close()

        # Determine the name of the ski file for this simulation
        skifilename = simulation.prefix() + ".ski"

        # Determine the output path of this simulation
        outputpath = simulation.outpath()

        # Extract the timeline information
        extract(skifilename, outputpath, timelinefilepath)

        # Determine the path to the plot file
        plotfilepath = "timeline_" + simulation.prefix() + ".pdf"

        # Plot the timeline for this simulation
        plottimeline(timelinefilepath, plotfilepath)

# Inform the user
log.info("Finished plottimeline.")

# -----------------------------------------------------------------
