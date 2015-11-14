#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package plotting.plotprogress

# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS modules
from pts.skirtsimulation import createsimulations
from performance.plotprogress import plotprogress

# *****************************************************************

class ProgressPlotter(object):
    
    """
    This class ...
    """
    
    def __init__(self):
        
        """
        This function ...
        """

    # *****************************************************************

    def run(self, simulations, phase):
        
        """
        This function ...
        """
        
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
                resdirname = os.path.basename(os.path.dirname(directory))
                plotpath = os.path.join(vispath, os.path.basename(directory)) if resdirname == "res" else os.getcwd()

                # Create this directory, if it didn't already exist
                try: os.makedirs(plotpath)
                except OSError: pass

                # Determine the path to the plot file
                plotfilepath = os.path.join(plotpath, os.path.splitext(filename)[0] + "_" + phase + ".pdf")

                # If the file already exists, skip the plotting procedure
                if os.path.isfile(plotfilepath): continue

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
                from ..extract.extractprogress import extract

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
        
# *****************************************************************

