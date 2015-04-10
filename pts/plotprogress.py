#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plotprogress Plot progress for the various phases of a SKIRT simulation
#
# The function in this module plots the progress in function of time for certain phases of a SKIRT simulation,
# based on the log messages. A seperate PDF plot is created for each of the following phases:
# - shooting photons for stellar emission;
# - calculating dust emission spectra;
# - shooting photons for dust emission.
#

# -----------------------------------------------------------------

# Import standard modules
import os.path
import numpy as np

# Use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# Import relevant PTS modules
from pts.log import Log

# -----------------------------------------------------------------

# Create a logger
log = Log()

## This function plots the progress in function of time for certain phases of a SKIRT simulation, based on the progress
#  information extracted from its log files. The plots are saved in PDF format and are placed next to the original
#  file(s) with a similar name.
def plotprogress(filepath, plotpath, phase, figsize=(10,6)):

    # Try to get the values from this progress file. If no data could be found in the file, we skip it.
    try:
        ranks, phases, times, percentages = np.loadtxt(filepath, usecols=(0,1,2,3), unpack=True, dtype=str)
    except ValueError:
        log.warning("The file " + filepath + " does not contain any data")
        return

    # Convert the appropriate columns to integers or floats
    ranks = ranks.astype(int, copy=False)
    times = times.astype(float, copy=False)
    percentages = percentages.astype(float, copy=False)

    # Setup the figure
    plt.figure(figsize=figsize)

    # Determine the number of processes
    processes = ranks[-1]+1

    numplots = 0

    # For each process, plot the progress to the figure
    for rank in range(processes):

        # Create lists for the runtimes and percentages for the current process
        times_process = []
        percentages_process = []

        # Loop over the entries of the progress data
        for i in range(len(ranks)):

            # If we are below the current rank, continue to the next entry
            if ranks[i] < rank: continue

            # If we passed the current rank, stop searching for new entries
            elif ranks[i] > rank: break

            # Else, check whether this entry correponds with the specified phase
            elif phases[i] == phase:

                times_process.append(times[i])
                percentages_process.append(percentages[i])

        # Name of the current process
        process = "P" + str(rank)

        # If progress data is available for this process
        if len(percentages_process) > 1:

            # Add the progress of the current process to the figure
            plt.plot(times, percentages, label=process)
            numplots += 1

    # If we actually plotted something, generate the figure and save it
    if numplots > 0:

        plt.xlim(0)
        plt.ylim(0,100)
        plt.grid('on')
        plt.xlabel("time (s)", fontsize='large')
        plt.ylabel("progress (%)", fontsize='large')
        plt.title("Progress of emitting {} photons".format(phase))
        plt.legend(loc='lower right', ncol=4, prop={'size':8})

        # Save the figure
        plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
        log.info("Created PDF progress plot file " + plotpath)

    # Close the plotting environment
    plt.close()

# -----------------------------------------------------------------
