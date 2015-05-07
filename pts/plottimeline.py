#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plottimeline Plot a timeline of a simulation, for the different parallel processes, from
#  a text file containing the timeline information extracted from the simulation's log files.
#

# -----------------------------------------------------------------

# Import standard modules
import numpy as np

# Use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# Import relevant PTS modules
from pts.log import Log

# -----------------------------------------------------------------

# Ignore warnings, otherwise Canopy would give a UserWarning on top of the error encountered when a progress
# file does not contain any data (an error which is catched an produces an error message).
import warnings
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------

# Define the indices used to identify the different simulation phases
phaseindices = {'setup': 0, 'stellar': 1, 'comm': 2, 'spectra': 3, 'dust': 4, 'writing': 5, 'waiting': 6}

# Define the names of the different simulation phases (for the plot's legend)
phasenames = ['Setup', 'Stellar emission', 'Communication', 'Dust spectra calculation', 'Dust emission', 'Writing', 'Waiting']

# Define the colors for the different simulation phases in the plot
# (setup = red, stellar = green, comm = blue, spectra = magenta, dust = cyan, writing = yellow, waiting = black)
colors = list('rgbmcyk')

# -----------------------------------------------------------------

# Create a logger
log = Log()

## This function plots the timeline of a simulation, for the different parallel processes, from a text file
#  containing the timeline information extracted from the simulation's log files.
def plottimeline(filepath, plotpath, figsize=(12,8), percentages=False):

    # Try to get the values from the txt file. If no data could be found in the file, we skip it.
    try:
        ranks, phases, starttimes, endtimes = np.loadtxt(filepath, usecols=(0,1,2,3), unpack=True)
    except ValueError:
        log.warning("The file " + filepath + " does not contain any data")
        return

    # Determine the number of processors
    nprocs = int(ranks[-1]) + 1
    procranks = range(nprocs)

    # Initialize a data structure to contain the start times and endtimes for the different processes, indexed on the phase
    data = []

    # Iterate over the different rows found in the txt file and copy the start- and endtimes in the data list
    for i in range(len(ranks)):

        if int(ranks[i]) == 0:

            data.append([int(phases[i]), [], []])
            data[len(data)-1][1].append(starttimes[i])
            data[len(data)-1][2].append(endtimes[i])

        else:

            nphases = len(data)
            data[i%nphases][1].append(starttimes[i])
            data[i%nphases][2].append(endtimes[i])

    # Initialize a figure with the appropriate size
    plt.figure(figsize=figsize)
    ax = plt.gca()

    legendEntries = []
    legendNames = []
    unique_phases = []

    # Make the timeline plot, consisting of a set of bars of the same color for each simulation phase
    for phase, starttimes, endtimes in data:

        durations = np.array(endtimes) - np.array(starttimes)
        patch_handle = ax.barh(procranks, durations, color=colors[phase], align='center', left=starttimes, alpha=0.8)

        if phase not in unique_phases and not (phase == phaseindices['comm'] and nprocs == 1):

            unique_phases.append(phase)
            legendEntries.append(patch_handle)
            legendNames.append(phasenames[phase])

        if percentages:

            for rank, patch in enumerate(patch_handle.get_children()):

                bl = patch.get_xy()
                x = 0.5*patch.get_width() + bl[0]
                y = 0.5*patch.get_height() + bl[1]
                ax.text(x, y, "%d%%" % (rank), ha='center')

    # Format the axis ticks and labels
    ax.set_yticks(procranks)
    ax.set_yticklabels(procranks)
    ax.set_xlabel('Time (s)', fontsize='large')
    ax.set_ylabel('Process rank', fontsize='large')
    ax.yaxis.grid(True)

    if nprocs == 1:

        ax.set_frame_on(False)
        fig = plt.gcf()
        fig.set_size_inches(10,3)
        ax.xaxis.tick_bottom()
        ax.yaxis.set_visible(False)

    # Shrink current axis's height by 20% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])

    # Set the plot title
    plt.title("Timeline of the different simulation phases")

    # Put a legend below current axis
    ax.legend(legendEntries, legendNames, loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=False, ncol=3)

    # Save the figure
    plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.40)
    log.info("Created PDF timeline " + plotpath)

    # Close the plotting environment
    plt.close()

# -----------------------------------------------------------------
