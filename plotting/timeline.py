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

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os.path
import numpy as np
from collections import defaultdict

# Use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# Import relevant PTS modules
from ..pts.log import Log
from ..pts.skirtsimulation import createsimulations

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
# (setup = red, stellar = green, comm = orange, spectra = magenta, dust = cyan, writing = yellow, waiting = blue)
colors = ['r', 'g', '#FF7626', 'm', 'c', 'y', 'b']

# -----------------------------------------------------------------

## An instance of the TimelinePlotter is used to create timeline diagrams for the different simulation phases
#
class TimelinePlotter(object):

    ## The constructor accepts the following arguments:
    #
    #  - directory: the path of the directory where the result files are located and where the plots should be created
    #  - system: the system for which to plot the statistics. If no system is given, plots are created with different
    #            curves for all the different systems (and modes) for which a results file is found.
    #
    def __init__(self, directory, simulations):

        # Create a logger
        self._log = Log()

        # Set the results and visualization paths
        self._respath = os.path.join(directory, "res")
        self._vispath = os.path.join(directory, "vis")

        # Initialize an empty list to contain the paths of timeline files
        self._filepaths = []

        # Find .dat files that contain extracted timeline information
        for root, dirs, filenames in os.walk(directory):

            # For each file in the current (sub)directory
            for filename in filenames:

                # Check if this is a timeline data file
                if filename.endswith('.dat') and not filename.startswith(".") and "timeline" in filename:

                    # Determine the full path to the timeline data file
                    timelinefilepath = os.path.join(root, filename)

                    # Add this file to the list of timeline files
                    self._filepaths.append(timelinefilepath)

        # If no extracted timeline information could be found, create a list of simulations from the current directory or
        # a string given as a command-line argument and first extract the timeline information for these simulations into
        # a temporary timeline data file
        if not self._filepaths:

            # Construct the list of simulation objects and create the timeline data files
            for simulation in createsimulations(simulations):

                # Load the extract function to extract timeline information from simulation log files
                from ..extract.timeline import extract

                # Create a temporaray timeline file
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

                # Add this file to the list of timeline files
                self._filepaths.append(timelinefilepath)

    ## This function should be called to invoke the plotting routine, which will create a timeline plot for each
    #  file that is found with valid timeline data
    def plot(self, force=False):

        # Set the force parameter
        self._force = force

        # Create a data structure to contain the ...
        cpudata = dict()

        # Keep track of the number of processes
        nprocs_dict = defaultdict(list)

        # keep track of which scaling test was weak
        weaktests = dict()

        # Loop over all the timeline data files
        for filepath in self._filepaths:

            # Try to get the values from the data file. If no data could be found in the file, we skip it.
            try:
                ranks, phases, starttimes, endtimes = np.loadtxt(filepath, usecols=(0,1,2,3), unpack=True)
            except ValueError:
                self._log.warning("The file " + filepath + " does not contain any data")
                continue

            # Determine the number of processors and add it to the list
            nprocs = int(ranks[-1]) + 1

            # Make a list of the different process ranks
            procranks = range(nprocs)

            # Initialize a data structure to contain the start times and endtimes for the different processes,
            # indexed on the phase
            data = []

            # Iterate over the different rows found in the data file and copy the start- and endtimes in the data list
            for i in range(len(ranks)):

                if int(ranks[i]) == 0:

                    data.append([int(phases[i]), [], []])
                    data[len(data)-1][1].append(starttimes[i])
                    data[len(data)-1][2].append(endtimes[i])

                else:

                    nphases = len(data)
                    data[i%nphases][1].append(starttimes[i])
                    data[i%nphases][2].append(endtimes[i])

            # Try to obtain the simulation name from the data file, if this works we use it for the name of
            # the timeline plot file
            with open(filepath) as datafile: first = datafile.readlines()[0]
            try:

                # Find the simulation name in the first line of the data file
                simulationname = first.split("simulation ")[1].split("\n")[0]

                # Set the path of the plot file
                plotfilepath = "timeline_" + simulationname + ".pdf"

            # If the simulation name is not in the header of the data file, this data file was created as part
            # of a scaling benchmark test, and the data files are probably under the results path defined by
            # 'self._respath'. The plot files should be placed under the visualisation path ('self._vispath')
            except IndexError:

                # Determine the path of the scaling run directory where this timeline file belongs to
                scalingrunpath = os.path.dirname(filepath)

                # Determine the name of the scaling run
                scalingrunname = os.path.basename(scalingrunpath)

                # Open the info file corresponding with this scaling run
                infofilepath = os.path.join(scalingrunpath, "info.txt")
                with open(infofilepath) as infofile:

                    # Find out whether the file was generated as part of a weak or a strong scaling test
                    second = infofile.readlines()[1]
                    weak = "weak" in second

                # Determine the path of the plotting directory (a subdirectory of the visualization directory)
                plotpath = os.path.join(self._vispath, scalingrunname)

                # Create this directory, if it didn't already exist
                try: os.makedirs(plotpath)
                except OSError: pass

                # Set the path of the plot file
                filename = os.path.basename(filepath)
                plotfilepath = os.path.join(plotpath, os.path.splitext(filename)[0] + ".pdf")

                # Give the cpudata list the correct dimension (the number of different phases encountered in the data)
                if not scalingrunname in cpudata:

                    cpudata[scalingrunname] = [[entry[0], [], []] for entry in data]
                    weaktests[scalingrunname] = weak

                # TODO: FIX THIS FOR ONE PROCESS
                if not nprocs == 1:

                    # Add this number of processes to the list for this scaling test run
                    nprocs_dict[scalingrunname].append(nprocs)

                    # For each seperate phase
                    for i in range(len(data)):

                        starttimes = data[i][1]
                        endtimes = data[i][2]

                        # If the scaling test was weak, add the times for the process with rank zero
                        if weak:

                            cpudata[scalingrunname][i][1].append(starttimes[0])
                            cpudata[scalingrunname][i][2].append(endtimes[0])

                        # Else, if the scaling test was strong, sum the runtimes of each process to obtain
                        # the total CPU time
                        else:

                            if i == 0:

                                starttime = 0.0

                            else:

                                starttime = cpudata[scalingrunname][i-1][2][-1]

                            durations = np.array(endtimes) - np.array(starttimes)
                            totalduration = np.sum(durations)

                            cpudata[scalingrunname][i][1].append(starttime)
                            cpudata[scalingrunname][i][2].append(starttime+totalduration)

            # Create a plot for this timeline data file
            self._createplot(data, plotfilepath, procranks)

        # Plot a graph comparing the timeline for the runs with a different number of processes
        for scalingtest, data in cpudata.items():

            plotpath = os.path.join(self._vispath, scalingtest)
            plotfilepath = os.path.join(plotpath, "timeline.pdf")

            # Check whether this scaling test was weak
            weak = weaktests[scalingtest]

            # Get the list of number of processes for this scaling test
            nprocs_list = nprocs_dict[scalingtest]

            # Create the plot
            self._createplot(data, plotfilepath, nprocs_list, percentages=True, totals=True, unordered=True, numberofproc=True, cpu=(not weak))

    ## This function actually plots the timeline based on a data structure containing the starttimes and endtimes
    #  for the different simulation phases
    def _createplot(self, data, plotfilepath, procranks, figsize=(12,8), percentages=False, totals=False, unordered=False, numberofproc=False, cpu=False):

        # If the file already exists, skip the plotting procedure
        if os.path.isfile(plotfilepath) and not self._force: return

        # Initialize a figure with the appropriate size
        plt.figure(figsize=figsize)
        ax = plt.gca()

        legendEntries = []
        legendNames = []
        unique_phases = []

        # Determine the number of processes
        nprocs = len(procranks)

        # Get the ordering
        if unordered:
            yticks = np.array(procranks).argsort().argsort()
        else:
            yticks = procranks

        durations_list = []
        totaldurations = np.zeros(nprocs)
        patch_handles = []

        # Make the timeline plot, consisting of a set of bars of the same color for each simulation phase
        for phase, starttimes, endtimes in data:

            durations = np.array(endtimes) - np.array(starttimes)
            durations_list.append(durations)

            totaldurations += durations

            patch_handle = ax.barh(yticks, durations, color=colors[phase], align='center', left=starttimes, alpha=0.8, lw=0)
            patch_handles.append(patch_handle)

            if phase not in unique_phases and not (phase == phaseindices['comm'] and nprocs == 1):

                unique_phases.append(phase)
                legendEntries.append(patch_handle)
                legendNames.append(phasenames[phase])

        if percentages:

            # For the different phases
            for phase, patch_handle in enumerate(patch_handles):

                durations = durations_list[phase]

                for sorting_number, rectangle in enumerate(patch_handle.get_children()):

                    duration = durations[sorting_number]
                    percentage = float(duration) / float(totaldurations[sorting_number]) * 100.0

                    x = 0.5*rectangle.get_width() + rectangle.get_x()
                    y = 0.5*rectangle.get_height() + rectangle.get_y()

                    if rectangle.get_width() > 2000:

                        plt.text(x, y, "%d%%" % (percentage), ha='center', va='center', fontsize=10)

        if totals:

            for sorting_number, rectangle in enumerate(patch_handles[-1].get_children()):

                width = rectangle.get_width()
                label_text = str(int(totaldurations[sorting_number]))
                plt.text(rectangle.get_x() + width + 0.02*rectangle.get_x(), rectangle.get_y() + rectangle.get_height() / 2., label_text, ha="left", va="center", fontsize=10)

        if unordered:

            plt.yticks(yticks, procranks)

        else:

            ax.set_yticks(procranks)
            ax.set_yticklabels(procranks)

        # Format the axis ticks and labels
        if cpu:
            ax.set_xlabel('CPU time per thread (s)', fontsize='large')
        else:
            ax.set_xlabel('Time (s)', fontsize='large')

        if numberofproc:
            ax.set_ylabel('Number of processes', fontsize='large')
        else:
            ax.set_ylabel('Process rank', fontsize='large')

        #ax.yaxis.grid(True)

        if nprocs == 1:

            ax.set_frame_on(False)
            fig = plt.gcf()
            fig.set_size_inches(10,2)
            ax.xaxis.tick_bottom()
            ax.yaxis.set_visible(False)

        # Shrink current axis's height by 20% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])

        # Set the plot title
        plt.title("Timeline of the different simulation phases")

        # Put a legend below current axis
        ax.legend(legendEntries, legendNames, loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=False, ncol=4, prop={'size':12})

        # Save the figure
        plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.40)
        self._log.info("Created PDF timeline " + plotfilepath)

        # Close the plotting environment
        plt.close()

# -----------------------------------------------------------------
