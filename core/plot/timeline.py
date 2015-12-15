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
import os
import numpy as np
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

# -----------------------------------------------------------------

# Define the colors for the different simulation phases in the plot
colors = {"setup": 'r',         # setup -> red
          "stellar": 'g',       # stellar emission -> green
          "comm": '#FF7626',    # communication -> orange
          "spectra": 'm',       # spectra calculation -> magenta
          "dust": 'c',          # dust emission -> cyan
          "write": 'y',         # writing -> yellow
          "wait": 'b',          # waiting -> blue
          "other": 'k',         # other -> black
          None: 'k'}            # None -> black

# -----------------------------------------------------------------

class TimeLinePlotter(object):

    """
    An instance of the TimeLinePlotter class is used to create timeline diagrams for the different simulation phases
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        ## Attributes

        # The path to the output directory
        self.output_path = None

        # The table containing the timeline information
        self.table = None

    # -----------------------------------------------------------------

    def run(self, input, output_path):

        """
        This function should be called to invoke the plotting routine, which will create a timeline plot for each
        file that is found with valid timeline data
        :return:
        """

        # If the input is a Table object
        if isinstance(input, Table): self.table = input

        # If the input is a string
        elif isinstance(input, basestring):

            fill_values = ('--', '0', 'Simulation phase')
            self.table = ascii.read(input, fill_values=fill_values)

        # Invalid input
        else: raise ValueError("Input must be either an Astropy Table object or a filename (e.g. memory.dat)")

        # Set the path to the output directory
        self.output_path = output_path

        # Create the plots
        self.plot()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Get a list of the different process ranks
        ranks = np.unique(self.table["Process rank"])

        # Initialize a data structure to contain the start times and endtimes for the different processes,
        # indexed on the phase
        data = []

        # Iterate over the different entries in the timeline table
        for i in range(len(self.table)):

            if self.table["Process rank"][i] == 0:

                phase = self.table["Simulation phase"][i]
                phase = "other" if isinstance(phase, np.ma.core.MaskedConstant) else phase
                if phase == "start": phase = "other"

                data.append([phase, [], []])
                data[len(data) - 1][1].append(self.table["Start time"][i])
                data[len(data) - 1][2].append(self.table["End time"][i])

            else:

                nphases = len(data)
                data[i % nphases][1].append(self.table["Start time"][i])
                data[i % nphases][2].append(self.table["End time"][i])

        # Create the plot
        plot_path = os.path.join(self.output_path, "timeline.pdf")
        create_timeline_plot(data, plot_path, ranks)

# -----------------------------------------------------------------

def create_timeline_plot(data, path, procranks, figsize=(12, 8), percentages=False, totals=False, unordered=False, numberofproc=False, cpu=False):

    """
    This function actually plots the timeline based on a data structure containing the starttimes and endtimes
    for the different simulation phases
    :param data:
    :param path:
    :param procranks:
    :param figsize:
    :param percentages:
    :param totals:
    :param unordered:
    :param numberofproc:
    :param cpu:
    :return:
    """

    # Initialize figure
    plt.figure(figsize=figsize)
    plt.clf()

    ax = plt.gca()

    legend_entries = []
    legend_names = []
    unique_phases = []   # A LIST OF THE UNIQUE PHASE NAMES

    # Determine the number of processes
    nprocs = len(procranks)

    # Get the ordering
    if unordered: yticks = np.array(procranks).argsort().argsort()
    else: yticks = procranks

    #print("yticks=", yticks)
    #print("durations=", durations)

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

        if phase not in unique_phases and not (phase == "comm" and nprocs == 1):

            unique_phases.append(phase)
            legend_entries.append(patch_handle)
            legend_names.append(phase)

    if percentages:

        # For the different phases
        for phase, patch_handle in enumerate(patch_handles):

            durations = durations_list[phase]

            for sorting_number, rectangle in enumerate(patch_handle.get_children()):

                duration = durations[sorting_number]
                percentage = float(duration) / float(totaldurations[sorting_number]) * 100.0

                x = 0.5 * rectangle.get_width() + rectangle.get_x()
                y = 0.5 * rectangle.get_height() + rectangle.get_y()

                if rectangle.get_width() > 2000:

                    plt.text(x, y, "%d%%" % percentage, ha='center', va='center', fontsize=10)

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
    if cpu: ax.set_xlabel('CPU time (s)', fontsize='large')
    else: ax.set_xlabel('Time (s)', fontsize='large')

    if numberofproc: ax.set_ylabel('Number of processes', fontsize='large')
    else: ax.set_ylabel('Process rank', fontsize='large')

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
    ax.legend(legend_entries, legend_names, loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=False, ncol=4, prop={'size':12})

    # Save the figure
    plt.savefig(path, bbox_inches="tight", pad_inches=0.40)
    plt.close()

# -----------------------------------------------------------------
