#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.timeline Contains the TimeLinePlotter class, which is used to create timeline diagrams
#  of the different phases of a SKIRT simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .plotter import Plotter
from ..tools.logging import log
from ..tools import filesystem as fs
from ..extract.timeline import TimeLineExtractor
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

# Define the colors for the different simulation phases in the plot
colors = {"setup": 'r',         # setup -> red
          "stellar": 'g',       # stellar emission -> green
          "comm": '#FF7626',    # communication -> orange
          "spectra": 'm',       # spectra calculation -> magenta
          "dust": 'c',          # dust emission -> cyan
          "write": 'y',         # writing -> yellow
          "wait": 'b',          # waiting -> blue
          "other": 'k'}         # other -> black

# Define the names identifying the different phases in the plot
phase_label_names = {"setup": "setup",
                     "stellar": "stellar",
                     "comm": "communication",
                     "spectra": "spectra",
                     "dust": "dust",
                     "write": "write",
                     "wait": "waiting",
                     "other": "other"}

# -----------------------------------------------------------------

class BatchTimeLinePlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(BatchTimeLinePlotter, self).__init__(config)

        # The simulations
        self.simulations = []

        # The timelines
        self.timelines = []

        # The data
        self.single_data = []
        self.multi_data = []

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Extract
        self.extract()

        # Prepare
        self.prepare()

        # 3. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchTimeLinePlotter, self).setup(**kwargs)

        # Load simulations from working directory if none have been added
        if len(self.simulations) == 0: self.simulations = load_simulations(self.config.path)

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the timelines ...")

        # Loop over the simulations
        for simulation in self.simulations:

            # Create a TimeLineExtractor instance
            extractor = TimeLineExtractor()

            # Run the timeline extractor
            timeline = extractor.run(simulation)

            # Add the timeline
            self.timelines.append(timeline)

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the data for plotting ...")

        # Prepare
        self.prepare_single()

    # -----------------------------------------------------------------

    def prepare_single(self):

        """
        This function ...
        :return:
        """

        # Loop over the timelines
        for timeline in self.timelines:

            # Get a list of the different process ranks
            ranks = np.unique(timeline["Process rank"])

            # Initialize the data structure to contain the start times and endtimes for the different processes,
            # indexed on the phase
            data = []

            # Iterate over the different entries in the timeline table
            for i in range(len(timeline)):

                if timeline["Process rank"][i] == 0:

                    phase = timeline["Simulation phase"][i]

                    # Few special cases where we want the phase indicator to just say 'other'
                    if phase is None or phase == "start" or isinstance(phase, np.ma.core.MaskedConstant): phase = "other"

                    # Add the data
                    data.append([phase, [], []])
                    data[len(data) - 1][1].append(timeline["Start time"][i])
                    data[len(data) - 1][2].append(timeline["End time"][i])

                else:

                    nphases = len(data)
                    data[i % nphases][1].append(timeline["Start time"][i])
                    data[i % nphases][2].append(timeline["End time"][i])

            # Add the data
            self.single_data.append((ranks, data))

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        self.plot_single()

    # -----------------------------------------------------------------

    def plot_single(self):

        """
        This function ...
        :return:
        """

        # Loop over the data
        for ranks, data in self.single_data:

            # Determine the path
            #plot_path = fs.join(self.output_path, "timeline.pdf")

            plot_path = None

            # Create the plot
            create_timeline_plot(data, plot_path, ranks)

# -----------------------------------------------------------------

class TimeLinePlotter(Plotter):

    """
    An instance of the TimeLinePlotter class is used to create timeline diagrams for the different simulation phases
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(TimeLinePlotter, self).__init__()

        # -- Attributes --

        # A list of the process ranks
        self.ranks = None

    # -----------------------------------------------------------------

    @staticmethod
    def default_input():

        """
        This function ...
        :return:
        """

        return "timeline.dat"

    # -----------------------------------------------------------------

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        # Get a list of the different process ranks
        self.ranks = np.unique(self.table["Process rank"])

        # Initialize the data structure to contain the start times and endtimes for the different processes,
        # indexed on the phase
        self.data = []

        # Iterate over the different entries in the timeline table
        for i in range(len(self.table)):

            if self.table["Process rank"][i] == 0:

                phase = self.table["Simulation phase"][i]

                # Few special cases where we want the phase indicator to just say 'other'
                if phase is None or phase == "start" or isinstance(phase, np.ma.core.MaskedConstant): phase = "other"

                # Add the data
                self.data.append([phase, [], []])
                self.data[len(self.data) - 1][1].append(self.table["Start time"][i])
                self.data[len(self.data) - 1][2].append(self.table["End time"][i])

            else:

                nphases = len(self.data)
                self.data[i % nphases][1].append(self.table["Start time"][i])
                self.data[i % nphases][2].append(self.table["End time"][i])

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Making the plots...")

        # Create the plot
        plot_path = fs.join(self.output_path, "timeline.pdf")
        create_timeline_plot(self.data, plot_path, self.ranks)

# -----------------------------------------------------------------

def create_timeline_plot(data, path, procranks, figsize=(12, 8), percentages=False, totals=False, unordered=False, numberofproc=False, cpu=False, title=None):

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
            legend_names.append(phase_label_names[phase])

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
    if title is None: plt.title("Timeline of the different simulation phases")
    else: plt.title(title)

    # Put a legend below current axis
    ax.legend(legend_entries, legend_names, loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=False, ncol=4, prop={'size': 12})

    # Save the figure
    if path is not None: plt.savefig(path, bbox_inches="tight", pad_inches=0.40)
    else: plt.show()
    plt.close()

# -----------------------------------------------------------------
