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
from collections import defaultdict

# Import the relevant PTS classes and modules
from .plotter import Plotter
from ..tools.logging import log
from ..tools import filesystem as fs
from ..extract.timeline import TimeLineExtractor
from ..basics.configurable import Configurable
from ..simulation.simulation import load_simulations
from ..basics.map import Map

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

        # Prepare
        #self.prepare_multi()

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

    def prepare_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the input data into plottable format...")

        # Initialize a data structure to contain the scaling information in plottable format
        data = defaultdict(lambda: defaultdict(lambda: Map({"processor_counts": [], "times": [], "errors": []})))

        # Create an attribute to store the serial runtimes
        serial = defaultdict(lambda: Map({"time": None, "error": None}))

        sigma_level = 1.0

        # Create a dictionary to store the runtimes for the different simulation phases of the simulations
        # performed on one processor
        serial_times = defaultdict(list)

        # Keep track of the different processor counts encountered for the different parallelization modes
        modes = defaultdict(set)

        # Create dictionaries to contain the data before it is averaged over the different simulations
        total_times = defaultdict(lambda: defaultdict(list))
        setup_times = defaultdict(lambda: defaultdict(list))
        stellar_times = defaultdict(lambda: defaultdict(list))
        spectra_times = defaultdict(lambda: defaultdict(list))
        dust_times = defaultdict(lambda: defaultdict(list))
        writing_times = defaultdict(lambda: defaultdict(list))
        waiting_times = defaultdict(lambda: defaultdict(list))
        communication_times = defaultdict(lambda: defaultdict(list))
        memory = defaultdict(lambda: defaultdict(list))

        # Loop over the different entries in the scaling table
        for i in range(len(self.table)):

            # Get the parallelization mode
            mode = self.table["Parallelization mode"][i]

            # Get the number of processes and threads
            processes = self.table["Processes"][i]
            threads = self.table["Threads"][i]
            processors = processes * threads

            # If the number of processors is 1, add the runtimes for the different simulation phases to the
            # dictionary that contains the serial runtimes
            if processors == 1:
                serial_times["total"].append(self.table["Total time"][i])
                serial_times["setup"].append(self.table["Setup time"][i])
                serial_times["stellar"].append(self.table["Stellar emission time"][i])
                serial_times["spectra"].append(self.table["Spectra calculation time"][i])
                serial_times["dust"].append(self.table["Dust emission time"][i])
                serial_times["writing"].append(self.table["Writing time"][i])
                serial_times["waiting"].append(self.table["Waiting time"][i])
                serial_times["communication"].append(self.table["Communication time"][i])

            # Add the processor count for this parallelization mode
            modes[mode].add(processors)

            # Fill in the runtimes and memory usage at the appropriate place in the dictionaries
            total_times[mode][processors].append(self.table["Total time"][i])
            setup_times[mode][processors].append(self.table["Setup time"][i])
            stellar_times[mode][processors].append(self.table["Stellar emission time"][i])
            spectra_times[mode][processors].append(self.table["Spectra calculation time"][i])
            dust_times[mode][processors].append(self.table["Dust emission time"][i])
            writing_times[mode][processors].append(self.table["Writing time"][i])
            waiting_times[mode][processors].append(self.table["Waiting time"][i])
            communication_times[mode][processors].append(self.table["Communication time"][i])
            memory[mode][processors].append(self.table["Peak memory usage"][i])

        # Average the serial runtimes, loop over each phase
        for phase in serial_times:

            self.serial[phase].time = np.mean(serial_times[phase])
            self.serial[phase].error = sigma_level * np.std(serial_times[phase])

        # Loop over all encountered parallelization modes
        for mode in modes:

            # Loop over all processor counts encountered for this mode
            for processors in modes[mode]:
                # Average the runtimes for the different simulation phases and the memory usage for the different
                # runs for a certain parallelization mode and number of processors
                self.data["total"][mode].processor_counts.append(processors)
                self.data["total"][mode].times.append(np.mean(total_times[mode][processors]))
                self.data["total"][mode].errors.append(sigma_level * np.std(total_times[mode][processors]))

                self.data["setup"][mode].processor_counts.append(processors)
                self.data["setup"][mode].times.append(np.mean(setup_times[mode][processors]))
                self.data["setup"][mode].errors.append(sigma_level * np.std(setup_times[mode][processors]))

                self.data["stellar"][mode].processor_counts.append(processors)
                self.data["stellar"][mode].times.append(np.mean(stellar_times[mode][processors]))
                self.data["stellar"][mode].errors.append(sigma_level * np.std(stellar_times[mode][processors]))

                self.data["spectra"][mode].processor_counts.append(processors)
                self.data["spectra"][mode].times.append(np.mean(spectra_times[mode][processors]))
                self.data["spectra"][mode].errors.append(sigma_level * np.std(spectra_times[mode][processors]))

                self.data["dust"][mode].processor_counts.append(processors)
                self.data["dust"][mode].times.append(np.mean(dust_times[mode][processors]))
                self.data["dust"][mode].errors.append(sigma_level * np.std(dust_times[mode][processors]))

                self.data["writing"][mode].processor_counts.append(processors)
                self.data["writing"][mode].times.append(np.mean(writing_times[mode][processors]))
                self.data["writing"][mode].errors.append(sigma_level * np.std(writing_times[mode][processors]))

                self.data["waiting"][mode].processor_counts.append(processors)
                self.data["waiting"][mode].times.append(np.mean(waiting_times[mode][processors]))
                self.data["waiting"][mode].errors.append(sigma_level * np.std(waiting_times[mode][processors]))

                self.data["communication"][mode].processor_counts.append(processors)
                self.data["communication"][mode].times.append(np.mean(communication_times[mode][processors]))
                self.data["communication"][mode].errors.append(
                    sigma_level * np.std(communication_times[mode][processors]))

                self.data["memory"][mode].processor_counts.append(processors)
                self.data["memory"][mode].times.append(np.mean(memory[mode][processors]))
                self.data["memory"][mode].errors.append(sigma_level * np.std(memory[mode][processors]))

        # Determine the name of the system used for the scaling test
        self.system_name = os.path.basename(self.output_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        self.plot_single()

        #self.plot_multi()

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

    def plot_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the scaling timeline...")

        # Loop over the different parallelization modes
        for mode in self.data["total"]:

            # Determine the path to the timeline plot file
            plot_file_path = fs.join(self.output_path, "timeline_" + mode + ".pdf")

            # Initialize a data structure to contain the start times and endtimes of the different simulation phases,
            # for the different processor counts (data is indexed on the simulation phase)
            data = []
            nprocs_list = []

            # Loop over the different processor counts
            for j in range(len(self.data["total"][mode].processor_counts)):

                # Get the processor count
                processors = self.data["total"][mode].processor_counts[j]

                # Get the average runtimes for the different phases corresponding to the current processor count
                setup_time = self.data["setup"][mode].times[j] * processors
                stellar_time = self.data["stellar"][mode].times[j] * processors
                spectra_time = self.data["spectra"][mode].times[j] * processors
                dust_time = self.data["dust"][mode].times[j] * processors
                writing_time = self.data["writing"][mode].times[j] * processors
                waiting_time = self.data["waiting"][mode].times[j] * processors
                communication_time = self.data["communication"][mode].times[j] * processors

                # Add the processor count
                nprocs_list.append(processors)

                total = 0.0

                # For the first processor count
                if j == 0:

                    data.append(["setup", [total], [total + setup_time]])
                    total += setup_time
                    data.append(["stellar", [total], [total + stellar_time]])
                    total += stellar_time
                    data.append(["spectra", [total], [total + spectra_time]])
                    total += spectra_time
                    data.append(["dust", [total], [total + dust_time]])
                    total += dust_time
                    data.append(["write", [total], [total + writing_time]])
                    total += writing_time
                    data.append(["wait", [total], [total + waiting_time]])
                    total += waiting_time
                    data.append(["comm", [total], [total + communication_time]])
                    total += communication_time

                else:

                    # Setup
                    data[0][1].append(total)
                    total += setup_time
                    data[0][2].append(total)

                    # Stellar
                    data[1][1].append(total)
                    total += stellar_time
                    data[1][2].append(total)

                    # Spectra
                    data[2][1].append(total)
                    total += spectra_time
                    data[2][2].append(total)

                    # Dust
                    data[3][1].append(total)
                    total += dust_time
                    data[3][2].append(total)

                    # Writing
                    data[4][1].append(total)
                    total += writing_time
                    data[4][2].append(total)

                    # Waiting
                    data[5][1].append(total)
                    total += waiting_time
                    data[5][2].append(total)

                    # Communication
                    data[6][1].append(total)
                    total += communication_time
                    data[6][2].append(total)

            # Set the plot title
            title = "Scaling timeline for " + self.system_name

            # Create the plot
            create_timeline_plot(data, plot_file_path, nprocs_list, percentages=True, totals=True, unordered=True,
                                 numberofproc=True, cpu=True, title=title)


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
