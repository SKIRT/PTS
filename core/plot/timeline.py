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
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..tools import filesystem as fs
from ..extract.timeline import extract_timeline
from ..basics.configurable import Configurable
from ..simulation.discover import SimulationDiscoverer
from ..simulation.parallelization import Parallelization
from ..tools.stringify import tostr

# -----------------------------------------------------------------

rc('text', usetex=True)

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

def plot_timeline(timeline, path=None, write_timelines=False, write_data=False):

    """
    This function ...
    :param timeline:
    :param path:
    :param write_timelines:
    :param write_data:
    :return:
    """

    # Create a TimeLinePlotter object
    plotter = TimeLinePlotter()

    # Set the output path
    plotter.config.output = path

    # Set options
    plotter.config.write_timelines = write_timelines
    plotter.config.write_data = write_data

    # Run the timeline plotter
    plotter.run(timeline=timeline)

# -----------------------------------------------------------------

class TimeLinePlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(TimeLinePlotter, self).__init__(*args, **kwargs)

        # The simulations
        self.simulations = []

        # The timelines
        self.timelines = dict()

        # The parallelization schemes
        self.parallelizations = dict()

        # The data
        self.single_data = dict()
        self.multi_data = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Extract
        if not self.is_extracted: self.extract()

        # 3. Prepare
        if not self.is_prepared: self.prepare()

        # 4. Writing
        self.write()

        # 5. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TimeLinePlotter, self).setup(**kwargs)

        # Get the timeline(s)
        if "timeline" in kwargs:
            output_path = self.config.output if self.config.output is not None else self.config.path
            self.timelines[output_path] = kwargs.pop("timeline")
        elif "timelines" in kwargs:
            timelines = kwargs.pop("timelines")
            if isinstance(timelines, list): raise ValueError("Timelines must be a dictionary")
            elif isinstance(timelines, dict): self.timelines = timelines
            else: raise ValueError("Timelines must be a dictionary")

        # Load simulations from working directory if none have been added
        if len(self.simulations) == 0 and len(self.timelines) == 0: self.load_simulations()

    # -----------------------------------------------------------------

    def load_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading simulations ...")

        # Create the simulation discoverer
        discoverer = SimulationDiscoverer()
        discoverer.config.path = self.config.path
        discoverer.config.list = False

        # Run the simulation discoverer
        discoverer.run()

        # Set the simulations
        if self.config.hetero: self.simulations = discoverer.simulations
        else: self.simulations = discoverer.simulations_single_ski

    # -----------------------------------------------------------------

    @property
    def simulation_prefix(self):

        """
        This function ...
        :return:
        """

        return self.simulations[0].prefix()

    # -----------------------------------------------------------------

    @property
    def is_extracted(self):

        """
        This function ...
        :return:
        """

        return len(self.timelines) > 0

    # -----------------------------------------------------------------

    @property
    def is_prepared(self):

        """
        This function ...
        :return:
        """

        return len(self.single_data) > 0

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

            # Get the timeline
            timeline = extract_timeline(simulation)

            # Get the simulation output path
            output_path = simulation.output_path

            # Check whether unique
            if output_path in self.timelines: raise RuntimeError("Multiple simulations have their output in the '" + output_path + "' directory")

            # Get the log file
            log_file = simulation.log_file

            # Get properties
            nprocesses = log_file.nprocesses
            nthreads = log_file.nthreads
            data_parallel = log_file.data_parallel

            # Create parallelization scheme
            parallelization = Parallelization.from_processes_and_threads(nprocesses, nthreads, data_parallel=data_parallel)

            # Set the parallelization
            self.parallelizations[output_path] = parallelization

            # Add the timeline
            self.timelines[output_path] = timeline

    # -----------------------------------------------------------------

    @property
    def has_multi(self):

        """
        This function ...
        :return:
        """

        return len(self.timelines) > 1

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the data for plotting ...")

        # Prepare data for single-simulation plots
        self.prepare_single()

        # Prepare for multi-simulation plot
        if self.has_multi: self.prepare_multi()

    # -----------------------------------------------------------------

    def prepare_single(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the timeline data for making single-simulation plots ...")

        # Loop over the timelines
        for output_path in self.timelines:

            # Get the timeline
            timeline = self.timelines[output_path]

            # Check the timeline
            check_timeline(timeline)

            # Get a list of the different process ranks
            ranks = np.unique(timeline["Process rank"])

            # Initialize the data structure to contain the start times and endtimes for the different processes,
            # indexed on the phase
            data = []

            # Skipped entries
            skipped = []

            nphases = 0

            counter = 0

            previous_rank = 0

            # Iterate over the different entries in the timeline table
            for i in range(len(timeline)):

                rank = timeline["Process rank"][i]

                if rank == 0:

                    nphases += 1

                    # Few special cases where we want the phase indicator to just say 'other'
                    phase = timeline["Phase"][i]
                    if phase is None or phase == "start" or isinstance(phase, np.ma.core.MaskedConstant): phase = "other"

                    # Don't plot 'other' phases
                    if phase == "other" and not self.config.other:
                        skipped.append(i)
                        continue

                    # Add the data
                    data.append([phase, [], []])
                    data[len(data) - 1][1].append(timeline["Start time"][i])
                    data[len(data) - 1][2].append(timeline["End time"][i])

                else:

                    if rank != previous_rank: counter = 0

                    # Few special cases where we want the phase indicator to just say 'other'
                    phase = timeline["Phase"][i]
                    if phase is None or phase == "start" or isinstance(phase, np.ma.core.MaskedConstant): phase = "other"

                    # Don't plot 'other' phases
                    if phase == "other" and not self.config.other:
                        skipped.append(i)
                        continue

                    #index = i % nphases
                    #if index in skipped:
                    #    continue # skip skipped entries

                    index = counter

                    data[index][1].append(timeline["Start time"][i])
                    data[index][2].append(timeline["End time"][i])

                    counter += 1
                    previous_rank = rank

            # Add the data
            self.single_data[output_path] = (ranks, data)

    # -----------------------------------------------------------------

    def prepare_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the timeline data for making multi-simulation plots ...")

        # Initialize data structures
        data = []
        nprocs_list = []

        # Initialize phases
        data.append(["setup", [], []])
        data.append(["stellar", [], []])
        data.append(["spectra", [], []])
        data.append(["dust", [], []])
        data.append(["write", [], []])
        data.append(["wait", [], []])
        data.append(["comm", [], []])

        # The simulation names
        simulation_names = []

        # Loop over the timelines
        for output_path in self.timelines:

            # Get the timeline
            timeline = self.timelines[output_path]

            # Get the number of processes
            nprocesses = timeline.nprocesses

            # Get the average runtimes for the different phases corresponding to the current processor count
            setup_time = timeline.setup * nprocesses
            stellar_time = timeline.stellar * nprocesses
            spectra_time = timeline.spectra * nprocesses
            dust_time = timeline.dust * nprocesses
            writing_time = timeline.writing * nprocesses
            waiting_time = timeline.waiting * nprocesses
            communication_time = timeline.communication * nprocesses

            total = 0.0

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

            # Add the process count
            nprocs_list.append(nprocesses)

            # Add the simulation name
            simulation_name = fs.name(output_path)
            simulation_names.append(simulation_name)

        # Set the data
        self.multi_data = (nprocs_list, simulation_names, data)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the extracted timelines
        if self.config.write_timelines: self.write_timelines()

        # Write the data
        if self.config.write_data: self.write_data()

    # -----------------------------------------------------------------

    def write_timelines(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the timelines ...")

        # Loop over the data
        for output_path in self.timelines:

            # Determine path
            if self.config.output is not None: path = fs.join(self.config.output, "timeline_" + fs.name(output_path) + ".dat")
            else: path = fs.join(output_path, "timeline.dat")

            # Write the timeline
            self.timelines[output_path].saveto(path)

    # -----------------------------------------------------------------

    def write_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the timeline data ...")

        # Write single data
        self.write_single_data()

        # Write multi data
        if self.has_multi: self.write_multi_data()

    # -----------------------------------------------------------------

    def write_single_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the single timeline data ...")

        # Determine the path
        if self.config.output is not None: path = fs.join(self.config.output, "single_data.dat")
        #else: path = fs.join(self.config.path, "single_data.dat")
        else: path = self.output_path_file("single_data.dat")

        # Write
        #write_dict(self.single_data, path)
        with open(path, "w") as datafile:

            # Loop over the data for the different timelines
            for output_path in self.single_data:

                ranks, data = self.single_data[output_path]
                print(output_path, file=datafile)
                print(tostr(ranks), file=datafile)
                #print(tostr(data), file=datafile)
                for phase, start_times, end_times in data:
                    print(phase + " " + str(start_times) + " " + str(end_times), file=datafile)
                    #print(tostr(data_phase), file=datafile)
                print("", file=datafile)

    # -----------------------------------------------------------------

    def write_multi_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the multi timeline data ...")

        # Determine the path
        if self.config.output is not None: path = fs.join(self.config.output, "multi_data.dat")
        #else: path = fs.join(self.config.path, "multi_data.dat")
        else: path = self.output_path_file("multi_data.dat")

        # Write
        #write_data_tuple(self.multi_data, path)
        with open(path, 'w') as datafile:

            # Loop over the data
            for nprocs_list, simulation_names, data in self.multi_data:

                print(tostr(nprocs_list), file=datafile)
                print(tostr(simulation_names), file=datafile)
                for phase, start_times, end_times in data:
                    print(phase + " " + str(start_times) + " " + str(end_times), file=datafile)
                print("", file=datafile)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot single
        self.plot_single()

        # Plot combined timelines
        if self.has_multi: self.plot_combined()

        # Plot multi
        if self.has_multi: self.plot_multi()

    # -----------------------------------------------------------------

    def plot_single(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting timelines of individual simulations ...")

        # Loop over the data
        for output_path in self.single_data:

            # Debugging
            log.debug("Plotting timeline for the " + fs.name(output_path) + " simulation ...")

            # Get the data
            ranks, data = self.single_data[output_path]

            # Determine path
            if self.config.output is not None: path = fs.join(self.config.output, "timeline_" + fs.name(output_path) + ".pdf")
            #else: path = fs.join(output_path, "timeline.pdf")
            else: path = None

            # Create the plot
            create_timeline_plot(data, ranks, path, figsize=self.config.figsize, percentages=self.config.percentages,
                                 totals=self.config.totals, unordered=False, cpu=False, title=self.config.title,
                                 ylabels=None, yaxis=None, rpc="r", add_border=self.config.add_border,
                                 show_ranks=self.config.show_ranks, label_fontsize=self.config.label_fontsize,
                                 ticks_fontsize=self.config.ticks_fontsize)

    # -----------------------------------------------------------------

    def plot_combined(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Combining timelines of different simulations with the same number of processes in the same plot ...")

        # Dictionary
        output_paths_per_nproc = defaultdict(list)

        # Loop over the data
        for output_path in self.single_data:

            # Get the number of processes
            ranks, data = self.single_data[output_path]
            nproc = len(ranks)

            # Add the output path to the dictionary
            output_paths_per_nproc[nproc].append(output_path)

        # Loop over the different nprocs
        for nproc in output_paths_per_nproc:

            # Determine path
            if self.config.output is not None: path = fs.join(self.config.output, "timelines_" + str(nproc) + "processes.pdf")
            #else: path = fs.join(self.config.path, "timelines_" + str(nproc) + "processes.pdf")
            else: path = None

            # Gather the data
            nproc_dict = dict()
            data_dict = dict()

            # Loop over the timelines
            for output_path in output_paths_per_nproc[nproc]:

                #print(output_path)
                #print(self.single_data)

                # Gather the data
                nproc_dict[output_path] = self.single_data[output_path][0]
                data_dict[output_path] = self.single_data[output_path][1]

            # Create the plot
            create_multiple_timelines_plot(data_dict, nproc_dict, self.parallelizations, path, show_ranks=False,
                                           title="Timelines for " + str(nproc) + " processes",
                                           figsize=self.config.figsize, ticks_fontsize=self.config.ticks_fontsize)

    # -----------------------------------------------------------------

    def plot_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a timeline of the CPU time of all simulations ...")

        # Set the plot title
        title = "Timeline of CPU time"

        # Get the data
        nprocs_list, simulation_names, data = self.multi_data

        # Determine the path
        if self.config.output is not None: path = fs.join(self.config.output, "timeline_cputime.pdf")
        #else: path = fs.join(self.config.path, "timeline_cputime.pdf")
        else: path = None

        # Create the plot
        create_timeline_plot(data, nprocs_list, path, percentages=True, totals=True, unordered=True, cpu=True,
                             title=title, ylabels=simulation_names, yaxis="Simulations")

# -----------------------------------------------------------------

def create_timeline_plot(data, procranks, path=None, figsize=(12, 8), percentages=False, totals=False, unordered=False,
                         cpu=False, title=None, ylabels=None, yaxis=None, rpc="r", add_border=False, show_ranks=True,
                         label_fontsize=18, title_fontsize=20, ticks_fontsize=12):

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
    :param cpu:
    :param title:
    :param ylabels:
    :param yaxis:
    :param rpc: 'rank', 'processes' or 'cores'
    :param add_border:
    :param show_ranks:
    :param label_fontsize:
    :param title_fontsize:
    :param ticks_fontsize:
    :return:
    """

    # Initialize figure
    plt.figure(figsize=figsize)
    plt.clf()

    ax = plt.gca()

    # Set x axis grid
    ax.xaxis.grid(linestyle="dotted", linewidth=2.0)

    legend_entries = []
    legend_names = []
    unique_phases = []   # A LIST OF THE UNIQUE PHASE NAMES

    # Determine the number of processes
    nprocs = len(procranks)

    # Get the ordering
    if unordered: yticks = np.array(procranks).argsort().argsort()
    else: yticks = procranks

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

        #print("YTICKS", yticks)
        #print("PROCRANKS", procranks)
        #plt.yticks(yticks, procranks)
        if rpc == 'r':
            if show_ranks:
                ax.set_yticks(yticks)
                ax.set_yticklabels(procranks)
            else:
                ax.set_yticks([])
                ax.set_yticklabels([])
        else:
            ax.set_yticks(yticks)
            ax.set_yticklabels(procranks)
    else:

        if rpc == 'r':
            if show_ranks:
                ax.set_yticks(procranks)
                ax.set_yticklabels(procranks)
            else:
                ax.set_yticks([])
                ax.set_yticklabels([])
        else:
            ax.set_yticks(procranks)
            ax.set_yticklabels(procranks)

    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    if not add_border:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis=u'both', which=u'both', length=0)

    # Format the axis ticks and labels
    if cpu: ax.set_xlabel('CPU time (s)', fontsize=label_fontsize)
    else: ax.set_xlabel('Time (s)', fontsize=label_fontsize)

    # Set y label
    if rpc == 'r':
        if show_ranks: ax.set_ylabel('Process rank', fontsize=label_fontsize)
        else: ax.set_ylabel("Processes", fontsize=label_fontsize)
    elif rpc == 'p': ax.set_ylabel('Number of processes', fontsize=label_fontsize)
    elif rpc == 'c': ax.set_ylabel('Number of cores', fontsize=label_fontsize)

    #ax.yaxis.grid(True)

    # Custom y labels
    if ylabels is not None:
        plt.yticks(yticks, ylabels)
        ax.set_ylabel("")

    # Custom y axis label
    if yaxis is not None: ax.set_ylabel(yaxis)

    # Set ticks fontsize
    plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=ticks_fontsize)
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=ticks_fontsize)

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
    #if title is None: plt.title("Timeline of the different simulation phases")
    #else: plt.title(title)
    if title is not None: plt.suptitle(title, fontsize=title_fontsize)

    # Put a legend below current axis
    legend = ax.legend(legend_entries, legend_names, loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=4, prop={'size': 12})

    # Change legend properties
    frame = legend.get_frame()
    frame.set_linewidth(0)
    frame.set_facecolor('0.85')
    legend.legendPatch.set_alpha(0.75)

    # Save the figure
    if path is not None: plt.savefig(path, bbox_inches="tight", pad_inches=0.40)
    else: plt.show()
    plt.close()

# -----------------------------------------------------------------

def create_multiple_timelines_plot(data_dict, procranks_dict, parallelization, path=None, figsize=(12, 8), unordered=False,
                                   cpu=False, title=None, ylabels=None, yaxis=None, add_border=False, show_ranks=True,
                                   label_fontsize=18, title_fontsize=20, ticks_fontsize=12):

    """
    This function ...
    :param data_dict:
    :param procranks_dict:
    :param path:
    :param figsize:
    :param unordered:
    :param cpu:
    :param title:
    :param ylabels:
    :param yaxis:
    :param rpc:
    :param add_border:
    :param show_ranks:
    :param label_fontsize:
    :param title_fontsize:
    :param ticks_fontsize:
    :return:
    """

    # Initialize figure
    fig = plt.figure(figsize=figsize)
    plt.clf()

    ax = plt.gca()
    fig.subplots_adjust(hspace=0.2)

    #   subplot(211)
    #produces a subaxes in a figure which represents the top plot (i.e. the
    #first) in a 2 row by 1 column notional grid

    ntimelines = len(data_dict)

    legend_entries = []
    legend_names = []
    unique_phases = []  # A LIST OF THE UNIQUE PHASE NAMES

    shared_axis = None
    last_axis = None
    for index in range(ntimelines):

        output_path = data_dict.keys()[index]
        data = data_dict[output_path]

        # Add subplot
        subplot_spec = ntimelines * 100 + 10 + index + 1
        ax = plt.subplot(subplot_spec, sharex=shared_axis)

        # Set x axis grid
        ax.xaxis.grid(linestyle="dotted", linewidth=2.0)

        # Determine the number of processes
        procranks = procranks_dict[output_path]
        nprocs = len(procranks)

        durations_list = []
        totaldurations = np.zeros(nprocs)
        patch_handles = []

        # Get the ordering
        if unordered: yticks = np.array(procranks).argsort().argsort()
        else: yticks = procranks

        # Make the timeline plot, consisting of a set of bars of the same color for each simulation phase
        for phase, starttimes, endtimes in data:

            durations = np.array(endtimes) - np.array(starttimes)
            durations_list.append(durations)

            totaldurations += durations

            patch_handle = ax.barh(yticks, durations, color=colors[phase], align='center', left=starttimes, alpha=0.8, lw=0)
            patch_handles.append(patch_handle)

            if index == 0:
                if phase not in unique_phases and not (phase == "comm" and nprocs == 1):
                    unique_phases.append(phase)
                    legend_entries.append(patch_handle)
                    legend_names.append(phase_label_names[phase])

        #plt.plot(t, s1)

        # Set axis limits
        #ax.set_xlim([xmin, xmax])
        ax.set_ylim([-0.5, nprocs-0.5])

        # Hide process ranks
        if not show_ranks:
            ax.set_yticks([])
            ax.set_yticklabels([])

        # Set x label
        if ax.is_last_row(): ax.set_xlabel('Time (s)', fontsize=label_fontsize)

        # Set y label
        ax.set_ylabel("Processes")

        # Set title
        if parallelization[output_path].nprocesses > 1:
            if parallelization[output_path].nthreads > 1:
                subplot_title = "Hybrid task+data parallelization" if parallelization[output_path].data_parallel else "Hybrid task parallelization"
            else:
                subplot_title = "Multiprocessing task+data parallelization" if parallelization[output_path].data_parallel else "Multiprocessing task parallelization"
        else: subplot_title = None
        ax.set_title(subplot_title)

        # Set axes tick formatter
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))

        # Remove border if requested
        if not add_border:

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis=u'both', which=u'both', length=0)

        if ax.is_last_row():
            # Set ticks fontsize
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=ticks_fontsize)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=ticks_fontsize)
            last_axis = ax
        else:
            # make these tick labels invisible
            plt.setp(ax.get_xticklabels(), visible=False)

        if index == 0: shared_axis = ax

    # Set the plot title
    if title is not None: plt.suptitle(title, fontsize=title_fontsize)

    # Put a legend below current axis
    legend = last_axis.legend(legend_entries, legend_names, loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=4, prop={'size': 12})

    # Change legend properties
    frame = legend.get_frame()
    frame.set_linewidth(0)
    frame.set_facecolor('0.85')
    legend.legendPatch.set_alpha(0.75)

    # Save the figure
    if path is not None: plt.savefig(path, bbox_inches="tight", pad_inches=0.40)
    else: plt.show()
    plt.close()

# -----------------------------------------------------------------

def check_timeline(timeline):

    """
    This function ...
    :param timeline:
    :return:
    """

    phases = []

    # Iterate over the different entries in the timeline table
    for i in range(len(timeline)):

        if timeline["Process rank"][i] == 0:

            phase = timeline["Phase"][i]

            phases.append(phase)

        else:

            nphases = len(phases)
            index = i % nphases

            phase = timeline["Phase"][i]

            #print(phases[index], phase)

            if phases[index] != phase: raise RuntimeError("Timeline is not consistent")

# -----------------------------------------------------------------
