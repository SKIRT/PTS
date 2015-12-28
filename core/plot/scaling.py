#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.scaling The class ScalingPlotter in this module makes plots of the results SKIRT
#  scaling benchmark tests performed with the scalingtest module.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import os.path
import matplotlib
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from collections import defaultdict

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from .plotter import Plotter
from ..basics.quantity import Quantity
from ..basics.map import Map
from .timeline import create_timeline_plot

# -----------------------------------------------------------------

phase_names = {"total": "total simulation", "setup": "simulation setup", "stellar": "stellar emission phase",
               "spectra": "calculation of dust emission spectra", "dust": "dust emission phase",
               "writing": "writing phase", "waiting": "waiting phases", "communication": "communication phases"}

phase_labels = {"total": "Total runtime", "setup": "Setup time", "stellar": "Stellar runtime",
                "spectra": "Runtime of dust spectra calculation", "dust": "Dust emission runtime",
                "writing": "Writing time", "waiting": "Waiting time", "communication": "Communication time"}

# -----------------------------------------------------------------

class ScalingPlotter(Plotter):

    """
    An instance of the ScalingPlotter is used to create plots of the runtimes, speedups and efficiencies as a function
    of the number of threads, based on the results of (multiple) SKIRT scaling benchmark tests.
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(ScalingPlotter, self).__init__()

        # -- Attributes --

        # A data structure to store the serial runtimes
        self.serial = None

        # The name of the system used for the scaling test
        self.system_name = None

    # -----------------------------------------------------------------

    @staticmethod
    def fill_values():

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    @staticmethod
    def default_input():

        """
        This function ...
        :return:
        """

        return "scaling.dat"

    # -----------------------------------------------------------------

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Preparing the input data into plottable format...")

        # Initialize a data structure to contain the scaling information in plottable format
        self.data = defaultdict(lambda: defaultdict(lambda: Map({"processor_counts": [], "times": [], "errors": []})))

        # Create an attribute to store the serial runtimes
        self.serial = defaultdict(lambda: Map({"time": None, "error": None}))

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
                self.data["communication"][mode].errors.append(sigma_level * np.std(communication_times[mode][processors]))

                self.data["memory"][mode].processor_counts.append(processors)
                self.data["memory"][mode].times.append(np.mean(memory[mode][processors]))
                self.data["memory"][mode].errors.append(sigma_level * np.std(memory[mode][processors]))

        # Determine the name of the system used for the scaling test
        self.system_name = os.path.basename(self.output_path)

    # -----------------------------------------------------------------

    @property
    def has_serial(self):

        """
        This function ...
        :return:
        """

        return len(self.serial) != 0

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Making the plots...")

        # Plot the runtimes, speedups and efficiencies for the total simulation
        self.plot_total()

        # Plot the runtimes, speedups and efficiencies for the different simulation phases
        self.plot_phases()

        # Plot the scaling relation of the memory usage
        self.plot_memory()

        # Plot a timeline of the CPU time spent in the different simulation phases
        self.plot_timeline()

    # -----------------------------------------------------------------

    def plot_total(self):

        """
        This function ...
        :return:
        """

        # Plot the total runtimes
        times_path = os.path.join(self.output_path, "times.pdf")
        self.plot_times("total", times_path)

        # Plot the speedups and efficiencies if serial runtimes are available
        if self.has_serial:

            # Plot the total speedups
            speedups_path = os.path.join(self.output_path, "speedups.pdf")
            self.plot_speedups("total", speedups_path)

            # Plot the total efficiencies
            efficiencies_path = os.path.join(self.output_path, "efficiencies.pdf")
            self.plot_efficiencies("total", efficiencies_path)

    # -----------------------------------------------------------------

    def plot_phases(self):

        """
        This function ...
        :return:
        """

        # Loop over the different simulation phases and plot the runtimes, speedups and efficiencies
        for phase in self.data:

            # Skip the total simulation runtimes and the memory information
            if phase == "total" or phase == "memory": continue

            # Make a seperate directory in the output directory to contain the plots for the current phase
            output_path_phase = os.path.join(self.output_path, phase)
            if not os.path.isdir(output_path_phase): os.mkdir(output_path_phase)

            # Plot the total runtimes
            times_path = os.path.join(output_path_phase, "times.pdf")
            self.plot_times(phase, times_path)

            # Plot the speedups and efficiencies if serial runtimes are available
            if self.has_serial:

                # Plot the total speedups
                speedups_path = os.path.join(output_path_phase, "speedups.pdf")
                self.plot_speedups("total", speedups_path)

                # Plot the total efficiencies
                efficiencies_path = os.path.join(output_path_phase, "efficiencies.pdf")
                self.plot_efficiencies("total", efficiencies_path)

    # -----------------------------------------------------------------

    def plot_times(self, phase, file_path, size=(12, 8)):

        """
        This function ...
        :param phase:
        :param file_path:
        :param size:
        :return:
        """

        # Inform the user
        self.log.info("Plotting the runtimes for the " + phase_names[phase] + "...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=size)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parallelization modes (the different curves)
        for mode in self.data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.data[phase][mode].processor_counts
            times = self.data[phase][mode].times
            errors = self.data[phase][mode].errors

            # Plot the data points for this mode
            plt.errorbar(processor_counts, times, errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        plt.grid(True)

        # Add axis labels and a legend
        plt.xlabel("Number of processors $n$", fontsize='large')
        plt.ylabel(phase_labels[phase] + " T (s)", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Scaling of the " + phase_labels[phase].lower() + " for " + self.system_name)

        # Save the figure
        plt.savefig(file_path)
        plt.close()

    # -----------------------------------------------------------------

    def plot_speedups(self, phase, file_path, size=(12, 8), plot_fit=True):

        """
        This function ...
        :param phase:
        :param file_path:
        :param size:
        :param plot_fit:
        :return:
        """

        # Inform the user of the fact that the speedups are being calculated and plotted
        self.log.info("Calculating and plotting the speedups for the " + phase_names[phase] + "...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=size)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Create a dictionary that stores the fitted parameters for each different mode
        parameters = dict()

        # Get the serial runtime (and error) for this phase (create a Quantity object)
        serial_time = self.serial[phase].time
        serial_error = self.serial[phase].error
        serial = Quantity(serial_time, serial_error)

        # Loop over the different parallelization modes (the different curves)
        for mode in self.data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.data[phase][mode].processor_counts
            times = self.data[phase][mode].times
            errors = self.data[phase][mode].errors

            # Calculate the speedups and the errors on the speedups
            speedups = []
            speedup_errors = []
            for i in range(len(processor_counts)):

                # Create a quantity for the current runtime
                time = Quantity(times[i], errors[i])

                # Calculate the speedup based on the current runtime and the serial runtime
                speedup = serial / time

                # Add the value and the propagated error of the speedup to the appropriate lists
                speedups.append(speedup.value)
                speedup_errors.append(speedup.error)

            # Plot the data points for this curve
            plt.errorbar(processor_counts, speedups, speedup_errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

            # Set the weights of the different speedup points for the fitting procedure
            speedup_weigths = speedup_errors if not np.any(np.isinf(speedup_errors)) else None

            # Fit (standard or modified) Amdahl's law to the speedups
            if len(processor_counts) < 10:

                popt, pcov = curve_fit(Amdahl, processor_counts, speedups, sigma=speedup_weigths, absolute_sigma=False)
                perr = np.sqrt(np.diag(pcov))
                parameters[mode] = Map({"p": popt[0], "p_error": perr[0], "a": 0.0, "a_error": 0.0, "b": 0.0, "b_error": 0.0, "c": 0.0, "c_error": 0.0})

            else:

                popt, pcov = curve_fit(modAmdahl, processor_counts, speedups, sigma=speedup_weigths, absolute_sigma=False)
                perr = np.sqrt(np.diag(pcov))
                parameters[mode] = Map({"p": popt[0], "p_error": perr[0], "a": popt[1], "a_error": perr[1], "b": popt[2], "b_error": perr[2], "c": popt[3], "c_error": perr[3]})

        # Use a logarithmic scale for both axes
        plt.xscale('log')
        plt.yscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Create a data file to contain the fitted parameters
        directory = os.path.dirname(file_path)
        parameter_file_path = os.path.join(directory, "parameters.dat")

        # Fit parameters for the speedups to Amdahl's law
        #  S_n = 1 / ( 1 - p + p/n + a + b*n + c*n^2 ) \n")
        mode_list = []
        p_list = []
        p_error_list = []
        a_list = []
        a_error_list = []
        b_list = []
        b_error_list = []
        c_list = []
        c_error_list = []

        for mode in parameters:

            mode_list.append(mode)
            p_list.append(parameters[mode].p)
            p_error_list.append(parameters[mode].p_error)
            a_list.append(parameters[mode].a)
            a_error_list.append(parameters[mode].a_error)
            b_list.append(parameters[mode].b)
            b_error_list.append(parameters[mode].b_error)
            c_list.append(parameters[mode].c)
            c_error_list.append(parameters[mode].c_error)

        # Create the parameters table and write to file
        data = [p_list, p_error_list, a_list, a_error_list, b_list, b_error_list, c_list, c_error_list]
        names = ["Parallel fraction p", "Error on p", "Parameter a", "Error on a", "Parameter b", "Error on b", "Parameter c", "Error on c"]
        table = Table(data, names)
        table.write(parameter_file_path, format="ascii.commented_header")

        # Plot the fitted speedup curves and write the parameters to the file
        fit_nthreads = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
        for mode in parameters:

            # Get the parameter values
            p = parameters[mode].p
            a = parameters[mode].a
            b = parameters[mode].b
            c = parameters[mode].c

            # Show the fitted curve
            if plot_fit:

                # Calculate the fitted speedups
                fit_speedups = [modAmdahl(n, p, a, b, c) for n in fit_nthreads]

                # Add the plot
                plt.plot(fit_nthreads, fit_speedups)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        plt.ylim(ticks[0], ticks[-1])
        plt.grid(True)

        # Plot a line that denotes linear scaling (speedup = nthreads)
        plt.plot(ticks, ticks, linestyle='--')

        # Add axis labels and a legend
        plt.xlabel("Number of processors $n$", fontsize='large')
        plt.ylabel(phase_labels[phase] + " speedup $S$", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Speedup of the " + phase_labels[phase].lower() + " for " + self.system_name)

        # Save the figure
        plt.savefig(file_path)
        plt.close()

    # -----------------------------------------------------------------

    def plot_efficiencies(self, phase, file_path, size=(12, 8), plot_fit=True):

        """
        This function creates a PDF plot showing the efficiency as a function of the number of threads.
        Efficiency is defined as T(1)/T(N)/N. It is a dimensionless quantity <= 1.
        The function takes the following (optional) arguments:
        :param phase:
        :param file_path:
        :param size:
        :param plot_fit:
        :return:
        """

        # Inform the user of the fact that the efficiencies are being calculated and plotted
        self.log.info("Calculating and plotting the efficiencies for the " + phase_names[phase] + "...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=size)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Get the serial runtime (and error) for this phase (create a Quantity object)
        serial_time = self.serial[phase].time
        serial_error = self.serial[phase].error
        serial = Quantity(serial_time, serial_error)

        # Loop over the different parallelization modes (the different curves)
        for mode in self.data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.data[phase][mode].processor_counts
            times = self.data[phase][mode].times
            errors = self.data[phase][mode].errors

            # Calculate the efficiencies and the errors on the efficiencies
            efficiencies = []
            efficiency_errors = []
            for i in range(len(processor_counts)):

                # Create a quantity for the current runtime
                time = Quantity(times[i], errors[i])

                # Calculate the efficiency based on the current runtime and the serial runtime
                speedup = serial / time
                efficiency = speedup.value / processor_counts[i]
                efficiency_error = speedup.error / processor_counts[i]

                # Add the value and the propagated error of the efficiency to the appropriate lists
                efficiencies.append(efficiency)
                efficiency_errors.append(efficiency_error)

            # Plot the data points for this curve
            plt.errorbar(processor_counts, efficiencies, efficiency_errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

        # Use a logaritmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        plt.ylim(0, 1.1)
        plt.grid(True)

        # Add axis labels and a legend
        plt.xlabel("Number of processors $n$", fontsize='large')
        plt.ylabel(phase_labels[phase] + " efficiency $\epsilon$", fontsize='large')
        plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Efficiency of the " + phase_labels[phase].lower() + " for " + self.system_name)

        # Save the figure
        plt.savefig(file_path)
        plt.close()

    # -----------------------------------------------------------------

    def plot_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the scaling timeline...")

        # Loop over the different parallelization modes
        for mode in self.data["total"]:

            # Determine the path to the timeline plot file
            plot_file_path = os.path.join(self.output_path, "timeline_" + mode + ".pdf")

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
            create_timeline_plot(data, plot_file_path, nprocs_list, percentages=True, totals=True, unordered=True, numberofproc=True, cpu=True, title=title)

    # -----------------------------------------------------------------

    def plot_memory(self, figsize=(12,8)):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the memory scaling...")

        # Determine the file path for this plot
        file_path = os.path.join(self.output_path, "memory.pdf")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parallelization modes (the different curves)
        for mode in self.data["memory"]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.data["memory"][mode].processor_counts
            memories = self.data["memory"][mode].times
            errors = self.data["memory"][mode].errors

            # Plot the data points for this mode
            plt.errorbar(processor_counts, memories, errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        plt.grid(True)

        # Add axis labels and a legend
        plt.xlabel("Number of processors $n$", fontsize='large')
        plt.ylabel("Memory usage (GB)", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Memory scaling for " + self.system_name)

        # Save the figure
        plt.savefig(file_path)
        plt.close()

# -----------------------------------------------------------------

## This function defines Amdahl's law for the speedup
def Amdahl(n, p): return 1.0 / (1 - p + p / n)

# This function defines a modified version of Amdahl's law, which accounts for different kinds of overhead
def modAmdahl(n, p, a, b, c): return 1.0 / (1 - p + p / n + a + b * n + c * n**2)

# -----------------------------------------------------------------
