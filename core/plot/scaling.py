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
from ..basics.quantity import Quantity
from ..basics.map import Map
from .timeline import create_timeline_plot
from ..tools.logging import log
from ..tools import filesystem as fs
from ..basics.configurable import Configurable
from ..launch.timing import TimingTable
from ..launch.memory import MemoryTable
from ..extract.timeline import TimeLineExtractor
from ..simulation.discover import SimulationDiscoverer

# -----------------------------------------------------------------

phase_names = {"total": "total simulation", "setup": "simulation setup", "stellar": "stellar emission phase",
               "spectra": "calculation of dust emission spectra", "dust": "dust emission phase",
               "writing": "writing phase", "waiting": "waiting phases", "communication": "communication phases"}

phase_labels = {"total": "Total runtime", "setup": "Setup time", "stellar": "Stellar runtime",
                "spectra": "Runtime of dust spectra calculation", "dust": "Dust emission runtime",
                "writing": "Writing time", "waiting": "Waiting time", "communication": "Communication time"}

# -----------------------------------------------------------------

scaling_properties = ["runtime", "speedup", "efficiency", "CPU time", "memory", "memory gain", "total memory", "timeline"]
simulation_phases = ["total", "setup", "stellar", "spectra", "dust", "writing", "waiting", "communication", "intermediate"]

# -----------------------------------------------------------------

class BatchScalingPlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(BatchScalingPlotter, self).__init__(config)

        # The list of simulations
        self.simulations = []

        # The timing and memory tables
        self.timing = None
        self.memory = None

        # The data
        self.timing_data = None
        self.memory_data = None

        # Serial
        self.serial_timing = None
        self.serial_memory = None

    # -----------------------------------------------------------------

    @property
    def has_serial_timing(self):

        """
        This function ...
        :return:
        """

        return len(self.serial_timing) != 0 if self.serial_timing is not None else False

    # -----------------------------------------------------------------

    @property
    def has_serial_memory(self):

        """
        This function ...
        :return:
        """

        return len(self.serial_memory) != 0 if self.serial_memory is not None else False

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Prepare data into plottable format
        self.prepare()

        # 3. Plot
        self.plot()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def add_simulation(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Add the simulation to the list
        self.simulations.append(simulation)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchScalingPlotter, self).setup(**kwargs)

        # Timing or memory specified
        self.timing = kwargs.pop("timing", None)
        self.memory = kwargs.pop("memory", None)

        # If either extracted timing or memory information is not passed
        if self.timing is None or self.memory is None:

            # If simulations are passed
            if "simulations" in kwargs: self.simulations = kwargs.pop("simulations")

            # If simulations have been added
            elif len(self.simulations) > 0: pass

            # Load simulations from working directory if none have been added
            else: self.load_simulations()

            # Do extraction
            self.extract()

        # Check the coverage of the timing and memory data
        self.check_coverage()

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
        self.simulations = discoverer.simulations_single_ski

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the timing and memory information ...")

        extract_timing = self.timing is None
        extract_memory = self.memory is None

        # Initialize a timing table
        if extract_timing: self.timing = TimingTable.initialize()

        # Initialize a memory table
        if extract_memory: self.memory = MemoryTable.initialize()

        # Loop over the simulations
        for simulation in self.simulations:

            # Load the log file
            log_file = simulation.log_file

            # Load the ski file or parameters map
            ski = simulation.parameters()
            parameters = None
            if isinstance(ski, Map):
                parameters = ski
                ski = None

            # If timing has to be extracted
            if extract_timing:

                # Create a TimeLineExtractor instance
                extractor = TimeLineExtractor()

                # Run the timeline extractor
                timeline = extractor.run(simulation)

                # Add an entry to the timing table
                self.timing.add_from_simulation(simulation, ski, log_file, timeline, parameters=parameters)

            if extract_memory:

                # Add an entry to the memory table
                self.memory.add_from_simulation(simulation, ski, log_file, parameters=parameters)

    # -----------------------------------------------------------------

    def check_coverage(self):

        """
        This function ...
        :return:
        """

        # Check if the data spans multiple number of cores
        if self.config.hybridisation:
            if np.min(self.timing["Processes"]) == np.max(self.timing["Processes"]): raise RuntimeError("All runtimes are generated with the same number of processes, you cannot run with --hybridisation")
            if np.min(self.memory["Processes"]) == np.max(self.timing["Processes"]): raise RuntimeError("All memory data are generated with the same number of processes, you cannot run with --hybridisation")
        else:
            if np.min(self.timing["Cores"]) == np.max(self.timing["Cores"]): raise RuntimeError("All runtimes are generated on the same number of cores. Run with --hybridisation")
            if np.min(self.memory["Cores"]) == np.max(self.memory["Cores"]): raise RuntimeError("All memory data are generated on the same number of cores. Run with --hybridisation")

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the data for plotting ...")

        # Initialize a data structure to contain the performance scaling information in plottable format
        self.timing_data = defaultdict(lambda: defaultdict(lambda: Map({"processor_counts": [], "times": [], "errors": []})))

        # Initialize a data structure to contain the memory scaling information in plottable format
        self.memory_data = defaultdict(lambda: defaultdict(lambda: Map({"processor_counts": [], "memory": [], "errors": []})))

        # Create an attribute to store the serial runtimes (one core)
        self.serial_timing = defaultdict(lambda: Map({"time": None, "error": None}))

        # Create an attribute to store the serial memory (one process)
        self.serial_memory = defaultdict(lambda: Map({"memory": None, "error": None}))

        sigma_level = 1.0

        # Create a dictionary to store the runtimes for the different simulation phases of the simulations
        # performed on one processor
        serial_times = defaultdict(list)
        serial_memory = defaultdict(list)

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
        intermediate_times = defaultdict(lambda: defaultdict(list))

        total_memory = defaultdict(lambda: defaultdict(list))
        setup_memory =  defaultdict(lambda: defaultdict(list))
        stellar_memory = defaultdict(lambda: defaultdict(list))
        spectra_memory =  defaultdict(lambda: defaultdict(list))
        dust_memory =  defaultdict(lambda: defaultdict(list))
        writing_memory =  defaultdict(lambda: defaultdict(list))

        # Loop over the different entries in the timing table
        for i in range(len(self.timing)):

            # Get the number of processes and threads
            processes = self.timing["Processes"][i]
            threads_per_core = self.timing["Threads per core"][i]
            cores_per_process = int(self.timing["Cores"][i] / processes)
            threads = threads_per_core * cores_per_process
            processors = processes * threads

            data_parallel = self.timing["Data-parallel"][i]

            #if processes > 1:
            #    if threads > 1:
            #        scaling_mode = "hybrid-" + str(threads)
            #    else: scaling_mode = "mpi"
            #else: scaling_mode = "threads"

            if self.config.hybridisation: mode = str(self.timing["Cores"][i]) + " cores"
            else:

                if processes > 1:
                    if threads > 1:
                        if data_parallel: mode = "hybrid task+data (" + str(threads) + ")"
                        else: mode = "hybrid task (" + str(threads) + ")"
                    else:
                        if data_parallel: mode = "mpi task+data"
                        else: mode = "mpi task"
                else: mode = "multithreading"

            # If the number of processors is 1, add the runtimes for the different simulation phases to the
            # dictionary that contains the serial runtimes
            if (self.config.hybridisation and processes == 1) or (not self.config.hybridisation and processors == 1):

                serial_times["total"].append(self.timing["Total runtime"][i])
                serial_times["setup"].append(self.timing["Setup time"][i])
                serial_times["stellar"].append(self.timing["Stellar emission time"][i])
                serial_times["spectra"].append(self.timing["Spectra calculation time"][i])
                serial_times["dust"].append(self.timing["Dust emission time"][i])
                serial_times["writing"].append(self.timing["Writing time"][i])
                serial_times["waiting"].append(self.timing["Waiting time"][i])
                serial_times["communication"].append(self.timing["Communication time"][i])

            # Number of processes = 1: equivalent to 'serial' in terms of memory consumption
            if processes == 1:

                serial_memory["total"].append(self.memory["Total peak memory"][i])
                serial_memory["setup"].append(self.memory["Setup peak memory"][i])
                serial_memory["stellar"].append(self.memory["Stellar emission peak memory"][i])
                serial_memory["spectra"].append(self.memory["Spectra calculation peak memory"][i])
                serial_memory["dust"].append(self.memory["Dust emission peak memory"][i])
                serial_memory["writing"].append(self.memory["Writing peak memory"][i])

            processes_or_processors = processes if self.config.hybridisation else processors

            # Add the processor count for this parallelization mode
            modes[mode].add(processes_or_processors)

            # Fill in the runtimes and memory usage at the appropriate place in the dictionaries
            total_times[mode][processes_or_processors].append(self.timing["Total runtime"][i])
            setup_times[mode][processes_or_processors].append(self.timing["Setup time"][i])
            stellar_times[mode][processes_or_processors].append(self.timing["Stellar emission time"][i])
            spectra_times[mode][processes_or_processors].append(self.timing["Spectra calculation time"][i])
            dust_times[mode][processes_or_processors].append(self.timing["Dust emission time"][i])
            writing_times[mode][processes_or_processors].append(self.timing["Writing time"][i])
            waiting_times[mode][processes_or_processors].append(self.timing["Waiting time"][i])
            communication_times[mode][processes_or_processors].append(self.timing["Communication time"][i])
            intermediate_times[mode][processes_or_processors].append(self.timing["Intermediate time"][i])

            # Fill in the memory usage at the appropriate place in the dictionaries
            total_memory[mode][processes_or_processors].append(self.memory["Total peak memory"][i])
            setup_memory[mode][processes_or_processors].append(self.memory["Setup peak memory"][i])
            stellar_memory[mode][processes_or_processors].append(self.memory["Stellar emission peak memory"][i])
            spectra_memory[mode][processes_or_processors].append(self.memory["Spectra calculation peak memory"][i])
            dust_memory[mode][processes_or_processors].append(self.memory["Dust emission peak memory"][i])
            writing_memory[mode][processes_or_processors].append(self.memory["Writing peak memory"][i])

        # Average the serial runtimes, loop over each phase
        for phase in serial_times:

            self.serial_timing[phase].time = np.mean(serial_times[phase])
            self.serial_timing[phase].error = sigma_level * np.std(serial_times[phase])

        # Average the serial memory usages, loop over each phase
        for phase in serial_memory:

            self.serial_memory[phase].memory = np.mean(serial_memory[phase])
            self.serial_memory[phase].error = sigma_level * np.std(serial_memory[phase])

        # Loop over all encountered parallelization modes
        for mode in modes:

            # Loop over all processor counts encountered for this mode
            for processors in modes[mode]:

                ## TIMING

                # Average the runtimes for the different simulation phases and the memory usage for the different
                # runs for a certain parallelization mode and number of processors
                self.timing_data["total"][mode].processor_counts.append(processors)
                self.timing_data["total"][mode].times.append(np.mean(total_times[mode][processors]))
                self.timing_data["total"][mode].errors.append(sigma_level * np.std(total_times[mode][processors]))

                self.timing_data["setup"][mode].processor_counts.append(processors)
                self.timing_data["setup"][mode].times.append(np.mean(setup_times[mode][processors]))
                self.timing_data["setup"][mode].errors.append(sigma_level * np.std(setup_times[mode][processors]))

                self.timing_data["stellar"][mode].processor_counts.append(processors)
                self.timing_data["stellar"][mode].times.append(np.mean(stellar_times[mode][processors]))
                self.timing_data["stellar"][mode].errors.append(sigma_level * np.std(stellar_times[mode][processors]))

                self.timing_data["spectra"][mode].processor_counts.append(processors)
                self.timing_data["spectra"][mode].times.append(np.mean(spectra_times[mode][processors]))
                self.timing_data["spectra"][mode].errors.append(sigma_level * np.std(spectra_times[mode][processors]))

                self.timing_data["dust"][mode].processor_counts.append(processors)
                self.timing_data["dust"][mode].times.append(np.mean(dust_times[mode][processors]))
                self.timing_data["dust"][mode].errors.append(sigma_level * np.std(dust_times[mode][processors]))

                self.timing_data["writing"][mode].processor_counts.append(processors)
                self.timing_data["writing"][mode].times.append(np.mean(writing_times[mode][processors]))
                self.timing_data["writing"][mode].errors.append(sigma_level * np.std(writing_times[mode][processors]))

                self.timing_data["waiting"][mode].processor_counts.append(processors)
                self.timing_data["waiting"][mode].times.append(np.mean(waiting_times[mode][processors]))
                self.timing_data["waiting"][mode].errors.append(sigma_level * np.std(waiting_times[mode][processors]))

                self.timing_data["communication"][mode].processor_counts.append(processors)
                self.timing_data["communication"][mode].times.append(np.mean(communication_times[mode][processors]))
                self.timing_data["communication"][mode].errors.append(sigma_level * np.std(communication_times[mode][processors]))

                self.timing_data["intermediate"][mode].processor_counts.append(processors)
                self.timing_data["intermediate"][mode].times.append(np.mean(intermediate_times[mode][processors]))
                self.timing_data["intermediate"][mode].errors.append(sigma_level * np.std(intermediate_times[mode][processors]))

                ## MEMORY

                self.memory_data["total"][mode].processor_counts.append(processors)
                self.memory_data["total"][mode].memory.append(np.mean(total_memory[mode][processors]))
                self.memory_data["total"][mode].errors.append(sigma_level * np.std(total_memory[mode][processors]))

                self.memory_data["setup"][mode].processor_counts.append(processors)
                self.memory_data["setup"][mode].memory.append(np.mean(setup_memory[mode][processors]))
                self.memory_data["setup"][mode].errors.append(sigma_level * np.std(setup_memory[mode][processors]))

                self.memory_data["stellar"][mode].processor_counts.append(processors)
                self.memory_data["stellar"][mode].memory.append(np.mean(stellar_memory[mode][processors]))
                self.memory_data["stellar"][mode].errors.append(sigma_level * np.std(stellar_memory[mode][processors]))

                self.memory_data["spectra"][mode].processor_counts.append(processors)
                self.memory_data["spectra"][mode].memory.append(np.mean(spectra_memory[mode][processors]))
                self.memory_data["spectra"][mode].errors.append(sigma_level * np.std(spectra_memory[mode][processors]))

                self.memory_data["dust"][mode].processor_counts.append(processors)
                self.memory_data["dust"][mode].memory.append(np.mean(dust_memory[mode][processors]))
                self.memory_data["dust"][mode].errors.append(sigma_level * np.std(dust_memory[mode][processors]))

                self.memory_data["writing"][mode].processor_counts.append(processors)
                self.memory_data["writing"][mode].memory.append(np.mean(writing_memory[mode][processors]))
                self.memory_data["writing"][mode].errors.append(sigma_level * np.std(writing_memory[mode][processors]))

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Runtime
        if "runtime" in self.config.properties: self.plot_runtimes()

        # Speedup
        if "speedup" in self.config.properties: self.plot_speedups()

        # Efficiency
        if "efficiency" in self.config.properties: self.plot_efficiencies()

        # CPU time
        if "CPU time" in self.config.properties: self.plot_cpu_times()

        # Memory
        if "memory" in self.config.properties: self.plot_memory()

        # Memory gain
        if "memory gain" in self.config.properties: self.plot_memory_gain()

        # Total memory
        if "total memory" in self.config.properties: self.plot_total_memory()

        # Timeline
        if "timeline" in self.config.properties: self.plot_timeline()

    # -----------------------------------------------------------------

    def plot_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Plot
            self.plot_runtimes_phase(phase)

    # -----------------------------------------------------------------

    def plot_speedups(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the speedups ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Plot
            self.plot_speedups_phase(phase)

    # -----------------------------------------------------------------

    def plot_efficiencies(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the efficiencies ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Plot
            self.plot_efficiencies_phase(phase)

    # -----------------------------------------------------------------

    def plot_cpu_times(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the CPU times ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Plot
            self.plot_cpu_times_phase(phase)

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the memory usage ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Plot
            self.plot_memory_phase(phase)

    # -----------------------------------------------------------------

    def plot_memory_gain(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the memory gain ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Plot
            self.plot_memory_gain_phase(phase)

    # -----------------------------------------------------------------

    def plot_total_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the total memory usage ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Plot
            self.plot_total_memory_phase(phase)

    # -----------------------------------------------------------------

    def plot_runtimes_phase(self, phase, file_path=None, figsize=(12, 8)):

        """
        This function ...
        :param phase:
        :param file_path:
        :param figsize:
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes for the " + phase_names[phase] + "...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parallelization modes (the different curves)
        for mode in self.timing_data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.timing_data[phase][mode].processor_counts
            times = self.timing_data[phase][mode].times
            errors = self.timing_data[phase][mode].errors

            # Sort the lists
            processor_counts, times, errors = sort_lists(processor_counts, times, errors)

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
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize='large')
        plt.ylabel(phase_labels[phase] + " T (s)", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Scaling of the " + phase_labels[phase].lower())

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_speedups_phase(self, phase, file_path=None, figsize=(12, 8), plot_fit=True):

        """
        This function ...
        :param phase:
        :param file_path:
        :param figsize:
        :param plot_fit:
        :return:
        """

        # Inform the user of the fact that the speedups are being calculated and plotted
        log.info("Plotting the speedups for the " + phase_names[phase] + "...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Create a dictionary that stores the fitted parameters for each different mode
        parameters = dict()

        # Get the serial runtime (and error) for this phase (create a Quantity object)
        serial_time = self.serial_timing[phase].time
        serial_error = self.serial_timing[phase].error
        serial = Quantity(serial_time, serial_error)

        # Loop over the different parallelization modes (the different curves)
        for mode in self.timing_data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.timing_data[phase][mode].processor_counts
            times = self.timing_data[phase][mode].times
            errors = self.timing_data[phase][mode].errors

            # Sort the lists
            processor_counts, times, errors = sort_lists(processor_counts, times, errors)

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

            if not self.config.hybridisation:

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

        if not self.config.hybridisation:

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

            if file_path is not None:

                # Create a data file to contain the fitted parameters
                directory = os.path.dirname(file_path)
                parameter_file_path = fs.join(directory, "parameters.dat")

                # Create the parameters table and write to file
                data = [p_list, p_error_list, a_list, a_error_list, b_list, b_error_list, c_list, c_error_list]
                names = ["Parallel fraction p", "Error on p", "Parameter a", "Error on a", "Parameter b", "Error on b", "Parameter c", "Error on c"]
                table = Table(data=data, names=names)
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
        if not self.config.hybridisation: ax.set_yticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        if not self.config.hybridisation: plt.ylim(ticks[0], ticks[-1])
        plt.grid(True)

        # Plot a line that denotes linear scaling (speedup = nthreads)
        if not self.config.hybridisation: plt.plot(ticks, ticks, linestyle='--')

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize='large')
        plt.ylabel(phase_labels[phase] + " speedup $S$", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Speedup of the " + phase_labels[phase].lower())

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_efficiencies_phase(self, phase, file_path=None, figsize=(12,8), plot_fit=True):

        """
        This function creates a PDF plot showing the efficiency as a function of the number of threads.
        Efficiency is defined as T(1)/T(N)/N. It is a dimensionless quantity <= 1.
        The function takes the following (optional) arguments:
        :param phase:
        :param file_path:
        :param figsize:
        :param plot_fit:
        :return:
        """

        # Inform the user of the fact that the efficiencies are being calculated and plotted
        log.info("Calculating and plotting the efficiencies for the " + phase_names[phase] + "...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Get the serial runtime (and error) for this phase (create a Quantity object)
        serial_time = self.serial_timing[phase].time
        serial_error = self.serial_timing[phase].error
        serial = Quantity(serial_time, serial_error)

        # Loop over the different parallelization modes (the different curves)
        for mode in self.timing_data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.timing_data[phase][mode].processor_counts
            times = self.timing_data[phase][mode].times
            errors = self.timing_data[phase][mode].errors

            # Sort the lists
            processor_counts, times, errors = sort_lists(processor_counts, times, errors)

            # Get array of number of used cores
            if self.config.hybridisation: ncores = np.ones(len(processor_counts)) * int(mode.split(" cores")[0])
            else: ncores = processor_counts

            # Calculate the efficiencies and the errors on the efficiencies
            efficiencies = []
            efficiency_errors = []
            for i in range(len(processor_counts)):

                # Create a quantity for the current runtime
                time = Quantity(times[i], errors[i])

                # Calculate the efficiency based on the current runtime and the serial runtime
                speedup = serial / time
                efficiency = speedup.value / ncores[i]
                efficiency_error = speedup.error / ncores[i]

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
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize="large")
        plt.ylabel(phase_labels[phase] + " efficiency $\epsilon$", fontsize='large')
        plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Efficiency of the " + phase_labels[phase].lower())

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_cpu_times_phase(self, phase, file_path=None, figsize=(12,8)):

        """
        This function ...
        :param phase:
        :param file_path:
        :param figsize:
        :return:
        """

        # Inform the user
        log.info("Plotting the CPU times for the " + phase_names[phase] + "...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parallelization modes (the different curves)
        for mode in self.timing_data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.timing_data[phase][mode].processor_counts
            times = self.timing_data[phase][mode].times
            errors = self.timing_data[phase][mode].errors

            # Sort the lists
            processor_counts, times, errors = sort_lists(processor_counts, times, errors)

            # Create numpy arrays
            processor_counts = np.array(processor_counts)
            times = np.array(times)
            errors = np.array(errors)

            # Get array of number of used cores
            if self.config.hybridisation: ncores = np.ones(len(processor_counts)) * int(mode.split(" cores")[0])
            else: ncores = processor_counts

            # Get list of process count
            #processes = nprocesses_from_mode(mode, processor_counts)

            # Multiply to get total
            times *= ncores
            errors *= ncores

            # Plot the data points for this mode
            plt.errorbar(processor_counts, times, errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale("log")
        plt.yscale("log")

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
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize='large')
        plt.ylabel(phase_labels[phase] + " T (s)", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Scaling of the total CPU time of the " + phase_labels[phase].lower())

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_memory_phase(self, phase, file_path=None, figsize=(12,8)):

        """
        This function ...
        :param phase:
        :param file_path:
        :param figsize:
        :return:
        """

        # Inform the user
        log.info("Plotting the memory scaling...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parallelization modes (the different curves)
        for mode in self.memory_data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.memory_data[phase][mode].processor_counts
            memories = self.memory_data[phase][mode].memory
            errors = self.memory_data[phase][mode].errors

            # Sort the lists
            processor_counts, memories, errors = sort_lists(processor_counts, memories, errors)

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
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize='large')
        plt.ylabel("Memory usage per process (GB)", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Scaling of the memory usage (per process) of the " + phase + " phase")

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_memory_gain_phase(self, phase, file_path=None, figsize=(12,8)):

        """
        This function ...
        :param phase:
        :param file_path:
        :param figsize:
        :return:
        """

        # Inform the user
        log.info("Plotting the memory scaling...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Get the serial memory (and error) for this phase (create a Quantity object)
        serial_memory = self.serial_memory[phase].memory
        serial_error = self.serial_memory[phase].error
        serial = Quantity(serial_memory, serial_error)

        # Loop over the different parallelization modes (the different curves)
        for mode in self.memory_data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.memory_data[phase][mode].processor_counts
            memories = self.memory_data[phase][mode].memory
            errors = self.memory_data[phase][mode].errors

            # Sort the lists
            processor_counts, memories, errors = sort_lists(processor_counts, memories, errors)

            # Calculate the gains and the errors on the gains
            gains = []
            gain_errors = []
            for i in range(len(processor_counts)):

                # Create a quantity for the current memory usage
                memory = Quantity(memories[i], errors[i])

                # Calculate the efficiency based on the current memory usage and the serial memory usage
                gain = serial / memory

                # Add the value and the propagated error of the gain to the appropriate lists
                gains.append(gain.value)
                gain_errors.append(gain.error)

            # Plot the data points for this mode
            plt.errorbar(processor_counts, gains, gain_errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale("log")
        plt.yscale("log")

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
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize='large')
        plt.ylabel("Memory gain", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Scaling of the memory gain (serial memory usage per process / memory usage per process) of the " + phase_labels[phase].lower())

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_total_memory_phase(self, phase, file_path=None, figsize=(12,8)):

        """
        This function ...
        :param phase:
        :param file_path:
        :param figsize:
        :return:
        """

        # Inform the user
        log.info("Plotting the total memory scaling (all processes combined)...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parallelization modes (the different curves)
        for mode in self.memory_data[phase]:

            # Get the list of processor counts, runtimes and errors
            processor_counts = self.memory_data[phase][mode].processor_counts
            memories = self.memory_data[phase][mode].memory
            errors = self.memory_data[phase][mode].errors

            # Sort the lists
            processor_counts, memories, errors = sort_lists(processor_counts, memories, errors)

            # Create arrays
            processor_counts = np.array(processor_counts)
            memories = np.array(memories)
            errors = np.array(errors)

            # Get list of process count
            if self.config.hybridisation: processes = processor_counts
            else: processes = nprocesses_from_mode(mode, processor_counts)

            # Multiply the memory usage for each processor count with the corresponding number of processes (to get the total)
            memories *= processes
            errors *= processes

            # Plot the data points for this mode
            plt.errorbar(processor_counts, memories, errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale("log")
        plt.yscale("log")

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
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize='large')
        plt.ylabel("Total memory usage (all processes) (GB)", fontsize='large')
        plt.legend(title="Modes")

        # Set the plot title
        plt.title("Memory scaling for " + phase_labels[phase])

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_timeline(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the scaling timeline...")

        # Loop over the different parallelization modes
        for mode in self.timing_data["total"]:

            # Determine the path to the timeline plot file
            #plot_file_path = fs.join(self.output_path, "timeline_" + mode + ".pdf")
            plot_file_path = None

            # Initialize a data structure to contain the start times and endtimes of the different simulation phases,
            # for the different processor counts (data is indexed on the simulation phase)
            data = []
            nprocs_list = []

            # Loop over the different processor counts
            for j in range(len(self.timing_data["total"][mode].processor_counts)):

                # Get the processor count
                if self.config.hybridisation:
                    processors = int(mode.split(" cores")[0])
                    processes = self.timing_data["total"][mode].processor_counts[j]
                else:
                    processors = self.timing_data["total"][mode].processor_counts[j]
                    processes = nprocesses_from_mode_single(mode, processors)

                # Get the average runtimes for the different phases corresponding to the current processor count
                setup_time = self.timing_data["setup"][mode].times[j] * processors
                stellar_time = self.timing_data["stellar"][mode].times[j] * processors
                spectra_time = self.timing_data["spectra"][mode].times[j] * processors
                dust_time = self.timing_data["dust"][mode].times[j] * processors
                writing_time = self.timing_data["writing"][mode].times[j] * processors
                waiting_time = self.timing_data["waiting"][mode].times[j] * processors
                communication_time = self.timing_data["communication"][mode].times[j] * processors

                # Add the process count
                nprocs_list.append(processes)

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
            title = "Scaling timeline"

            # Create the plot
            create_timeline_plot(data, nprocs_list, plot_file_path, percentages=True, totals=True, unordered=True, numberofproc=True, cpu=True, title=title)

    # -----------------------------------------------------------------

    def plot_dries(self):

        """
        This function ...
        :return:
        """

        dmpiTimes = np.loadtxt("dmpiTimes.dat", usecols=range(1, 8))
        smpiTimes = np.loadtxt("smpiTimes.dat", usecols=range(1, 8))
        tTimes = np.loadtxt("tTimes.dat", usecols=range(1, 8))

        def sortByCores(times):
            times[:] = times[times[:, 6].argsort()]

        sortByCores(tTimes)
        sortByCores(smpiTimes)
        sortByCores(dmpiTimes)

        singlecore = tTimes[0, 5]

        plt.figure(figsize=(8, 5))

        def plot(times, key):
            cores = times[:, 6]
            values = singlecore / times[:, 5] / cores
            plt.semilogx(cores, values, label=key, basex=2, marker='D')
            ax = plt.gca()
            for axis in [ax.xaxis, ax.yaxis]:
                axis.set_major_formatter(ScalarFormatter())

        plot(tTimes, "Threads")
        plot(smpiTimes, "MPI task-based")
        plot(dmpiTimes, "MPI data-parallel")
        plt.plot([1, 16], [1, 1], color='k', linewidth=3, linestyle='dashed')
        plt.title("Single node scaling", y=1.025)
        plt.xlabel("Number of cores used")
        plt.ylabel("Efficiency")
        plt.grid()
        plt.legend(loc="best")
        #plt.savefig("efficiency.pdf")

        plt.show()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def nprocesses_from_mode_single(mode, nprocessors):

    """
    This function ...
    :param mode:
    :param nprocessors:
    :return:
    """

    # Get number of processes
    if mode == "multithreading": processes = 1
    elif "mpi" in mode: processes = nprocessors
    elif "hybrid" in mode:
        threads = int(mode.split("(")[1].split(")")[0])
        processes = nprocessors / threads
    else: raise ValueError("Invalid mode: " + mode)

    # Return the number of processes
    return processes

# -----------------------------------------------------------------

def nprocesses_from_mode(mode, processor_counts):

    """
    This function ...
    :param mode:
    :param processor_counts:
    :return:
    """

    # Get number of processes
    if mode == "multithreading": processes = np.ones(len(processor_counts))
    elif "mpi" in mode: processes = processor_counts
    elif "hybrid" in mode:
        threads = int(mode.split("(")[1].split(")")[0])
        processes = processor_counts / threads
    else: raise ValueError("Invalid mode: " + mode)

    # Return the list of proces count
    return processes

# -----------------------------------------------------------------

## This function defines Amdahl's law for the speedup
def Amdahl(n, p): return 1.0 / (1 - p + p / n)

# This function defines a modified version of Amdahl's law, which accounts for different kinds of overhead
def modAmdahl(n, p, a, b, c): return 1.0 / (1 - p + p / n + a + b * n + c * n**2)

# -----------------------------------------------------------------

def sort_lists(*args): return [list(t) for t in zip(*sorted(zip(*args)))]

# -----------------------------------------------------------------
