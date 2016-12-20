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
import matplotlib
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from collections import defaultdict
from collections import Callable
from types import FunctionType
from string import ascii_lowercase

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
from ..simulation.discover import SimulationDiscoverer, comparison_parameters_from_ski
from ..basics.range import RealRange
from ..tools import tables

# -----------------------------------------------------------------

alphabet = list(ascii_lowercase)

# -----------------------------------------------------------------

phase_names = {"total": "total simulation", "setup": "simulation setup", "stellar": "stellar emission phase",
               "spectra": "calculation of dust emission spectra", "dust": "dust emission phase",
               "writing": "writing phase", "waiting": "waiting phases", "communication": "communication phases",
               "dust densities communication": "commmunication of the dust densities",
               "stellar absorption communication": "communication of the absorbed stellar luminosities",
               "dust absorption communication": "communication of the absorbed dust luminosities",
               "emission spectra communication": "communication of the dust emission spectra",
               "instruments communication": "communication of the instrument data",
               "intermediate": "intermediate procedures"}

phase_labels = {"total": "Total runtime", "setup": "Setup time", "stellar": "Stellar runtime",
                "spectra": "Runtime of dust spectra calculation", "dust": "Dust emission runtime",
                "writing": "Writing time", "waiting": "Waiting time", "communication": "Communication time",
                "dust densities communication": "Dust densities communication time",
                "stellar absorption communication": "Stellar absorption table communication time",
                "dust absorption communication": "Dust absorption table communication time",
                "emission spectra communication": "Emission spectra communication time",
                "instruments communication": "Instrument tables communication time",
                "intermediate": "Intermediate time"}

# -----------------------------------------------------------------

parallel_phases = ["total", "stellar", "spectra", "dust"]
overhead_phases = ["waiting", "communication"]
serial_phases = ["writing", "setup"]

# -----------------------------------------------------------------

# Scaling properties
scaling_properties = ["runtime", "speedup", "efficiency", "CPU-time", "memory", "memory-gain", "total-memory", "timeline"]

# Scaling phases
simulation_phases = ["total", "setup", "stellar", "spectra", "dust", "writing", "waiting", "communication", "intermediate"]

simulation_phases.append("dust densities communication")
simulation_phases.append("stellar absorption communication")
simulation_phases.append("dust absorption communication")
simulation_phases.append("emission spectra communication")
simulation_phases.append("instruments communication")

# Seperate types of properties
timing_properties = ["runtime", "speedup", "efficiency", "CPU-time", "timeline"]
memory_properties = ["memory", "memory-gain", "total-memory"]

# -----------------------------------------------------------------

# Communication things:
# dust densities communication
# stellar absorption communication
# dust absorption communication
# emission spectra communication
# instruments communication

# -----------------------------------------------------------------

pure_scaling_behaviour = dict()

# Setup
pure_scaling_behaviour["setup"] = [0]

# Stellar emission
pure_scaling_behaviour["stellar"] = [-1]

# Dust emission
pure_scaling_behaviour["dust"] = [-1]

# Spectra
pure_scaling_behaviour["spectra"] = [-1]

# Communication
pure_scaling_behaviour["dust densities communication"] = [0, 1, 2, "log"]
pure_scaling_behaviour["stellar absorption communication"] = [0, 1, 2, "log"]
pure_scaling_behaviour["dust absorption communication"] = [0, 1, 2, "log"]
pure_scaling_behaviour["emission spectra communication"] = [0, 1, 2, "log"]
pure_scaling_behaviour["instruments communication"] = [0, 1, 2, "log"]

# Waiting
pure_scaling_behaviour["waiting"] = [0, 1, 2]

# Writing
pure_scaling_behaviour["writing"] = [0]

# -----------------------------------------------------------------

composite_scaling_behaviour = dict()

# Communication
composite_scaling_behaviour["communication"] = ("dust densities communication", "stellar absorption communication",
                                                "dust absorption communication", "emission spectra communication",
                                                "instruments communication")

# Total simulation
composite_scaling_behaviour["total"] = ("setup", "stellar", "dust", "spectra", "communication", "waiting", "writing")

# -----------------------------------------------------------------

class ScalingPlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(ScalingPlotter, self).__init__(config)

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

        # Fit functions
        self.timing_fit_functions = defaultdict(dict)

        # Fit parameters
        self.timing_fit_parameters = defaultdict(dict)
        self.timing_fit_parameter_errors = defaultdict(dict)

        self.memory_fit_parameters = dict()

        # ...
        self.serial_timing_ncores = dict()
        self.serial_memory_ncores = dict()

        # Parameter set mapping for the simulations
        self.distinct_parameter_sets = []
        self.parameters_mapping = []

    # -----------------------------------------------------------------

    def has_serial_timing(self, phase):

        """
        This function ...
        :return:
        """

        return self.has_reference_timing and self.serial_timing_ncores[phase] == 1

    # -----------------------------------------------------------------

    @property
    def has_reference_timing(self):

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

        return self.has_reference_memory

    # -----------------------------------------------------------------

    @property
    def has_reference_memory(self):

        """
        This function ...
        :return:
        """

        return len(self.serial_memory) != 0 if self.serial_memory is not None else False

    # -----------------------------------------------------------------

    @property
    def needs_timing(self):

        """
        This function ...
        :return:
        """

        for property in self.config.properties:
            if property in timing_properties: return True
        return False

    # -----------------------------------------------------------------

    @property
    def needs_memory(self):

        """
        This function ...
        :return:
        """

        for property in self.config.properties:
            if property in memory_properties: return True
        return False

    # -----------------------------------------------------------------

    def has_timing_fit(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        return phase in self.timing_fit_parameters

    # -----------------------------------------------------------------

    def has_memory_fit(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        return phase in self.memory_fit_parameters

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

        # Do fitting
        if not self.config.hybridisation and self.config.fit: self.fit()

        # 3. Plot
        self.plot()

        # 4. Write
        if self.config.output is not None: self.write()

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
        super(ScalingPlotter, self).setup(**kwargs)

        # Check for 'all' flag
        if self.config.all:
            self.config.properties = scaling_properties
            self.config.phases = simulation_phases

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

        # Check the ski file parameters
        if self.has_multiski: self.check_parameters()
        else:
            parameters = comparison_parameters_from_ski(self.simulations[0].ski_path, self.simulations[0].input_path)
            self.distinct_parameter_sets.append(parameters)

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

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the timing and memory information ...")

        # Check whether timing and/or memory tables have to be extracted
        extract_timing = self.timing is None
        extract_memory = self.memory is None

        # Initialize a timing table
        if extract_timing: self.timing = TimingTable.initialize()

        # Initialize a memory table
        if extract_memory: self.memory = MemoryTable.initialize()

        # Loop over the simulations
        for simulation in self.simulations:

            # Debugging
            log.debug("Extracting timeline and memory information from simulation in '" + simulation.output_path + "' ...")

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

        # Debugging
        log.debug("Checking the coverage of the data ...")

        # Check if the data spans multiple number of cores
        if self.config.hybridisation:
            if np.min(self.timing["Processes"]) == np.max(self.timing["Processes"]): raise RuntimeError("All runtimes are generated with the same number of processes, you cannot run with --hybridisation")
            if np.min(self.memory["Processes"]) == np.max(self.timing["Processes"]): raise RuntimeError("All memory data are generated with the same number of processes, you cannot run with --hybridisation")
        else:
            if np.min(self.timing["Cores"]) == np.max(self.timing["Cores"]): raise RuntimeError("All runtimes are generated on the same number of cores. Run with --hybridisation")
            if np.min(self.memory["Cores"]) == np.max(self.memory["Cores"]): raise RuntimeError("All memory data are generated on the same number of cores. Run with --hybridisation")

    # -----------------------------------------------------------------

    def check_parameters(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Checking simulation parameters ...")

        # Loop over the simulations
        for simulation in self.simulations:

            # Get the parameters
            parameters = comparison_parameters_from_ski(simulation.ski_path, simulation.input_path)

            # Get index
            index = None
            for i, dps in enumerate(self.distinct_parameter_sets):
                if dps == parameters:
                    index = i
                    break

            if index is None:
                self.distinct_parameter_sets.append(parameters)
                index = len(parameters) - 1

            # Set index for the mapping
            self.parameters_mapping.append(index)

    # -----------------------------------------------------------------

    @property
    def has_multiski(self):

        """
        This function ...
        :return:
        """

        ski_path = None

        for simulation in self.simulations:

            # Cannot determine whether multiski
            if simulation.ski_path is None: return None

            if ski_path is None: ski_path = simulation.ski_path
            elif ski_path != simulation.ski_path: return True

        return True

    # -----------------------------------------------------------------

    @property
    def has_multiparameters(self):

        """
        This function ...
        :return:
        """

        return len(self.distinct_parameter_sets) > 1

    # -----------------------------------------------------------------

    @property
    def different_parameters_timing(self):

        """
        This function ...
        :return:
        """

        return self.timing.different_parameters()

    # -----------------------------------------------------------------

    @property
    def different_parameters_memory(self):

        """
        This function ...
        :return:
        """

        return self.memory.different_parameters()

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the data for plotting ...")

        # Initialize a data structure to contain the performance scaling information in plottable format
        self.timing_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: Map({"processor_counts": [], "times": [], "errors": []}))))

        # Initialize a data structure to contain the memory scaling information in plottable format
        self.memory_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: Map({"processor_counts": [], "memory": [], "errors": []}))))

        # Create an attribute to store the serial runtimes (one core)
        self.serial_timing = defaultdict(lambda: defaultdict(lambda: Map({"time": None, "error": None})))

        # Create an attribute to store the serial memory (one process)
        self.serial_memory = defaultdict(lambda: defaultdict(lambda: Map({"memory": None, "error": None})))

        # Create a dictionary to store the runtimes for the different simulation phases of the simulations
        # performed on one processor
        serial_times = defaultdict(lambda: defaultdict(list))
        serial_memory = defaultdict(lambda: defaultdict(list))

        # Keep track of the different processor counts encountered for the different parameter sets and parallelization modes
        parameters_modes_processor_counts_dict = defaultdict(lambda: defaultdict(set))

        # Create dictionaries to contain the data before it is averaged over the different simulations
        total_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        setup_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        stellar_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        spectra_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        dust_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        writing_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        waiting_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        communication_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        densities_communication_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        stellarabsorption_communication_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        dustabsorption_communication_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        emission_communication_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        instruments_communication_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        intermediate_times = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        # Create dictionaries for the memory data
        total_memorys = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        setup_memorys = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        stellar_memorys = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        spectra_memorys = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        dust_memorys = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        writing_memorys = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        # Number of cores for the serial memory runs
        serial_memory_ncores = dict()

        # Loop over the different entries in the timing table
        for i in range(len(self.timing)):

            # Get relevant simulation properties
            relevant_parameters = dict()
            for parameter in self.different_parameters_timing:
                relevant_parameters[parameter] = self.timing[parameter][i]
            if len(relevant_parameters) == 0: relevant_parameters = "empty"

            # Get the number of processes and threads
            processes = self.timing["Processes"][i]
            threads_per_core = self.timing["Threads per core"][i]
            cores_per_process = int(self.timing["Cores"][i] / processes)
            threads = threads_per_core * cores_per_process
            processors = processes * threads

            data_parallel = self.timing["Data-parallel"][i]

            # Determine the parallelization mode
            if self.config.hybridisation: mode = str(self.timing["Cores"][i]) + " cores"
            else: mode = get_parallelization_mode(processes, threads, data_parallel)

            #if len(relevant_parameters) > 0:
            #    mode_with_parameters = mode + " ["
            #    for parameter in relevant_parameters: mode += parameter + ": " + str(relevant_parameters[parameter]) + ", "
            #    mode_with_parameters += "]"
            #else: mode_with_parameters = mode
            #mode = mode_with_parameters

            # Get the times
            total_time = self.timing["Total runtime"][i]
            setup_time = self.timing["Setup time"][i]
            stellar_time = self.timing["Stellar emission time"][i]
            spectra_time = self.timing["Spectra calculation time"][i]
            dust_time = self.timing["Dust emission time"][i]
            writing_time = self.timing["Writing time"][i]
            waiting_time = self.timing["Waiting time"][i]
            communication_time = self.timing["Communication time"][i]
            densities_comm_time = self.timing["Dust densities communication time"][i]
            stellarabsorption_comm_time = self.timing["Stellar absorption communication time"][i]
            dustabsorption_comm_time = self.timing["Dust absorption communication time"][i]
            emission_comm_time = self.timing["Emission spectra communication time"][i]
            instruments_comm_time = self.timing["Instruments communication time"][i]
            intermediate_time = self.timing["Intermediate time"][i]

            # Get the memory usages
            total_memory = self.memory["Total peak memory"][i]
            setup_memory = self.memory["Setup peak memory"][i]
            stellar_memory = self.memory["Stellar emission peak memory"][i]
            spectra_memory = self.memory["Spectra calculation peak memory"][i]
            dust_memory = self.memory["Dust emission peak memory"][i]
            writing_memory = self.memory["Writing peak memory"][i]

            # If the number of processors is 1, add the runtimes for the different simulation phases to the
            # dictionary that contains the serial runtimes
            if (self.config.hybridisation and processes == 1) or (not self.config.hybridisation and processors == 1):

                serial_times[relevant_parameters]["total"].append(total_time)
                serial_times[relevant_parameters]["setup"].append(setup_time)
                serial_times[relevant_parameters]["stellar"].append(stellar_time)
                serial_times[relevant_parameters]["spectra"].append(spectra_time)
                serial_times[relevant_parameters]["dust"].append(dust_time)
                serial_times[relevant_parameters]["writing"].append(writing_time)
                serial_times[relevant_parameters]["waiting"].append(waiting_time)
                serial_times[relevant_parameters]["communication"].append(communication_time)

            # Number of processes = 1: equivalent to 'serial' in terms of memory consumption
            if processes == 1:

                serial_memory[relevant_parameters]["total"].append(total_memory)
                serial_memory[relevant_parameters]["setup"].append(setup_memory)
                serial_memory[relevant_parameters]["stellar"].append(stellar_memory)
                serial_memory[relevant_parameters]["spectra"].append(spectra_memory)
                serial_memory[relevant_parameters]["dust"].append(dust_memory)
                serial_memory[relevant_parameters]["writing"].append(writing_memory)

                serial_memory_ncores[relevant_parameters] = processors

            processes_or_processors = processes if self.config.hybridisation else processors

            # Add the processor count for this parallelization mode
            #modes[mode].add(processes_or_processors)

            parameters_modes_processor_counts_dict[relevant_parameters][mode].add(processes_or_processors)

            # Fill in the runtimes and memory usage at the appropriate place in the dictionaries
            total_times[relevant_parameters][mode][processes_or_processors].append(total_time)
            setup_times[relevant_parameters][mode][processes_or_processors].append(setup_time)
            stellar_times[relevant_parameters][mode][processes_or_processors].append(stellar_time)
            spectra_times[relevant_parameters][mode][processes_or_processors].append(spectra_time)
            dust_times[relevant_parameters][mode][processes_or_processors].append(dust_time)
            writing_times[relevant_parameters][mode][processes_or_processors].append(writing_time)
            waiting_times[relevant_parameters][mode][processes_or_processors].append(waiting_time)
            communication_times[relevant_parameters][mode][processes_or_processors].append(communication_time)
            densities_communication_times[relevant_parameters][mode][processes_or_processors].append(densities_comm_time)
            stellarabsorption_communication_times[relevant_parameters][mode][processes_or_processors].append(stellarabsorption_comm_time)
            dustabsorption_communication_times[relevant_parameters][mode][processes_or_processors].append(dustabsorption_comm_time)
            emission_communication_times[relevant_parameters][mode][processes_or_processors].append(emission_comm_time)
            instruments_communication_times[relevant_parameters][mode][processes_or_processors].append(instruments_comm_time)
            intermediate_times[relevant_parameters][mode][processes_or_processors].append(intermediate_time)

            # Fill in the memory usage at the appropriate place in the dictionaries
            total_memorys[relevant_parameters][mode][processes_or_processors].append(total_memory)
            setup_memorys[relevant_parameters][mode][processes_or_processors].append(setup_memory)
            stellar_memorys[relevant_parameters][mode][processes_or_processors].append(stellar_memory)
            spectra_memorys[relevant_parameters][mode][processes_or_processors].append(spectra_memory)
            dust_memorys[relevant_parameters][mode][processes_or_processors].append(dust_memory)
            writing_memorys[relevant_parameters][mode][processes_or_processors].append(writing_memory)

        # Average the serial runtimes, loop over each parameter set and phase
        for parameter_set in serial_times:
            for phase in serial_times[parameter_set]:

                self.serial_timing[parameter_set][phase].time = np.mean(serial_times[parameter_set][phase])
                self.serial_timing[parameter_set][phase].error = self.config.sigma_level * np.std(serial_times[parameter_set][phase])
                self.serial_timing_ncores[parameter_set][phase] = 1

        # Average the serial memory usages, loop over each phase
        for parameter_set in serial_memory:
            for phase in serial_memory[parameter_set]:

                self.serial_memory[parameter_set][phase].memory = np.mean(serial_memory[parameter_set][phase])
                self.serial_memory[parameter_set][phase].error = self.config.sigma_level * np.std(serial_memory[parameter_set][phase])
                self.serial_memory_ncores[parameter_set][phase] = serial_memory_ncores[parameter_set]

        # Loop over all parameter sets
        for parameter_set in parameters_modes_processor_counts_dict:

            # Loop over all encountered parallelization modes
            modes = parameters_modes_processor_counts_dict[parameter_set]
            for mode in modes:

                # Loop over all processor counts encountered for this mode
                for processors in modes[mode]:

                    ## TIMING

                    # Average the runtimes for the different simulation phases and the memory usage for the different
                    # runs for a certain parallelization mode and number of processors
                    self.timing_data["total"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["total"][parameter_set][mode].times.append(np.mean(total_times[parameter_set][mode][processors]))
                    self.timing_data["total"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(total_times[parameter_set][mode][processors]))

                    self.timing_data["setup"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["setup"][parameter_set][mode].times.append(np.mean(setup_times[parameter_set][mode][processors]))
                    self.timing_data["setup"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(setup_times[parameter_set][mode][processors]))

                    self.timing_data["stellar"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["stellar"][parameter_set][mode].times.append(np.mean(stellar_times[parameter_set][mode][processors]))
                    self.timing_data["stellar"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(stellar_times[parameter_set][mode][processors]))

                    self.timing_data["spectra"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["spectra"][parameter_set][mode].times.append(np.mean(spectra_times[parameter_set][mode][processors]))
                    self.timing_data["spectra"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(spectra_times[parameter_set][mode][processors]))

                    self.timing_data["dust"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["dust"][parameter_set][mode].times.append(np.mean(dust_times[parameter_set][mode][processors]))
                    self.timing_data["dust"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(dust_times[parameter_set][mode][processors]))

                    self.timing_data["writing"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["writing"][parameter_set][mode].times.append(np.mean(writing_times[parameter_set][mode][processors]))
                    self.timing_data["writing"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(writing_times[parameter_set][mode][processors]))

                    self.timing_data["waiting"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["waiting"][parameter_set][mode].times.append(np.mean(waiting_times[parameter_set][mode][processors]))
                    self.timing_data["waiting"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(waiting_times[parameter_set][mode][processors]))

                    self.timing_data["communication"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["communication"][parameter_set][mode].times.append(np.mean(communication_times[parameter_set][mode][processors]))
                    self.timing_data["communication"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(communication_times[parameter_set][mode][processors]))

                    self.timing_data["dust densities communication"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["dust densities communication"][parameter_set][mode].times.append(np.mean(densities_communication_times[parameter_set][mode][processors]))
                    self.timing_data["dust densities communication"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(densities_communication_times[parameter_set][mode][processors]))

                    self.timing_data["stellar absorption communication"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["stellar absorption communication"][parameter_set][mode].times.append(np.mean(stellarabsorption_communication_times[parameter_set][mode][processors]))
                    self.timing_data["stellar absorption communication"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(stellarabsorption_communication_times[parameter_set][mode][processors]))

                    self.timing_data["dust absorption communication"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["dust absorption communication"][parameter_set][mode].times.append(np.mean(dustabsorption_communication_times[parameter_set][mode][processors]))
                    self.timing_data["dust absorption communication"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(dustabsorption_communication_times[parameter_set][mode][processors]))

                    self.timing_data["emission spectra communication"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["emission spectra communication"][parameter_set][mode].times.append(np.mean(emission_communication_times[parameter_set][mode][processors]))
                    self.timing_data["emission spectra communication"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(emission_communication_times[parameter_set][mode][processors]))

                    self.timing_data["instruments communication"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["instruments communication"][parameter_set][mode].times.append(np.mean(instruments_communication_times[parameter_set][mode][processors]))
                    self.timing_data["instruments communication"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(instruments_communication_times[parameter_set][mode][processors]))

                    self.timing_data["intermediate"][parameter_set][mode].processor_counts.append(processors)
                    self.timing_data["intermediate"][parameter_set][mode].times.append(np.mean(intermediate_times[parameter_set][mode][processors]))
                    self.timing_data["intermediate"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(intermediate_times[parameter_set][mode][processors]))

                    ## MEMORY

                    self.memory_data["total"][parameter_set][mode].processor_counts.append(processors)
                    self.memory_data["total"][parameter_set][mode].memory.append(np.mean(total_memorys[parameter_set][mode][processors]))
                    self.memory_data["total"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(total_memorys[parameter_set][mode][processors]))

                    self.memory_data["setup"][parameter_set][mode].processor_counts.append(processors)
                    self.memory_data["setup"][parameter_set][mode].memory.append(np.mean(setup_memorys[parameter_set][mode][processors]))
                    self.memory_data["setup"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(setup_memorys[parameter_set][mode][processors]))

                    self.memory_data["stellar"][parameter_set][mode].processor_counts.append(processors)
                    self.memory_data["stellar"][parameter_set][mode].memory.append(np.mean(stellar_memorys[parameter_set][mode][processors]))
                    self.memory_data["stellar"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(stellar_memorys[parameter_set][mode][processors]))

                    self.memory_data["spectra"][parameter_set][mode].processor_counts.append(processors)
                    self.memory_data["spectra"][parameter_set][mode].memory.append(np.mean(spectra_memorys[parameter_set][mode][processors]))
                    self.memory_data["spectra"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(spectra_memorys[parameter_set][mode][processors]))

                    self.memory_data["dust"][parameter_set][mode].processor_counts.append(processors)
                    self.memory_data["dust"][parameter_set][mode].memory.append(np.mean(dust_memorys[parameter_set][mode][processors]))
                    self.memory_data["dust"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(dust_memorys[parameter_set][mode][processors]))

                    self.memory_data["writing"][parameter_set][mode].processor_counts.append(processors)
                    self.memory_data["writing"][parameter_set][mode].memory.append(np.mean(writing_memorys[parameter_set][mode][processors]))
                    self.memory_data["writing"][parameter_set][mode].errors.append(self.config.sigma_level * np.std(writing_memorys[parameter_set][mode][processors]))

        # Set missing serial timing data
        if self.needs_timing: self.set_missing_serial_timing(parameters_modes_processor_counts_dict)

        # Set missing serial memory data
        if self.needs_memory: self.set_missing_serial_memory(parameters_modes_processor_counts_dict)

    # -----------------------------------------------------------------

    def set_missing_serial_timing(self, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        # Check if serial data is found for the different parameters set
        for parameter_set in parameters_modes_processor_counts_dict:

            # Set missing serial timing through extrapolation of npackages
            if len(self.serial_timing[parameter_set]) == 0 and self.config.extrapolate_npackages:
                self.set_missing_serial_timing_npackages(parameter_set, parameters_modes_processor_counts_dict)

            # Set missing serial timing through extrapolation of ncores
            if len(self.serial_timing[parameter_set]) == 0 and self.config.extrapolate_ncores: pass

            # Set missing serial timing as timing that is not really serial ...
            if len(self.serial_timing[parameter_set]) == 0: self.set_missing_serial_timing_false(parameter_set)

    # -----------------------------------------------------------------

    def set_missing_serial_timing_npackages(self, parameter_set, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameter_set:
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        the_other_npackages = None
        the_other_parameter_set = None

        # Loop over the other parameter sets
        other_parameter_sets = parameters_modes_processor_counts_dict.keys()
        other_parameter_sets.remove(parameter_set)
        for other_parameter_set in other_parameter_sets:

            other_npackages = None

            for name in other_parameter_set:

                if name == "Packages":
                    other_npackages = other_parameter_set["Packages"]
                elif other_parameter_set[name] != parameter_set[name]:
                    break

            # If a break is not encountered
            else:

                # End the loop over the other parameter sets because we found a good one
                the_other_npackages = other_npackages
                the_other_parameter_set = other_parameter_set
                break

        # Only if a good match was found
        if the_other_parameter_set is None: return

        # Multiplication factor
        npackages_factor = float(the_other_npackages) / float(parameter_set["Packages"])

        the_mode = None
        the_serial_index = None

        # Check if a serial run is available for this parameter set (mode is not important)
        for mode in parameters_modes_processor_counts_dict[the_other_parameter_set]:

            processor_counts = self.timing_data["total"][the_other_parameter_set][mode].processor_counts
            if 1 in processor_counts:
                the_mode = mode
                the_serial_index = processor_counts.index(1)
                break

        # Only if a serial run was found for any mode
        if the_mode is None: return

        # Scale the runtimes of the emission phases
        stellar_time = self.timing_data["stellar"][parameter_set][the_mode].times[the_serial_index] * npackages_factor
        dust_time = self.timing_data["dust"][parameter_set][the_mode].times[the_serial_index] * npackages_factor

        # Other runtimes
        setup_time = self.timing_data["setup"][parameter_set][the_mode].times[the_serial_index]
        spectra_time = self.timing_data["spectra"][parameter_set][the_mode].times[the_serial_index]
        writing_time = self.timing_data["writing"][parameter_set][the_mode].times[the_serial_index]
        waiting_time = self.timing_data["waiting"][parameter_set][the_mode].times[the_serial_index]
        communication_time = 0.0  # serial run
        intermediate_time = self.timing_data["intermediate"][parameter_set][the_mode].times[the_serial_index]

        # Calculate other
        densities_comm_time = 0.0  # serial run
        stellarabsorption_comm_time = 0.0  # serial run
        dustabsorption_comm_time = 0.0  # serial run
        emission_comm_time = 0.0  # serial run
        instruments_comm_time = 0.0  # serial run
        total_time = setup_time + stellar_time + spectra_time + dust_time + writing_time + waiting_time + communication_time + intermediate_time

        # total
        self.serial_timing[parameter_set]["total"].time = total_time
        self.serial_timing[parameter_set]["total"].error = 0.0
        self.serial_timing_ncores[parameter_set]["total"] = 1

        # setup
        self.serial_timing[parameter_set]["setup"].time = setup_time
        self.serial_timing[parameter_set]["setup"].error = 0.0
        self.serial_timing_ncores[parameter_set]["setup"] = 1

        # stellar
        self.serial_timing[parameter_set]["stellar"].time = stellar_time
        self.serial_timing[parameter_set]["stellar"].error = 0.0
        self.serial_timing_ncores[parameter_set]["stellar"] = 1

        # spectra
        self.serial_timing[parameter_set]["spectra"].time = spectra_time
        self.serial_timing[parameter_set]["spectra"].error = 0.0
        self.serial_timing_ncores[parameter_set]["spectra"] = 1

        # dust
        self.serial_timing[parameter_set]["dust"].time = dust_time
        self.serial_timing[parameter_set]["dust"].error = 0.0
        self.serial_timing_ncores[parameter_set]["dust"] = 1

        # writing
        self.serial_timing[parameter_set]["writing"].time = writing_time
        self.serial_timing[parameter_set]["writing"].error = 0.0
        self.serial_timing_ncores[parameter_set]["writing"] = 1

        # waiting
        self.serial_timing[parameter_set]["waiting"].time = waiting_time
        self.serial_timing[parameter_set]["waiting"].error = 0.0
        self.serial_timing_ncores[parameter_set]["waiting"] = 1

        # communication
        self.serial_timing[parameter_set]["communication"].time = communication_time
        self.serial_timing[parameter_set]["communication"].error = 0.0
        self.serial_timing_ncores[parameter_set]["communication"] = 1

        # dust densities communication
        self.serial_timing[parameter_set]["dust densities communication"].time = densities_comm_time
        self.serial_timing[parameter_set]["dust densities communication"].error = 0.0
        self.serial_timing_ncores[parameter_set]["dust densities communication"] = 1

        # stellar absorption communication
        self.serial_timing[parameter_set]["stellar absorption communication"].time = stellarabsorption_comm_time
        self.serial_timing[parameter_set]["stellar absorption communication"].error = 0.0
        self.serial_timing_ncores[parameter_set]["stellar absorption communication"] = 1

        # dust absorption communication
        self.serial_timing[parameter_set]["dust absorption communication"].time = dustabsorption_comm_time
        self.serial_timing[parameter_set]["dust absorption communication"].error = 0.0
        self.serial_timing_ncores[parameter_set]["dust absorption communication"] = 1

        # emission spectra communication
        self.serial_timing[parameter_set]["emission spectra communication"].time = emission_comm_time
        self.serial_timing[parameter_set]["emission spectra communication"].error = 0.0
        self.serial_timing_ncores[parameter_set]["emission spectra communication"] = 1

        # instruments communication
        self.serial_timing[parameter_set]["instruments communication"].time = instruments_comm_time
        self.serial_timing[parameter_set]["instruments communication"].error = 0.0
        self.serial_timing_ncores[parameter_set]["instruments communication"] = 1

        # intermediate
        self.serial_timing[parameter_set]["intermediate"].time = intermediate_time
        self.serial_timing[parameter_set]["intermediate"].error = 0.0
        self.serial_timing_ncores[parameter_set]["intermediate"] = 1

    # -----------------------------------------------------------------

    def set_missing_serial_timing_false(self, parameter_set):

        """
        This function ...
        :param parameter_set:
        :return:
        """

        # Give warning
        log.warning("Serial (one core) timing data not found, searching for longest runtime with the least amount of cores for each phase (any parallelization mode)")

        # Loop over the phases
        for phase in self.timing_data[parameter_set]:

            max_time = 0.0
            max_time_error = None
            max_time_ncores = float("+inf")

            # Loop over the modes
            for mode in self.timing_data[parameter_set][phase]:

                index = np.argmax(self.timing_data[parameter_set][phase][mode].times)
                max_time_mode = self.timing_data[parameter_set][phase][mode].times[index]
                max_time_error_mode = self.timing_data[parameter_set][phase][mode].errors[index]
                ncores = self.timing_data[parameter_set][phase][mode].processor_counts[index]

                # Check condition and adapt accordingly
                if ncores < max_time_ncores and max_time_mode > max_time:

                    max_time = max_time_mode
                    max_time_error = max_time_error_mode
                    max_time_ncores = ncores

            # Set the time and error
            self.serial_timing[parameter_set][phase].time = max_time
            self.serial_timing[parameter_set][phase].error = max_time_error
            self.serial_timing_ncores[parameter_set][phase] = max_time_ncores

    # -----------------------------------------------------------------

    def set_missing_serial_memory(self, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        # Check if serial memory data is found for the different parameter sets
        for parameter_set in parameters_modes_processor_counts_dict:

            # Check
            if len(self.serial_memory[parameter_set]) == 0:

                # Give warning
                log.warning("Serial (one core) memory data not found, searching for highest memory usage with the least amount of cores for each phase (any parallelization mode)")

                # Loop over the phases
                for phase in self.memory_data[parameter_set]:

                    max_memory = 0.0
                    max_memory_error = None
                    max_memory_ncores = float("+inf")

                    # Loop over the modes
                    for mode in self.memory_data[parameter_set][phase]:

                        index = np.argmax(self.memory_data[parameter_set][phase][mode].memory)
                        max_memory_mode = self.memory_data[parameter_set][phase][mode].memory[index]
                        max_memory_error_mode = self.memory_data[parameter_set][phase][mode].errors[index]
                        ncores = self.memory_data[parameter_set][phase][mode].processor_counts[index]

                        # Check condition and adapt accordingly
                        if ncores < max_memory_ncores and max_memory_mode > max_memory:
                            max_memory = max_memory_mode
                            max_memory_error = max_memory_error_mode
                            max_memory_ncores = ncores

                    # Set the memory usage and error
                    self.serial_memory[parameter_set][phase].memory = max_memory
                    self.serial_memory[parameter_set][phase].error = max_memory_error
                    self.serial_memory_ncores[parameter_set][phase] = max_memory_ncores

    # -----------------------------------------------------------------

    @property
    def parameter_sets(self):

        """
        This function ...
        :return:
        """

        return self.timing_data["total"].keys()

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting ...")

        # Fit timing
        if self.needs_timing: self.fit_timing()

        # Fit memory
        if self.needs_memory: self.fit_memory()

    # -----------------------------------------------------------------

    def fit_timing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting the timing data ...")

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets:

            # Loop over the phases
            for phase in self.config.phases:

                # Skip phases for which a serial timing is not present
                #if not self.has_serial_timing(phase):  # continue
                #    # Give warning
                #    log.warning("Serial (one core) timing data not found, using longest runtime (any parallelization mode) for normalizing the speedups (and efficiencies) for the fit")

                # Debugging
                log.debug("Fitting timing data for " + phase_names[phase] + " ...")

                # Get the serial runtime (and error) for this phase (create a Quantity object)
                #serial_time = self.serial_timing[phase].time
                #serial_error = self.serial_timing[phase].error
                #serial = Quantity(serial_time, serial_error)

                # Create a dictionary that stores the fitted parameters for each different mode
                parameters = dict()
                parameter_errors = dict()

                # Generate the fit function
                behaviour = generate_scaling_behaviour(phase)
                fit_function, nparameters = generate_fit_function(behaviour)

                # Debugging
                log.debug("Expected scaling behaviour: " + str(behaviour))

                # Loop over the different parallelization modes (the different curves)
                for mode in self.timing_data[parameter_set][phase]:

                    # Get the list of processor counts, runtimes and errors
                    processor_counts = self.timing_data[parameter_set][phase][mode].processor_counts
                    times = self.timing_data[parameter_set][phase][mode].times
                    errors = self.timing_data[parameter_set][phase][mode].errors

                    # Sort the lists
                    processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                    # Set the weights of the different timing points for the fitting procedure
                    weigths = errors if not np.any(np.isinf(errors)) else None
                    if np.count_nonzero(errors) == 0: weights = None

                    # Calculate the normalized processor counts (relative to the number of processors used for the serial run)
                    normalized_processor_counts = processor_counts / self.serial_timing_ncores[parameter_set][phase]

                    from scipy.optimize.minpack import _initialize_feasible, prepare_bounds
                    n = nparameters
                    bounds = (-np.inf, np.inf)
                    lb, ub = prepare_bounds(bounds, n)
                    p0 = _initialize_feasible(lb, ub)

                    # Fit parameters for the speedups to Amdahl's law
                    popt, pcov = curve_fit(fit_function, normalized_processor_counts, times, sigma=weigths, absolute_sigma=False, p0=p0)
                    perr = np.sqrt(np.diag(pcov))

                    parameters[mode] = popt
                    parameter_errors[mode] = perr

                # If output path is specified, write parameter files
                if self.config.output is not None:

                    mode_list = []
                    other_columns = defaultdict(list)

                    # Fill columns
                    for mode in parameters:
                        mode_list.append(mode)
                        for i in range(len(parameters[mode])):
                            other_columns[alphabet[i]].append(parameters[mode][i])
                            other_columns[alphabet[i] + "_error"].append(parameter_errors[mode][i])

                    # Create a data file to contain the fitted parameters
                    directory = self.config.output
                    parameter_file_path = fs.join(directory, "parameters_timing_" + phase + ".dat")

                    # Create the parameters table and write to file
                    data = [mode_list]
                    names = ["Parallelization mode"]

                    for column_name in other_columns:
                        data.append(other_columns[column_name])
                        names.append(column_name)

                    table = Table(data=data, names=names)
                    tables.write(table, parameter_file_path)

                # Add the parameters
                self.timing_fit_functions[parameter_set][phase] = fit_function
                self.timing_fit_parameters[parameter_set][phase] = parameters
                self.timing_fit_parameter_errors[parameter_set][phase] = parameter_errors

    # -----------------------------------------------------------------

    def fit_timing_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting the timing data ...")

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets:

            # Loop over the phases
            for phase in self.config.phases:

                # Skip certain phases
                if phase not in parallel_phases: continue

                # Skip phases for which a serial timing is not present
                if not self.has_serial_timing(phase): #continue

                    # Give warning
                    log.warning("Serial (one core) timing data not found, using longest runtime (any parallelization mode) for normalizing the speedups (and efficiencies) for the fit")

                # Get the serial runtime (and error) for this phase (create a Quantity object)
                serial_time = self.serial_timing[phase].time
                serial_error = self.serial_timing[phase].error
                serial = Quantity(serial_time, serial_error)

                # Create a dictionary that stores the fitted parameters for each different mode
                parameters = dict()

                # Loop over the different parallelization modes (the different curves)
                for mode in self.timing_data[phase]:

                    # Get the list of processor counts, runtimes and errors
                    processor_counts = self.timing_data[phase][mode].processor_counts
                    times = self.timing_data[phase][mode].times
                    errors = self.timing_data[phase][mode].errors

                    # Sort the lists
                    processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

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

                    # Set the weights of the different speedup points for the fitting procedure
                    speedup_weigths = speedup_errors if not np.any(np.isinf(speedup_errors)) else None
                    if np.count_nonzero(speedup_weigths) == 0: speedup_weigths = None

                    # Calculate the normalized processor counts (relative to the number of processors used for the serial run)
                    normalized_processor_counts = processor_counts / self.serial_timing_ncores[phase]

                    # Fit (standard or modified) Amdahl's law to the speedups
                    if len(processor_counts) < 10:

                        # Fit parameters for the speedups to Amdahl's law
                        popt, pcov = curve_fit(amdahl_law, normalized_processor_counts, speedups, sigma=speedup_weigths, absolute_sigma=False)
                        perr = np.sqrt(np.diag(pcov))
                        parameters[mode] = Map({"p": popt[0], "p_error": perr[0], "a": 0.0, "a_error": 0.0, "b": 0.0, "b_error": 0.0, "c": 0.0, "c_error": 0.0})

                    else:

                        # Fit parameters for the speedups to Amdahl's law
                        popt, pcov = curve_fit(modified_amdahl_law, normalized_processor_counts, speedups, sigma=speedup_weigths, absolute_sigma=False)
                        perr = np.sqrt(np.diag(pcov))
                        parameters[mode] = Map({"p": popt[0], "p_error": perr[0], "a": popt[1], "a_error": perr[1], "b": popt[2], "b_error": perr[2], "c": popt[3], "c_error": perr[3]})

                # If output path is specified, write parameter files
                if self.config.output is not None:

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

                    # Create a data file to contain the fitted parameters
                    directory = self.config.output
                    parameter_file_path = fs.join(directory, "parameters_timing_" + phase + ".dat")

                    # Create the parameters table and write to file
                    data = [mode_list, p_list, p_error_list, a_list, a_error_list, b_list, b_error_list, c_list,
                            c_error_list]
                    names = ["Parallelization mode", "Parallel fraction p", "Error on p", "Parameter a",
                             "Error on a", "Parameter b", "Error on b", "Parameter c", "Error on c"]
                    table = Table(data=data, names=names)
                    table.write(parameter_file_path, format="ascii.commented_header")

                # Add the parameters
                self.timing_fit_parameters[phase] = parameters

            ## COMMUNICATION

            # Get the serial runtime (and error) for this phase (create a Quantity object)
            #serial_time = self.serial_timing["communication"].time
            #serial_error = self.serial_timing["communication"].error
            #serial = Quantity(serial_time, serial_error)

            # Create a dictionary that stores the fitted parameters for each different mode
            parameters = dict()

            # Loop over the different parallelization modes (the different curves)
            for mode in self.timing_data["communication"]:

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.timing_data[phase][mode].processor_counts
                times = self.timing_data[phase][mode].times
                errors = self.timing_data[phase][mode].errors

                # Sort the lists
                processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                # Set the weights of the different runtime points for the fitting procedure
                weights = errors if not np.any(np.isinf(errors)) else None
                if np.count_nonzero(weights) == 0: weights = None

                # Calculate the normalized processor counts (relative to the number of processors used for the serial run)
                normalized_processor_counts = processor_counts / self.serial_timing_ncores["communication"]

                # Fit logarithmic + linear curve to the data
                popt, pcov = curve_fit(communication_time_scaling, normalized_processor_counts, times, sigma=weights, absolute_sigma=False)
                perr = np.sqrt(np.diag(pcov))
                parameters[mode] = Map({"a": popt[0], "a_error": perr[0], "b": popt[1], "b_error": perr[1], "c": popt[2], "c_error": perr[2]})

            # If output path is specified, write parameter files
            if self.config.output is not None:

                mode_list = []
                a_list = []
                a_error_list = []
                b_list = []
                b_error_list = []
                c_list = []
                c_error_list = []

                for mode in parameters:

                    mode_list.append(mode)
                    a_list.append(parameters[mode].a)
                    a_error_list.append(parameters[mode].a_error)
                    b_list.append(parameters[mode].b)
                    b_error_list.append(parameters[mode].b_error)
                    c_list.append(parameters[mode].c)
                    c_error_list.append(parameters[mode].c_error)

                # Create a data file to contain the fitted parameters
                directory = self.config.output
                parameter_file_path = fs.join(directory, "parameters_timing_communication.dat")

                # Create the parameters table and write to file
                data = [mode_list, a_list, a_error_list, b_list, b_error_list, c_list, c_error_list]
                names = ["Parallelization mode", "Parameter a", "Error on a", "Parameter b", "Error on b", "Parameter c", "Error on c"]
                table = Table(data=data, names=names)
                table.write(parameter_file_path, format="ascii.commented_header")

            # Add the parameters
            self.timing_fit_parameters["communication"] = parameters

    # -----------------------------------------------------------------

    def fit_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting the memory data ...")

        # Loop over the phases
        for phase in self.config.phases:

            # Get the serial (1 process) memory consumption (and error) for this phase (create a Quantity object)
            #serial_memory = self.serial_memory[phase].memory
            #serial_error = self.serial_memory[phase].error
            #serial = Quantity(serial_memory, serial_error)

            # Create a dictionary that stores the fitted parameters for each different mode
            parameters = dict()

            # Loop over the different parallelization modes (the different curves)
            for mode in self.memory_data[phase]:

                # Skip modes that are not data parallel
                if not "task+data" in mode: continue

                # Get the list of processor counts, memory and errors
                processor_counts = self.memory_data[phase][mode].processor_counts
                memories = self.memory_data[phase][mode].memory
                errors = self.memory_data[phase][mode].errors

                # Sort the lists
                processor_counts, memories, errors = sort_lists(processor_counts, memories, errors, to_arrays=True)

                # Set the weights of the different memory points for the fitting procedure
                weights = errors if not np.any(np.isinf(errors)) else None
                if np.count_nonzero(weights) == 0: weights = None

                # Check data
                if np.any(np.isinf(memories)) or np.any(np.isnan(memories)):
                    log.warning("Fitting not possible to the memory data of " + mode + " for " + phase_names[phase].lower() + " (nans and/or infs)")
                    continue

                # Get list of nprocesses
                nprocesses = nprocesses_from_mode(mode, processor_counts)

                # Fit (standard or modified) memory scaling law
                if len(processor_counts) < 5:

                    popt, pcov = curve_fit(memory_scaling, nprocesses, memories, sigma=weights, absolute_sigma=False)
                    perr = np.sqrt(np.diag(pcov))
                    parameters[mode] = Map({"a": popt[0], "a_error": perr[0], "b": popt[1], "b_error": perr[1], "c": 0.0, "c_error": 0.0})

                else:

                    popt, pcov = curve_fit(modified_memory_scaling, nprocesses, memories, sigma=weights, absolute_sigma=False)
                    perr = np.sqrt(np.diag(pcov))
                    parameters[mode] = Map({"a": popt[0], "a_error": perr[0], "b": popt[1], "b_error": perr[1], "c": popt[2], "c_error": perr[2]})

            # If output path is specified, write parameter files
            if self.config.output is not None:

                mode_list = []
                a_list = []
                a_error_list = []
                b_list = []
                b_error_list = []
                c_list = []
                c_error_list = []

                for mode in parameters:

                    mode_list.append(mode)
                    a_list.append(parameters[mode].a)
                    a_error_list.append(parameters[mode].a_error)
                    b_list.append(parameters[mode].b)
                    b_error_list.append(parameters[mode].b_error)
                    c_list.append(parameters[mode].c)
                    c_error_list.append(parameters[mode].c_error)

                # Create a data file to contain the fitted parameters
                directory = self.config.output
                parameter_file_path = fs.join(directory, "parameters_memory_" + phase + ".dat")

                # Create the parameters table and write to file
                data = [mode_list, a_list, a_error_list, b_list, b_error_list, c_list, c_error_list]
                names = ["Parallelization mode", "Parameter a", "Error on a", "Parameter b", "Error on b", "Parameter c", "Error on c"]
                table = Table(data=data, names=names)
                table.write(parameter_file_path, format="ascii.commented_header")

            # Add the parameters
            self.memory_fit_parameters[phase] = parameters

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
        if "CPU-time" in self.config.properties: self.plot_cpu_times()

        # Memory
        if "memory" in self.config.properties: self.plot_memory()

        # Memory gain
        if "memory-gain" in self.config.properties: self.plot_memory_gain()

        # Total memory
        if "total-memory" in self.config.properties: self.plot_total_memory()

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

    def plot_runtimes_phase(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes for the " + phase_names[phase] + " ...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.config.figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.timing_data:

            # Loop over the different parallelization modes (the different curves)
            for mode in self.timing_data[parameter_set][phase]:

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.timing_data[parameter_set][phase][mode].processor_counts
                times = self.timing_data[parameter_set][phase][mode].times
                errors = self.timing_data[parameter_set][phase][mode].errors

                # Sort the lists
                processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                # Determine the label
                if parameter_set == "empty": label = mode
                else: label = mode + " [" + ", ".join([name + ": " + str(parameter_set[name]) for name in parameter_set]) + "]"

                # Plot the data points for this mode
                plt.errorbar(processor_counts, times, errors, marker='.', label=label)

                # Add the appropriate ticks
                ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Plot the fit
        if not self.config.hybridisation and self.config.fit and self.config.plot_fit:

            # Loop over the parameter sets
            for parameter_set in self.timing_fit_functions:

                # Get the fit function
                fit_function = self.timing_fit_functions[parameter_set][phase]

                # Get the number of processers taken as the reference for normalization, and thus calculation of the speedups and as reference for the fit
                #reference_ncores = self.serial_timing_ncores[phase]

                fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                for mode in self.timing_fit_parameters[parameter_set][phase]:

                    # Get the parameter values
                    parameters = self.timing_fit_parameters[parameter_set][phase][mode]

                    # Get the parameter errors
                    parameter_errors = self.timing_fit_parameter_errors[parameter_set][phase][mode]

                    # Calculate the fitted times
                    fit_times = [fit_function(n, *parameters) for n in fit_ncores]

                    # Add the plot
                    plt.plot(fit_ncores, fit_times, color="grey")

        # Plot curve of communication times
        #if not self.config.hybridisation and self.config.fit and self.config.plot_fit and phase == "communication":

            # Get the fit parameters
            #parameters = self.timing_fit_parameters[phase]

            # Get the number of processers taken as the reference for normalization, and thus calculation of the speedups and as reference for the fit
            #reference_ncores = self.serial_timing_ncores[phase]

            # Plot the fitted speedup curves and write the parameters to the file
            #fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
            #for mode in parameters:

                # Get the parameter values
                #a = parameters[mode].a
                #b = parameters[mode].b
                #c = parameters[mode].c

                # Calculate the fitted times
                #fit_times = [communication_time_scaling(n / float(reference_ncores), a, b, c) for n in fit_ncores]

                # Add the plot
                #plt.plot(fit_ncores, fit_times, color="grey")

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
        if self.config.hybridisation: plt.legend(title="Number of cores")
        else: plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Scaling of the " + phase_labels[phase].lower())

        # Set file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "runtimes_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_speedups_phase(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        # Inform the user of the fact that the speedups are being calculated and plotted
        log.info("Plotting the speedups for the " + phase_names[phase] + " ...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.config.figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Get the serial runtime (and error) for this phase (create a Quantity object)
        serial_time = self.serial_timing[phase].time
        serial_error = self.serial_timing[phase].error
        serial = Quantity(serial_time, serial_error)

        # Keep track of the minimal and maximal speedup
        speedup_range = RealRange.zero()

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.timing_data:

            # Loop over the different parallelization modes (the different curves)
            for mode in self.timing_data[parameter_set][phase]:

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.timing_data[parameter_set][phase][mode].processor_counts
                times = self.timing_data[parameter_set][phase][mode].times
                errors = self.timing_data[parameter_set][phase][mode].errors

                # Sort the lists
                processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                # Calculate the speedups and the errors on the speedups
                speedups = []
                speedup_errors = []
                for i in range(len(processor_counts)):

                    # Create a quantity for the current runtime
                    time = Quantity(times[i], errors[i])

                    # Calculate the speedup based on the current runtime and the serial runtime
                    speedup = serial / time

                    # Adjust the speedup range
                    speedup_range.adjust(speedup.value)

                    # Add the value and the propagated error of the speedup to the appropriate lists
                    speedups.append(speedup.value)
                    speedup_errors.append(speedup.error)

                # Determine the label
                if parameter_set == "empty": label = mode
                else: label = mode + " [" + ", ".join([name + ": " + str(parameter_set[name]) for name in parameter_set]) + "]"

                # Plot the data points for this curve
                plt.errorbar(processor_counts, speedups, speedup_errors, marker='.', label=label)

                # Add the appropriate ticks
                ticks |= set(processor_counts)

        # Use a logarithmic scale for both axes
        plt.xscale('log')
        plt.yscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Plot the fit
        if not self.config.hybridisation and self.config.fit and self.config.plot_fit:

            # Loop over the parameter sets
            for parameter_set in self.timing_fit_functions:

                # Get the fit function
                fit_function = self.timing_fit_functions[parameter_set][phase]

                fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                for mode in self.timing_fit_parameters[parameter_set][phase]:

                    # Get the parameter values
                    parameters = self.timing_fit_parameters[parameter_set][phase][mode]

                    # Get the parameter errors
                    parameter_errors = self.timing_fit_parameter_errors[parameter_set][phase][mode]

                    # Calculate the fitted speedups
                    fit_speedups = [serial.value / fit_function(n, *parameters) for n in fit_ncores]

                    # Add the plot
                    plt.plot(fit_ncores, fit_speedups, color="grey")

        #if not self.config.hybridisation and self.config.fit and self.config.plot_fit:

            # Get the fit parameters
            #parameters = self.timing_fit_parameters[phase]

            # Get the number of processers taken as the reference for normalization, and thus calculation of the speedups and as reference for the fit
            #reference_ncores = self.serial_timing_ncores[phase]

            # Plot the fitted speedup curves and write the parameters to the file
            #fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
            #for mode in parameters:

                # Get the parameter values
                #p = parameters[mode].p
                #a = parameters[mode].a
                #b = parameters[mode].b
                #c = parameters[mode].c

                # Calculate the fitted speedups
                #fit_speedups = [modified_amdahl_law(n/reference_ncores, p, a, b, c) for n in fit_ncores]

                # Add the plot
                #plt.plot(fit_ncores, fit_speedups, color="grey")

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        if not self.config.hybridisation: ax.set_yticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        #if not self.config.hybridisation: plt.ylim(ticks[0], ticks[-1])
        plt.ylim(speedup_range.min, speedup_range.max)
        plt.grid(True)

        # Plot a line that denotes linear scaling (speedup = nthreads)
        if not self.config.hybridisation: plt.plot(ticks, ticks, linestyle='--')

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize='large')
        plt.ylabel(phase_labels[phase] + " speedup $S$", fontsize='large')
        if self.config.hybridisation: plt.legend(title="Number of cores")
        else: plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Speedup of the " + phase_labels[phase].lower())

        # Set file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "speedups_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_efficiencies_phase(self, phase):

        """
        This function creates a PDF plot showing the efficiency as a function of the number of threads.
        Efficiency is defined as T(1)/T(N)/N. It is a dimensionless quantity <= 1.
        The function takes the following (optional) arguments:
        :param phase:
        :return:
        """

        # Inform the user of the fact that the efficiencies are being calculated and plotted
        log.info("Calculating and plotting the efficiencies for the " + phase_names[phase] + " ...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.config.figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Get the serial runtime (and error) for this phase (create a Quantity object)
        serial_time = self.serial_timing[phase].time
        serial_error = self.serial_timing[phase].error
        serial = Quantity(serial_time, serial_error)
        serial_ncores = self.serial_timing_ncores[phase] if phase in self.serial_timing_ncores else 1

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.timing_data:

            # Loop over the different parallelization modes (the different curves)
            for mode in self.timing_data[parameter_set][phase]:

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.timing_data[parameter_set][phase][mode].processor_counts
                times = self.timing_data[parameter_set][phase][mode].times
                errors = self.timing_data[parameter_set][phase][mode].errors

                # Sort the lists
                processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

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
                    efficiency = speedup.value / (ncores[i]/serial_ncores)
                    efficiency_error = speedup.error / (ncores[i]/serial_ncores)

                    # Add the value and the propagated error of the efficiency to the appropriate lists
                    efficiencies.append(efficiency)
                    efficiency_errors.append(efficiency_error)

                # Determine the label
                if parameter_set == "empty": label = mode
                else: label = mode + " [" + ", ".join([name + ": " + str(parameter_set[name]) for name in parameter_set]) + "]"

                # Plot the data points for this curve
                plt.errorbar(processor_counts, efficiencies, efficiency_errors, marker='.', label=label)

                # Add the appropriate ticks
                ticks |= set(processor_counts)

        # Use a logaritmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Plot fit
        if not self.config.hybridisation and self.config.fit and self.config.plot_fit:

            # Loop over the parameter sets
            for parameter_set in self.timing_fit_functions:

                # Get the fit function
                fit_function = self.timing_fit_functions[parameter_set][phase]

                fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                for mode in self.timing_fit_parameters[parameter_set][phase]:

                    # Get the parameter values
                    parameters = self.timing_fit_parameters[parameter_set][phase][mode]

                    # Get the parameter errors
                    parameter_errors = self.timing_fit_parameter_errors[parameter_set][phase][mode]

                    # Calculate the fitted efficiencies
                    fit_efficiencies = [serial.value / fit_function(n, *parameters) / n for n in fit_ncores]

                    # Add the plot
                    plt.plot(fit_ncores, fit_efficiencies, color="grey")

        # Plot fit
        #if not self.config.hybridisation and self.config.fit and self.config.plot_fit and self.has_timing_fit(phase):

            # Get the number of processers taken as the reference for normalization, and thus calculation of the speedups and as reference for the fit
            #reference_ncores = self.serial_timing_ncores[phase]

            # Plot the fitted speedup curves
            #fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
            #for mode in self.timing_fit_parameters[phase]:

                # Get the parameter values
                #p = self.timing_fit_parameters[phase][mode].p
                #a = self.timing_fit_parameters[phase][mode].a
                #b = self.timing_fit_parameters[phase][mode].b
                #c = self.timing_fit_parameters[phase][mode].c

                # Calculate the fitted efficiencies
                #fit_efficiencies = [modified_amdahl_law(n/float(reference_ncores), p, a, b, c) / (n*reference_ncores) for n in fit_ncores]

                # Add the plot
                #plt.plot(fit_ncores, fit_efficiencies, color="grey")

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        #plt.ylim(0, 1.1)
        plt.grid(True)

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes $N_p$", fontsize='large')
        else: plt.xlabel("Number of cores $N_c$", fontsize="large")
        plt.ylabel(phase_labels[phase] + " efficiency $\epsilon$", fontsize='large')
        if self.config.hybridisation: plt.legend(title="Number of cores")
        else: plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Efficiency of the " + phase_labels[phase].lower())

        # Determine the path
        if self.config.output is not None: file_path = fs.join(self.config.output, "efficiencies_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_cpu_times_phase(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        # Inform the user
        log.info("Plotting the CPU times for the " + phase_names[phase] + " ...")

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.config.figsize)
        plt.clf()

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.timing_data:

            # Loop over the different parallelization modes (the different curves)
            for mode in self.timing_data[parameter_set][phase]:

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.timing_data[parameter_set][phase][mode].processor_counts
                times = self.timing_data[parameter_set][phase][mode].times
                errors = self.timing_data[parameter_set][phase][mode].errors

                # Sort the lists
                processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                # Get array of number of used cores
                if self.config.hybridisation: ncores = np.ones(len(processor_counts)) * int(mode.split(" cores")[0])
                else: ncores = processor_counts

                # Get list of process count
                #processes = nprocesses_from_mode(mode, processor_counts)

                # Multiply to get total
                times *= ncores
                errors *= ncores

                # Determine the label
                if parameter_set == "empty": label = mode
                else: label = mode + " [" + ", ".join([name + ": " + str(parameter_set[name]) for name in parameter_set]) + "]"

                # Plot the data points for this mode
                plt.errorbar(processor_counts, times, errors, marker='.', label=label)

                # Add the appropriate ticks
                ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale("log")
        plt.yscale("log")

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1] * 2)

        # Plot curve of communication times
        #if not self.config.hybridisation and self.config.fit and self.config.plot_fit and phase == "communication":

            # Get the fit parameters
            #parameters = self.timing_fit_parameters[phase]

            # Get the number of processers taken as the reference for normalization, and thus calculation of the speedups and as reference for the fit
            #reference_ncores = self.serial_timing_ncores[phase]

            # Plot the fitted speedup curves and write the parameters to the file
            #fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
            #for mode in parameters:

                # Get the parameter values
                #a = parameters[mode].a
                #b = parameters[mode].b
                #c = parameters[mode].c

                # Calculate the fitted times
                #fit_times = [communication_time_scaling(n/reference_ncores, a, b, c) * n for n in fit_ncores]

                # Add the plot
                #plt.plot(fit_ncores, fit_times, color="grey")

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
        if self.config.hybridisation: plt.legend(title="Number of cores")
        else: plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Scaling of the total CPU time of the " + phase_labels[phase].lower())

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "cpu_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_memory_phase(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        # Inform the user
        log.info("Plotting the memory scaling for the " + phase_names[phase] + " ...")

        if len(self.memory_data[phase]) == 0:
            log.warning("No memory data for the " + phase_names[phase])
            return

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.config.figsize)
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
            processor_counts, memories, errors = sort_lists(processor_counts, memories, errors, to_arrays=True)

            # Plot the data points for this mode
            plt.errorbar(processor_counts, memories, errors, marker='.', label=mode)

            # Add the appropriate ticks
            ticks |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        #ticks.append(ticks[-1] * 2)

        # Plot fit
        if not self.config.hybridisation and self.config.fit and self.config.plot_fit:

            # Plot the fitted curves
            fit_ncores = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
            for mode in self.memory_fit_parameters[phase]:

                # Get the parameter values
                a = self.memory_fit_parameters[phase][mode].a
                b = self.memory_fit_parameters[phase][mode].b
                c = self.memory_fit_parameters[phase][mode].c

                # Calculate the fitted memory usages
                fit_memories = [modified_memory_scaling(nprocesses_from_mode_single(mode, ncores), a, b, c) for ncores in fit_ncores]

                # Add the plot
                plt.plot(fit_ncores, fit_memories, color="grey")

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
        if self.config.hybridisation: plt.legend(title="Number of cores")
        else: plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Scaling of the memory usage (per process) of the " + phase + " phase")

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "memory_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_memory_gain_phase(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        # Inform the user
        log.info("Plotting the memory gain scaling for the " + phase_names[phase] + " ...")

        if len(self.memory_data[phase]) == 0:
            log.warning("No memory data for the " + phase_names[phase])
            return

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.config.figsize)
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
            processor_counts, memories, errors = sort_lists(processor_counts, memories, errors, to_arrays=True)

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
        if self.config.hybridisation: plt.legend(title="Number of cores")
        else: plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Scaling of the memory gain (serial memory usage per process / memory usage per process) of the " + phase_labels[phase].lower())

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "memorygain_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path)
        else: plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_total_memory_phase(self, phase):

        """
        This function ...
        :param phase:
        :return:
        """

        # Inform the user
        log.info("Plotting the total memory scaling (all processes combined) for the " + phase_names[phase] + " ...")

        if len(self.memory_data[phase]) == 0:
            log.warning("No memory data for the " + phase_names[phase])
            return

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.config.figsize)
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
            processor_counts, memories, errors = sort_lists(processor_counts, memories, errors, to_arrays=True)

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
        if self.config.hybridisation: plt.legend(title="Number of cores")
        else: plt.legend(title="Parallelization modes")

        # Set the plot title
        plt.title("Memory scaling for " + phase_labels[phase])

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "totalmemory_" + phase + ".pdf")
        else: file_path = None

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

            if self.config.output is not None:
                plot_file_path = fs.join(self.config.output, "timeline_" + mode + ".pdf")
            else: plot_file_path = None

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the timing table
        self.write_timing()

        # Write the memory table
        self.write_memory()

    # -----------------------------------------------------------------

    def write_timing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing timing table ...")

        table_path = fs.join(self.config.output, "timing.dat")
        self.timing.saveto(table_path)

    # -----------------------------------------------------------------

    def write_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing memory table ...")

        table_path = fs.join(self.config.output, "memory.dat")
        self.memory.saveto(table_path)

# -----------------------------------------------------------------

def ncores_from_mode_single(mode, nprocesses):

    """
    This function ...
    :param mode:
    :param nprocesses:
    :return:
    """

    if "[" in mode: mode = mode.split(" [")[0]

    # Get number of cores
    if mode == "multithreading": raise RuntimeError("Number of cores cannot be determined from nprocesses = " + str(nprocesses) + " in " + mode + " mode")
    elif "multiprocessing" in mode: ncores = nprocesses
    elif "hybrid" in mode:
        threads = int(mode.split("(")[1].split(" threads)")[0])
        ncores = nprocesses * threads
    else: raise ValueError("Invalid mode: " + mode)

    # Return the number of cores
    return ncores

# -----------------------------------------------------------------

def nprocesses_from_mode_single(mode, nprocessors):

    """
    This function ...
    :param mode:
    :param nprocessors:
    :return:
    """

    if "[" in mode: mode = mode.split(" [")[0]

    # Get number of processes
    if mode == "multithreading": processes = 1
    elif "multiprocessing" in mode: processes = nprocessors
    elif "hybrid" in mode:
        threads = int(mode.split("(")[1].split(" threads)")[0])
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

    if "[" in mode: mode = mode.split(" [")[0]

    # Get number of processes
    if mode == "multithreading": processes = np.ones(len(processor_counts))
    elif "multiprocessing" in mode: processes = processor_counts
    elif "hybrid" in mode:
        threads = int(mode.split("(")[1].split(" threads)")[0])
        processes = processor_counts / threads
    else: raise ValueError("Invalid mode: " + mode)

    # Return the list of proces count
    return processes

# -----------------------------------------------------------------

def amdahl_law(n, p):

    """
    This function defines Amdahl's law for the speedup
    :param n:
    :param p:
    :return:
    """

    return 1.0 / (1 - p + p / n)

# -----------------------------------------------------------------

def modified_amdahl_law(n, p, a, b, c):

    """
    This function defines a modified version of Amdahl's law, which accounts for different kinds of overhead
    :param n:
    :param p:
    :param a:
    :param b:
    :param c:
    :return:
    """

    return 1.0 / (1 - p + p / n + a + b * n + c * n**2)

# -----------------------------------------------------------------

def memory_scaling(n, a, b):

    """
    This function ...
    :param n: number of processes
    :param a:
    :param b:
    :return:
    """

    return a / n + b

# -----------------------------------------------------------------

def modified_memory_scaling(n, a, b, c):

    """
    This function ...
    :param n: number of processes
    :param a:
    :param b:
    :param c:
    :return:
    """

    return a / n + b + c * n

# -----------------------------------------------------------------

def serial_time_scaling(n, a):

    """
    This function ...
    :param n:
    :param a:
    :return:
    """

    return a

# -----------------------------------------------------------------

def parallel_time_scaling(n, a):

    """
    This function ...
    :param n:
    :param a:
    :return:
    """

    return a/n

# -----------------------------------------------------------------

def waiting_time_scaling(n, ):

    """
    This function ...
    :param n:
    :return:
    """

    return a * n

# -----------------------------------------------------------------

def communication_time_scaling(n, a, b, c):

    """
    This function ...
    :param n:
    :param a:
    :param b:
    :param c:
    :return:
    """

    return a + b * n + c * np.log10(n)

# -----------------------------------------------------------------

def sort_lists(*args, **kwargs):

    """
    This function ...
    :param args:
    :return:
    """

    to_arrays = kwargs.pop("to_arrays", False)

    if to_arrays: return [np.array(list(t)) for t in zip(*sorted(zip(*args)))]
    else: return [list(t) for t in zip(*sorted(zip(*args)))]

# -----------------------------------------------------------------

def generate_scaling_behaviour(phase):

    """
    This function ...
    :param phase:
    :return:
    """

    behaviour = set()

    # Pure scaling term
    if phase in pure_scaling_behaviour:
        for component in pure_scaling_behaviour[phase]: behaviour.add(component)

    # Composite scaling term
    elif phase in composite_scaling_behaviour:
        for composite_phase in composite_scaling_behaviour[phase]:
            for component in generate_scaling_behaviour(composite_phase): behaviour.add(component)

    # Not recognized
    else: raise ValueError("Phase '" + phase + "' not recognized")

    # Return the scaling behaviour
    return list(behaviour)

# -----------------------------------------------------------------

class FitFunction(Callable):
#class FitFunction(FunctionType):

    def __init__(self, powers):

        self.powers = powers
        self.n = len(self.powers)

    def __call__(self, *args):
        x = args[0]
        result = 0.
        for i in range(self.n):
            a = args[i+1]
            if self.powers[i] == "log":
                result += a * np.log10(x)
            else: result += a * x**self.powers[i]
        return result

# -----------------------------------------------------------------

def generate_fit_function(behaviour):

    """
    This function ...
    :param behaviour:
    :return:
    """

    function = FitFunction(behaviour)

    nparameters = len(behaviour)

    return function.__call__, nparameters

# -----------------------------------------------------------------

def get_parallelization_mode(processes, threads, data_parallel):

    """
    This function ...
    :param processes:
    :param threads:
    :param data_parallel:
    :return:
    """

    if processes > 1:
        if threads > 1:
            if data_parallel: mode = "hybrid task+data (" + str(threads) + " threads)"
            else: mode = "hybrid task (" + str(threads) + " threads)"
        else:
            if data_parallel: mode = "multiprocessing task+data"
            else: mode = "multiprocessing task"
    else: mode = "multithreading"

    # Return the mode
    return mode

# -----------------------------------------------------------------
