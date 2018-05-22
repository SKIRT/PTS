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
import copy
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from collections import defaultdict
from collections import Callable
from string import ascii_lowercase
from matplotlib import rc

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..basics.measurement import Measurement
from ..basics.map import Map
from .timeline import create_timeline_plot
from ..basics.log import log
from ..tools import filesystem as fs
from ..basics.configurable import Configurable
from ..launch.timing import TimingTable
from ..launch.memory import MemoryTable
from ..extract.timeline import TimeLineExtractor, extract_timeline
from ..simulation.discover import SimulationDiscoverer
from ..basics.range import RealRange
from ..tools import tables
from ..tools import stringify
from ..tools.serialization import write_dict, load_dict, write_list, load_list
from ..tools.strings import alphabet, split_in_lines
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

rc('text', usetex=True)

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

phase_labels_timing = {"total": "Total runtime", "setup": "Setup time", "stellar": "Stellar runtime",
                    "spectra": "Runtime of dust spectra calculation", "dust": "Dust emission runtime",
                    "writing": "Writing time", "waiting": "Waiting time", "communication": "Communication time",
                    "dust densities communication": "Dust densities communication time",
                    "stellar absorption communication": "Stellar absorption table communication time",
                    "dust absorption communication": "Dust absorption table communication time",
                    "emission spectra communication": "Emission spectra communication time",
                    "instruments communication": "Instrument tables communication time",
                    "intermediate": "Intermediate time"}

# -----------------------------------------------------------------

# For legend, be concise
phase_names_for_legend = dict()
phase_names_for_legend["total"] = "simulation"
phase_names_for_legend["setup"] = "setup"
phase_names_for_legend["stellar"] = "stellar emission"
phase_names_for_legend["spectra"] = "dust emission spectra calculation"
phase_names_for_legend["dust"] = "dust emission"
phase_names_for_legend["writing"] = "writing"
phase_names_for_legend["waiting"] = "waiting"
phase_names_for_legend["communication"] = "communication"
phase_names_for_legend["dust densities communication"] = "dust densities"
phase_names_for_legend["stellar absorption communication"] = "absorbed stellar luminosities"
phase_names_for_legend["dust absorption communication"] = "absorbed dust luminosities"
phase_names_for_legend["emission spectra communication"] = "dust emission spectra"
phase_names_for_legend["instruments communication"] = "instruments"
phase_names_for_legend["intermediate"] = "intermediate"

# -----------------------------------------------------------------

# Phases that are parallel or serial in nature
parallel_phases = ["total", "stellar", "spectra", "dust"]
overhead_phases = ["waiting", "communication"]
serial_phases = ["writing", "setup"]

# -----------------------------------------------------------------

# Scaling properties
scaling_properties = ["runtime", "speedup", "efficiency", "CPU-time", "memory", "memory-gain", "total-memory", "timeline"]
scaling_properties_timing = ["runtime", "speedup", "efficiency", "CPU-time", "timeline"]
scaling_properties_memory = ["memory", "memory-gain", "total-memory"]

# -----------------------------------------------------------------

# Whether properties are expected to scale generally up or down
properties_behaviour = dict()
properties_behaviour["runtime"] = "descending"
properties_behaviour["speedup"] = "ascending"
properties_behaviour["efficiency"] = "descending"
properties_behaviour["CPU-time"] = "ascending"
properties_behaviour["memory"] = "descending"
properties_behaviour["memory-gain"] = "ascending"
properties_behaviour["total-memory"] = "ascending"

# Timing phases
simulation_phases_timing = ["total", "setup", "stellar", "spectra", "dust", "writing", "waiting", "communication", "intermediate"]

# Communication timing phases
communication_phases = ["dust densities communication", "stellar absorption communication", "dust absorption communication", "emission spectra communication", "instruments communication"]

# Memory phases
simulation_phases_memory = ["total", "setup", "stellar", "spectra", "dust", "writing"]

# All simulation phases
simulation_phases = list(set(simulation_phases_timing) | set(simulation_phases_memory))

# -----------------------------------------------------------------

phases_table_mapping_timing = dict()
phases_table_mapping_timing["total"] = "Total runtime"
phases_table_mapping_timing["setup"] = "Setup time"
phases_table_mapping_timing["stellar"] = "Stellar emission time"
phases_table_mapping_timing["spectra"] = "Spectra calculation time"
phases_table_mapping_timing["dust"] = "Dust emission time"
phases_table_mapping_timing["writing"] = "Writing time"
phases_table_mapping_timing["waiting"] = "Waiting time"
phases_table_mapping_timing["communication"] = "Communication time"
phases_table_mapping_timing["intermediate"] = "Intermediate time"
phases_table_mapping_timing["dust densities communication"] = "Dust densities communication time"
phases_table_mapping_timing["stellar absorption communication"] = "Stellar absorption communication time"
phases_table_mapping_timing["dust absorption communication"] = "Dust absorption communication time"
phases_table_mapping_timing["emission spectra communication"] = "Emission spectra communication time"
phases_table_mapping_timing["instruments communication"] = "Instruments communication time"

phases_table_mapping_memory = dict()
phases_table_mapping_memory["total"] = "Total peak memory"
phases_table_mapping_memory["setup"] = "Setup peak memory"
phases_table_mapping_memory["stellar"] = "Stellar emission peak memory"
phases_table_mapping_memory["spectra"] = "Spectra calculation peak memory"
phases_table_mapping_memory["dust"] = "Dust emission peak memory"
phases_table_mapping_memory["writing"] = "Writing peak memory"

# -----------------------------------------------------------------

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

dont_fit_phases = ["setup", "intermediate", "writing", "waiting"]

# -----------------------------------------------------------------

multiprocessing_phases = ["communication", "waiting"]
phases_not_relevant_for_multithreading = ["communication", "waiting"] # start from nprocesses = 2 or higher (otherwise we have runtimes of zero, doesn't work on a log-log plot)
phases_not_logaritmic_runtimes = ["communication", "waiting"]

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
pure_scaling_behaviour["communication"] = [0, 1, 2, "log"]

# Waiting
pure_scaling_behaviour["waiting"] = [0, 1, 2]

# Writing
pure_scaling_behaviour["writing"] = [0]

# Intermediate
pure_scaling_behaviour["intermediate"] = [0]

# -----------------------------------------------------------------

derived_scaling_behaviour = dict()

derived_scaling_behaviour["dust densities communication"] = "communication"
derived_scaling_behaviour["stellar absorption communication"] = "communication"
derived_scaling_behaviour["dust absorption communication"] = "communication"
derived_scaling_behaviour["emission spectra communication"] = "communication"
derived_scaling_behaviour["instruments communication"] = "communication"

# -----------------------------------------------------------------

composite_scaling_behaviour = dict()

# Total simulation
composite_scaling_behaviour["total"] = ("setup", "stellar", "dust", "spectra", "communication", "waiting", "writing")

# -----------------------------------------------------------------

class ScalingPlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ScalingPlotter, self).__init__(*args, **kwargs)

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

        # Fit parameters for timing data
        self.timing_fit_parameters = defaultdict(dict)
        self.timing_fit_parameter_errors = defaultdict(dict)

        # Fit parameters for memory data
        self.memory_fit_parameters = defaultdict(dict)

        # The number of cores used for the 'serial' timing and memory data
        self.serial_timing_ncores = defaultdict(dict)
        self.serial_memory_ncores = defaultdict(dict)

        # Parameter set mapping for the simulations
        self.distinct_parameter_sets_timing = []
        self.distinct_parameter_sets_memory = []

        self.parameters_mapping_timing = []
        self.parameters_mapping_memory = []

        # Entries in the timing table to be ignored
        self.ignore_timing_entries = []

        # Entries in the memory table to be ignored
        self.ignore_memory_entries = []

        # Parameter sets to be ignored
        self.ignore_parameter_sets_timing = set() # set of parameter sets (tuples)
        self.ignore_parameter_sets_memory = set()

        # Flag
        self.is_prepared = False

        # Different ski parameters in timing and memory data
        self._different_parameters_timing = None
        self._different_parameters_memory = None

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

    @lazyproperty
    def needs_timing(self):

        """
        This function ...
        :return:
        """

        for property in self.config.properties:
            if property in timing_properties: return True
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def needs_memory(self):

        """
        This function ...
        :return:
        """

        for property in self.config.properties:
            if property in memory_properties: return True
        return False

    # -----------------------------------------------------------------

    def has_timing_fit(self, phase, parameter_set, mode):

        """
        This function ...
        :param phase:
        :param parameter_set:
        :param mode:
        :return:
        """

        if phase in self.timing_fit_parameters:
            if parameter_set in self.timing_fit_parameters[phase]:
                if mode in self.timing_fit_parameters[phase][parameter_set]: return True
                else: return False
            else: return False
        else: return False

    # -----------------------------------------------------------------

    def has_memory_fit(self, phase, parameter_set, mode):

        """
        This function ...
        :param phase:
        :param parameter_set:
        :param mode:
        :return:
        """

        if phase in self.memory_fit_parameters:
            if parameter_set in self.memory_fit_parameters[phase]:
                if mode in self.memory_fit_parameters[phase][parameter_set]: return True
                else: return False
            else: return False
        else: return False

    # -----------------------------------------------------------------

    @property
    def do_prepare(self):

        """
        This function ...
        :return:
        """

        return not self.is_prepared

    # -----------------------------------------------------------------

    @property
    def do_fit(self):

        """
        This function ...
        :return:
        """

        return not self.config.hybridisation and self.config.fit

    # -----------------------------------------------------------------

    @property
    def do_write(self):

        """
        This function ...
        :return:
        """

        return self.config.output is not None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Prepare data into plottable format
        if self.do_prepare: self.prepare()

        # 3. Do fitting
        if self.do_fit: self.fit()

        # 4. Write
        if self.do_write: self.write()

        # 5. Plot
        self.plot()

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

        # All timing
        if self.config.all_timing:
            self.config.properties = scaling_properties_timing
            self.config.phases = simulation_phases_timing

        # All memory
        if self.config.all_memory:
            self.config.properties = scaling_properties_memory
            self.config.phases = simulation_phases_memory

        # Check for 'all' flag
        if self.config.all:
            self.config.properties = scaling_properties
            self.config.phases = simulation_phases

        # If timeline is asked as property, enable all phases
        if "timeline" in self.config.properties: self.config.all_phases = True

        # Check for 'all_phases' flag
        if self.config.all_phases: self.config.phases = simulation_phases

        # Check for 'all_properties' flag
        if self.config.all_properties: self.config.properties = scaling_properties

        # Get timing table
        self.get_timing(**kwargs)

        # Get memory table
        self.get_memory(**kwargs)

        # If data input path is given
        if self.config.data_input is not None: self.load_data_input()

        # No input data is given and timing or memory table not obtained
        elif (not self.has_timing) or (not self.has_memory): self.extract_data(**kwargs)

        # Even out properties in the timing and memory tables
        self.even_out_ski_properties()

        # Check ski parameters
        if len(self.simulations) > 0:

            # Check the ski file parameters
            if self.has_multiski: self.check_parameters()
            else:
                #parameters = comparison_parameters_from_ski(self.simulations[0].ski_path, self.simulations[0].input_path)
                parameters = self.timing.ski_parameters_for_simulation(self.timing.simulation_names[0])
                self.distinct_parameter_sets_timing.append(parameters)

                parameters = self.memory.ski_parameters_for_simulation(self.memory.simulation_names[0])
                self.distinct_parameter_sets_memory.append(parameters)

        # Check consistency between the timing and memory tables
        if self.config.tolerant: self.check_consistency_tolerant()
        else: self.check_consistency_strict()

        # Check the coverage of the timing and memory data
        self.check_coverage()

    # -----------------------------------------------------------------

    def get_timing(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get table
        if "timing" in kwargs: self.timing = kwargs.pop("timing")

        # Load table
        elif self.config.timing_table is not None: self.timing = TimingTable.from_file(self.config.timing_table)

    # -----------------------------------------------------------------

    @property
    def has_timing(self):

        """
        This function ...
        :return:
        """

        return self.timing is not None

    # -----------------------------------------------------------------

    def get_memory(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get table
        if "memory" in kwargs: self.memory = kwargs.pop("memory")

        # If path of memory table is given
        elif self.config.memory_table is not None: self.memory = MemoryTable.from_file(self.config.memory_table)

    # -----------------------------------------------------------------

    @property
    def has_memory(self):

        """
        This function ...
        :return:
        """

        return self.memory is not None

    # -----------------------------------------------------------------

    def load_data_input(self):

        """
        This function ...
        :return:
        """

        # Load timing data
        timing_data_path = fs.join(self.config.data_input, "timing_data.dat")
        if fs.is_file(timing_data_path):
            self.timing_data = load_dict(timing_data_path)
            self.serial_timing = load_dict(fs.join(self.config.data_input, "serial_timing_data.dat"))
            self.serial_timing_ncores = load_dict(fs.join(self.config.data_input, "serial_timing_ncores.dat"))
            self._different_parameters_timing = load_list(fs.join(self.config.data_input, "different_parameters_timing.dat"))
            self.ignore_parameter_sets_timing = load_list(fs.join(self.config.data_input, "ignore_parameter_sets_timing.dat"))

            # Set flag
            self.is_prepared = True

        # Load memory data
        memory_data_path = fs.join(self.config.data_input, "memory_data.dat")
        if fs.is_file(memory_data_path):
            self.memory_data = load_dict(fs.join(self.config.data_input, "memory_data.dat"))
            self.serial_memory = load_dict(fs.join(self.config.data_input, "serial_memory_data.dat"))
            self.serial_memory_ncores = load_dict(fs.join(self.config.data_input, "serial_memory_ncores.dat"))
            self._different_parameters_memory = load_list(fs.join(self.config.data_input, "different_parameters_memory.dat"))
            self.ignore_parameter_sets_memory = load_list(fs.join(self.config.data_input, "ignore_parameter_sets_memory.dat"))

            # Set flag
            self.is_prepared = True

        # Data not prepared yet?
        if self.is_prepared: return
        else:

            timing_path = fs.join(self.config.data_input, "timing.dat")
            memory_path = fs.join(self.config.data_input, "memory.dat")

            # Load timing and/or memory table
            if fs.is_file(timing_path): self.timing = TimingTable.from_file(timing_path)
            if fs.is_file(memory_path): self.memory = MemoryTable.from_file(memory_path)

    # -----------------------------------------------------------------

    def extract_data(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # If simulations are passed
        if "simulations" in kwargs: self.simulations = kwargs.pop("simulations")

        # If simulations have been added
        elif len(self.simulations) > 0: pass

        # Load simulations from working directory if none have been added
        else: self.load_simulations()

        # Do extraction
        self.extract()

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
        discoverer.config.directories = self.config.input
        discoverer.config.list = self.config.report_simulations

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
        if extract_timing: self.timing = TimingTable()

        # Initialize a memory table
        if extract_memory: self.memory = MemoryTable()

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

                # Get the timeline
                timeline = extract_timing(simulation)

                # Add an entry to the timing table
                unique_name = self.timing.add_from_simulation(simulation, ski, log_file, timeline, parameters=parameters)

                # Change the simulation name
                simulation.name = unique_name

            # If memory has to be extracted
            if extract_memory:

                # Add an entry to the memory table
                unique_name = self.memory.add_from_simulation(simulation, ski, log_file, parameters=parameters)

                # Shouldn't happen
                if unique_name != simulation.name: raise RuntimeError("Something went wrong")

    # -----------------------------------------------------------------

    def even_out_ski_properties(self):

        """
        This function ...
        :return:
        """

        # Even the number of packages if they are equivalent
        self.timing.even_npackages()

    # -----------------------------------------------------------------

    def check_parameters(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Checking simulation parameters ...")

        # Loop over the simulations
        for simulation_name in self.timing.simulation_names:

            # Get the parameters relevant for the timing
            parameters = self.timing.ski_parameters_for_simulation(simulation_name)

            # Get index
            index = None
            for i, dps in enumerate(self.distinct_parameter_sets_timing):
                if dps == parameters:
                    index = i
                    break

            if index is None:
                self.distinct_parameter_sets_timing.append(parameters)
                index = len(self.distinct_parameter_sets_timing) - 1

            # Set index for the mapping: NOT USED
            self.parameters_mapping_timing.append(index)

        # Loop over the simulations
        for simulation_name in self.memory.simulation_names:

            # Get the parameters relevant for the memory consumption
            parameters = self.memory.ski_parameters_for_simulation(simulation_name)

            # Get index
            index = None
            for i, dps in enumerate(self.distinct_parameter_sets_memory):
                if dps == parameters:
                    index = i
                    break

            if index is None:
                self.distinct_parameter_sets_memory.append(parameters)
                index = len(self.distinct_parameter_sets_memory) - 1

            # Set index for the mapping: NOT USED
            self.parameters_mapping_memory.append(index)

    # -----------------------------------------------------------------

    def check_consistency_tolerant(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Checking whether the timing and memory tables are consistent ...")

        # Check whether both tables have the same length
        if len(self.timing) != len(self.memory): log.warning("Timing and memory tables do not have the same length")

        # Loop over the simulation names in the timing table
        for simulation_name in self.timing.simulation_names:

            # Check whether this simulation is in the memory table as well
            if simulation_name not in self.memory["Simulation name"]:
                log.warning("Simulation '" + simulation_name + "' not present in the memory table")
                continue

            # Get the timing parameters
            timing_parameters = self.timing.ski_parameters_for_simulation(simulation_name, which="all")

            # "Wavelengths", "Packages", "Dust cells", "Grid type", "Min level", "Max level",
            # "Search method", "Sample count", "Max optical depth", "Max mass fraction",
            # "Max density dispersion", "Self-absorption", "Transient heating"

            # Get the memory parameters
            memory_parameters = self.memory.ski_parameters_for_simulation(simulation_name, which="all")

            # "Wavelengths", "Dust cells", "Grid type", "Min level", "Max level", "Search method",
            # "Sample count", "Max optical depth", "Max mass fraction", "Max density dispersion",
            # "Self-absorption", "Transient heating", "Number of pixels"

            #print(timing_parameters)
            #print(memory_parameters)

            # Compare the number of wavelengths
            if timing_parameters["Wavelengths"] != memory_parameters["Wavelengths"]: raise RuntimeError("Number of wavelengths not equal for simulation '" + simulation_name + "'")

            # Compare the number of dust cells
            if timing_parameters["Dust cells"] != memory_parameters["Dust cells"]: raise RuntimeError("Number of dust cells not equal for simulation '" + simulation_name + "'")

            # Compare the min level, max level, search method, sample count, max optical depth, max mass fraction, max density dispersion

            # Compare self-absorption flag

            # Compare transient heating flag

    # -----------------------------------------------------------------

    def check_consistency_strict(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Checking whether the timing and memory tables are consistent ...")

        # Get simulation names in timing and memory table
        timing_simulation_names = sorted(self.timing.simulation_names)
        memory_simulation_names = sorted(self.memory.simulation_names)

        # Get the differences
        result = cmp(timing_simulation_names, memory_simulation_names)

        # Check
        if result != 0: raise RuntimeError("The timing and memory tables are not consistent")

    # -----------------------------------------------------------------

    def check_coverage(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Checking the coverage of the data ...")

        # Timing
        self.check_coverage_timing()

        # Memory
        self.check_coverage_memory()

    # -----------------------------------------------------------------

    def check_coverage_timing(self):

        """
        This function ...
        :return:
        """

        # If simulations with different parameters are used
        if self.has_multiparameters_timing:

            # Loop over the distinct parameter sets
            for parameter_set in self.distinct_parameter_sets_timing:

                # Find the entries in the timing table that correspond to this parameter set
                indices = self.timing.indices_for_parameters(parameter_set)

                remove_indices = False

                if self.config.hybridisation:
                    if np.min(self.timing["Processes"][indices] == np.max(self.timing["Processes"][indices])):
                        log.warning("All runtimes are generated with the same number of processes: skipping this parameter set for the timing data")
                        remove_indices = True

                else:
                    if np.min(self.timing["Cores"][indices] == np.max(self.timing["Cores"][indices])):
                        log.warning("All runtimes are generated on the same number of cores: skipping this parameter set for the timing data")
                        remove_indices = True

                # Set to ignore
                if remove_indices: self.ignore_timing_entries += indices

        #  Only simulations with the same parameters (of those relevant for the timings)
        else:

            # Check if the data spans multiple number of cores
            if self.config.hybridisation:
                if np.min(self.timing["Processes"]) == np.max(self.timing["Processes"]): raise RuntimeError("All runtimes are generated with the same number of processes, you cannot run with --hybridisation")
            else:
                if np.min(self.timing["Cores"]) == np.max(self.timing["Cores"]): raise RuntimeError("All runtimes are generated on the same number of cores. Run with --hybridisation")

    # -----------------------------------------------------------------

    def check_coverage_memory(self):

        """
        This function ...
        :return:
        """

        # If simulations with different parameters are used
        if self.has_multiparameters_memory:

            # Loop over the distinct parameter sets
            for parameter_set in self.distinct_parameter_sets_memory:

                # Find the entries in the memory table that correspond to this parameter set
                indices = self.memory.indices_for_parameters(parameter_set)

                remove_indices = False

                if self.config.hybridisation:
                    if np.min(self.memory["Processes"][indices] == np.max(self.memory["Processes"][indices])):
                        log.warning("All memory data are generated with the same number of processes: skipping this parameter set for the memory data")
                        remove_indices = True
                else:
                    if np.min(self.memory["Cores"][indices] == np.max(self.memory["Cores"][indices])):
                        log.warning("All runtimes are generated on the same number of cores: skipping this parameter set for the memory data")
                        remove_indices = True

                # Set to ignore
                if remove_indices: self.ignore_memory_entries += indices

        # Only simulations with the same parameters (of those relevant for the memory usage)
        else:

            # Check if the data spans multiple number of cores
            if self.config.hybridisation:
                if np.min(self.memory["Processes"]) == np.max(self.memory["Processes"]): raise RuntimeError("All memory data are generated with the same number of processes, you cannot run with --hybridisation")
            else:
                if np.min(self.memory["Cores"]) == np.max(self.memory["Cores"]): raise RuntimeError("All memory data are generated on the same number of cores. Run with --hybridisation")

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
    def has_multiparameters_timing(self):

        """
        This function ...
        :return:
        """

        return len(self.distinct_parameter_sets_timing) > 1

    # -----------------------------------------------------------------

    @property
    def has_multiparameters_memory(self):

        """
        This function ...
        :return:
        """

        return len(self.distinct_parameter_sets_memory) > 1

    # -----------------------------------------------------------------

    @lazyproperty
    def different_parameters_timing(self):

        """
        This function ...
        :return:
        """

        if self._different_parameters_timing is None: self._different_parameters_timing = self.timing.different_ski_parameters()
        return self._different_parameters_timing

    # -----------------------------------------------------------------

    @lazyproperty
    def different_parameters_memory(self):

        """
        This function ...
        :return:
        """

        if self._different_parameters_memory is None: self._different_parameters_memory = self.memory.different_ski_parameters()
        return self._different_parameters_memory

    # -----------------------------------------------------------------

    @property
    def different_parameters_timing_with_ignored(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def different_parmeters_memory_with_ignored(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def parameters_timing_left_after_ignored(self):

        """
        This function ...
        :return:
        """

        parameter_values = defaultdict(list)

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets_timing:

            if parameter_set in self.ignore_parameter_sets_timing: continue

            if parameter_set == "empty": continue

            if isinstance(parameter_set, tuple):
                for index, name in enumerate(self.different_parameters_timing):
                    parameter_values[name].append(parameter_set[index])
            else: parameter_values[self.different_parameters_timing[0]].append(parameter_set)

        return [key for key in parameter_values if not all_equal(parameter_values[key])]

    # -----------------------------------------------------------------

    @lazyproperty
    def parameters_memory_left_after_ignored(self):

        """
        This function ...
        :return:
        """

        parameter_values = defaultdict(list)

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets_memory:

            if parameter_set in self.ignore_parameter_sets_memory: continue

            if parameter_set == "empty": continue

            if isinstance(parameter_set, tuple):
                for index, name in enumerate(self.different_parameters_memory):
                    parameter_values[name].append(parameter_set[index])
            else: parameter_values[self.different_parameters_memory[0]].append(parameter_set)

        return [key for key in parameter_values if not all_equal(parameter_values[key])]

    # -----------------------------------------------------------------

    def get_parameter_set_timing(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get relevant simulation properties
        #relevant_parameters = Map()
        #for parameter in self.different_parameters_timing:
        #    relevant_parameters[parameter] = self.timing[parameter][i]
        #if len(relevant_parameters) == 0: relevant_parameters = "empty"

        values = []
        for parameter in self.different_parameters_timing:
            value = self.timing[parameter][index]
            if hasattr(value, "dtype") and np.issubdtype(str(value.dtype), np.string_): value = str(value) # dtype('S..') to str
            values.append(value)

        if len(values) == 0: values = "empty"
        elif len(values) == 1: values = values[0]
        else: values = tuple(values)

        return values

    # -----------------------------------------------------------------

    def get_parameter_set_memory(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        values = []

        for parameter in self.different_parameters_memory:
            value = self.memory[parameter][index]
            if hasattr(value, "dtype") and np.issubdtype(str(value.dtype), np.string_): value = str(value) # dtype('S..') to str
            values.append(value)

        if len(values) == 0: values = "empty"
        elif len(values) == 1: values = values[0]
        else: values = tuple(values)

        return values

    # -----------------------------------------------------------------

    def get_parameter_index_timing(self, parameter_name):

        """
        This function ...
        :param parameter_name:
        :return:
        """

        return self.different_parameters_timing.index(parameter_name)

    # -----------------------------------------------------------------

    def get_parameter_index_memory(self, parameter_name):

        """
        This function ...
        :param parameter_name:
        :return:
        """

        return self.different_parameters_memory.index(parameter_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_phases_timing(self):

        """
        This function ...
        :return:
        """

        return [phase for phase in simulation_phases_timing if phase in self.config.phases]

    # -----------------------------------------------------------------

    @lazyproperty
    def parallel_simulation_phases_timing(self):

        """
        This function ...
        :return:
        """

        return [phase for phase in self.simulation_phases_timing if phase in parallel_phases]

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_phases_and_subphases_timing(self):

        """
        This function ...
        :return:
        """

        phases = copy.deepcopy(self.simulation_phases_timing)
        if "communication" in phases and self.config.split_communication:
            phases += self.config.communication_subphases
        return phases

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_phases_memory(self):

        """
        This function ...
        :return:
        """

        return [phase for phase in simulation_phases_memory if phase in self.config.phases]

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_phases_and_subphases_memory(self):

        """
        This function ...
        :return:
        """

        phases = copy.deepcopy(self.simulation_phases_memory)
        return phases

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
        parameters_modes_processor_counts_dict_timing = defaultdict(lambda: defaultdict(set))

        # Keep track of the different processor counts encountered for the different parameter sets and parallelization modes
        parameters_modes_processor_counts_dict_memory = defaultdict(lambda: defaultdict(set))

        # Create dictionaries to contain the data before it is averaged over the different simulations
        times = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

        # Create dictionaries for the memory data
        memorys = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

        # Number of cores for the serial memory runs
        serial_memory_ncores = dict()

        # Debugging
        #print(self.distinct_parameter_sets_timing)
        #print(self.different_parameters_timing)
        if self.has_multiparameters_timing: log.debug("The parameters that differ between the simulations that are relevant for the timing are: " + ", ".join(self.different_parameters_timing))
        else: log.debug("There are no parameters that differ between the simulations that are relevant for the timing")

        # Debugging
        #print(self.distinct_parameter_sets_memory)
        #print(self.different_parameters_memory)
        if self.has_multiparameters_memory: log.debug("The parameters that differ between the simulations that are relevant for the memory usage are: " + ", ".join(self.different_parameters_memory))
        else: log.debug("There are no parameters that differ between the simulations that are relevant for the memory data")

        # Loop over the different entries in the timing table
        for i in range(len(self.timing)):

            # Get the simulation name and the corresponding index in the memory table
            simulation_name = self.timing["Simulation name"][i]
            j = self.memory.index_for_simulation(simulation_name)
            if j is None: continue
            #print(i, j)

            # Get the relevant timing parameters
            timing_parameters = self.get_parameter_set_timing(i)
            #print(timing_parameters)

            # Ignore parameter set
            if i in self.ignore_timing_entries: self.ignore_parameter_sets_timing.add(timing_parameters)

            # Get the relevant memory parameters
            memory_parameters = self.get_parameter_set_memory(j)
            if j in self.ignore_memory_entries: self.ignore_parameter_sets_memory.add(memory_parameters)
            #print(memory_parameters)

            # Get the number of processes and threads
            processes = int(self.timing["Processes"][i])
            threads_per_core = int(self.timing["Threads per core"][i])
            cores_per_process = int(self.timing["Cores"][i] / processes)
            threads = threads_per_core * cores_per_process
            processors = processes * threads

            # Data parallelization enabled or not
            data_parallel = self.timing["Data-parallel"][i]

            # Determine the parallelization mode
            if self.config.hybridisation: mode = get_hybridisation_mode(self.timing["Cores"][i], data_parallel)
            else: mode = get_parallelization_mode(processes, threads, data_parallel)

            # Skip certain modes if requested
            if not self.config.hybridisation:

                # Skip certain modes if requested
                simple_mode = mode.split(" ")[0]
                if simple_mode not in self.config.modes: continue

                # Skip task+data or task parallel mode if requested
                if mode.startswith("multiprocessing") or mode.startswith("hybrid"):

                    if "task+data" in mode:
                        if not self.config.use_task_data_parallel: continue
                    else:
                        if not self.config.use_task_parallel: continue

            # Get the runtimes per phase
            time_per_phase = dict()
            if self.needs_timing:
                for phase in self.simulation_phases_and_subphases_timing: time_per_phase[phase] = self.timing[phases_table_mapping_timing[phase]][i]

            # Get the memory usage per phase
            memory_per_phase = dict()
            if self.needs_memory:
                for phase in self.simulation_phases_and_subphases_memory: memory_per_phase[phase] = self.memory[phases_table_mapping_memory[phase]][j]

            if self.needs_timing:

                # If the number of processors is 1, add the runtimes for the different simulation phases to the
                # dictionary that contains the serial runtimes
                if (self.config.hybridisation and processes == 1) or (not self.config.hybridisation and processors == 1):

                    for phase in self.simulation_phases_and_subphases_timing: serial_times[timing_parameters][phase].append(time_per_phase[phase])

            if self.needs_memory:

                # Number of processes = 1: equivalent to 'serial' in terms of memory consumption
                if processes == 1:

                    for phase in self.simulation_phases_and_subphases_memory: serial_memory[memory_parameters][phase].append(memory_per_phase[phase])
                    serial_memory_ncores[memory_parameters] = processors

            # Determine variable
            processes_or_processors = processes if self.config.hybridisation else processors

            # Add the processor count for this parallelization mode
            parameters_modes_processor_counts_dict_timing[timing_parameters][mode].add(processes_or_processors)
            parameters_modes_processor_counts_dict_memory[memory_parameters][mode].add(processes_or_processors)

            # Fill in the runtimes at the appropriate place in the dictionaries
            for phase in time_per_phase: times[phase][timing_parameters][mode][processes_or_processors].append(time_per_phase[phase])

            # Fill in the memory usage at the appropriate place in the dictionaries
            for phase in memory_per_phase: memorys[phase][memory_parameters][mode][processes_or_processors].append(memory_per_phase[phase])

        if self.needs_timing:

            # Average the serial runtimes, loop over each parameter set and phase
            for parameter_set in serial_times:
                for phase in serial_times[parameter_set]:

                    self.serial_timing[parameter_set][phase].time = np.mean(serial_times[parameter_set][phase])
                    self.serial_timing[parameter_set][phase].error = self.config.sigma_level * np.std(serial_times[parameter_set][phase])
                    self.serial_timing_ncores[parameter_set][phase] = 1

        if self.needs_memory:

            # Average the serial memory usages, loop over each phase
            for parameter_set in serial_memory:
                for phase in serial_memory[parameter_set]:

                    self.serial_memory[parameter_set][phase].memory = np.mean(serial_memory[parameter_set][phase])
                    self.serial_memory[parameter_set][phase].error = self.config.sigma_level * np.std(serial_memory[parameter_set][phase])
                    self.serial_memory_ncores[parameter_set][phase] = serial_memory_ncores[parameter_set]

        ## TIMING

        if self.needs_timing:

            # Loop over all parameter sets of the timing data
            for parameter_set in parameters_modes_processor_counts_dict_timing:

                # Loop over all encountered parallelization modes
                modes = parameters_modes_processor_counts_dict_timing[parameter_set]
                for mode in modes:

                    # Loop over all processor counts encountered for this mode
                    for processors in modes[mode]:

                        # Average the runtimes for the different simulation phases for the different
                        # runs for a certain parallelization mode and number of processors (or processes)
                        for phase in times:

                            self.timing_data[phase][parameter_set][mode].processor_counts.append(processors)
                            self.timing_data[phase][parameter_set][mode].times.append(np.mean(times[phase][parameter_set][mode][processors]))
                            self.timing_data[phase][parameter_set][mode].errors.append(self.config.sigma_level * np.std(times[phase][parameter_set][mode][processors]))

        ## MEMORY

        if self.needs_memory:

            # Loop over all parameter sets of the memory data
            for parameter_set in parameters_modes_processor_counts_dict_memory:

                # Loop over all encountered parallelization modes
                modes = parameters_modes_processor_counts_dict_memory[parameter_set]
                for mode in modes:

                    # Loop over all processor counts encountered for this mode
                    for processors in modes[mode]:

                        # Average the memory usage for the different simulation phases for the different
                        # runs for a certain parallelization mode and number of processors (or processes)
                        for phase in memorys:

                            self.memory_data[phase][parameter_set][mode].processor_counts.append(processors)
                            self.memory_data[phase][parameter_set][mode].memory.append(np.mean(memorys[phase][parameter_set][mode][processors]))
                            self.memory_data[phase][parameter_set][mode].errors.append(self.config.sigma_level * np.std(memorys[phase][parameter_set][mode][processors]))

        # Check whether we have data
        self.check_data()

        # Set equivalent data
        self.set_equivalent_data(parameters_modes_processor_counts_dict_timing, parameters_modes_processor_counts_dict_memory)

        # Set missing data
        self.set_missing_data(parameters_modes_processor_counts_dict_timing, parameters_modes_processor_counts_dict_memory)

        # Check coverage of data in the different modes
        self.check_coverage_modes()

        # Add communication and waiting times of zero for 1 process
        if not self.config.hybridisation and self.needs_timing: self.add_fixed_multiprocessing_phases_runtimes()

        # Set flag
        self.is_prepared = True

    # -----------------------------------------------------------------

    def check_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking whether data is present ...")

        # Check whether we have timing data
        if self.needs_timing: self.check_timing_data()

        # Check whether we have memory data
        if self.needs_memory: self.check_memory_data()

    # -----------------------------------------------------------------

    def check_timing_data(self):

        """
        This function ...
        :return:
        """

        if len(self.timing_data) == 0: raise RuntimeError("We have no timing data")

    # -----------------------------------------------------------------

    def check_memory_data(self):

        """
        This function ...
        :return:
        """

        if len(self.memory_data) == 0: raise RuntimeError("We have no memory data")

    # -----------------------------------------------------------------

    def set_equivalent_data(self, parameters_modes_processor_counts_dict_timing, parameters_modes_processor_counts_dict_memory):

        """
        This function ...
        :param parameters_modes_processor_counts_dict_timing:
        :param parameters_modes_processor_counts_dict_memory:
        :return:
        """

        # Inform the user
        log.info("Setting equivalent data ...")

        # Set equivalent timing data
        if self.needs_timing: self.set_equivalent_timing_data(parameters_modes_processor_counts_dict_timing)

        # Set equivalent memory data
        if self.needs_memory: self.set_equivalent_memory_data(parameters_modes_processor_counts_dict_memory)

    # -----------------------------------------------------------------

    def set_equivalent_timing_data(self, parameters_modes_processor_counts_dict_timing):

        """
        This function ...
        :param parameters_modes_processor_counts_dict_timing:
        :return:
        """

        # Inform the user
        log.info("Setting equivalent timing data ...")

        # Between hybrid task and hybrid task+data
        # (in hybridisation plotting mode the single-process end of the curve is not labeled as
        # 'multithreading mode' but just as xx cores [task/task+data])
        if self.config.hybridisation: self.set_equivalent_timing_data_hybrid(parameters_modes_processor_counts_dict_timing)

        # From multithreading to multiprocessing and hybrid
        else: self.set_equivalent_timing_data_from_multithreading(parameters_modes_processor_counts_dict_timing)

    # -----------------------------------------------------------------

    def set_equivalent_timing_data_from_multithreading(self, parameters_modes_processor_counts_dict_timing):

        """
        This function ...
        :param parameters_modes_processor_counts_dict_timing:
        :return:
        """

        # Debugging
        log.debug("Setting equivalent timing data from multithreading data ...")

        # Loop over all parameter sets of the timing data
        for parameter_set in parameters_modes_processor_counts_dict_timing:

            # Get all parallelization modes
            modes = parameters_modes_processor_counts_dict_timing[parameter_set]

            #print(parameter_set, modes)

            # Does multithreading mode exist in the modes?
            if not "multithreading" in modes: continue

            # Get a list of the other modes
            other_modes = [mode for mode in modes if mode != "multithreading"]

            # If there are no other modes, there is no point in this function
            if len(other_modes) == 0: continue

            # Get the different processor counts for multithreading mode
            processor_counts = list(modes["multithreading"])

            # SEARCH FOR SINGLE-THREADING
            # Check whether serial (singlethreading) is in the multithreading data
            if 1 in processor_counts:

                # Get the index of the serial run in the times and errors lists
                index = processor_counts.index(1)
                # Should be the same as:
                # self.timing_data["total"][parameter_set][mode].processor_counts.index(1)
                # with "total" replaced by any phase
                # Asser this:
                assert index == self.timing_data["total"][parameter_set]["multithreading"].processor_counts.index(1)

                # Loop over all phases
                for phase in self.simulation_phases_and_subphases_timing:

                    # Get the serial time and error from the multithreading mode
                    time = self.timing_data[phase][parameter_set]["multithreading"].times[index]
                    error = self.timing_data[phase][parameter_set]["multithreading"].errors[index]

                    # Set this time and error for the pure multiprocessing modes
                    for mode in other_modes:

                        # Skip modes that are not purily multiprocessing
                        if not mode.startswith("multiprocessing"): continue

                        # Fill in the processor count of 1, the time and the error on the time
                        self.timing_data[phase][parameter_set][mode].processor_counts.append(1)
                        self.timing_data[phase][parameter_set][mode].times.append(time)
                        self.timing_data[phase][parameter_set][mode].errors.append(error)

            # ALSO FILL IN MULT-THREADING DATA POINTS FOR HYBRID CURVES
            # Loop over the hybrid modes
            for mode in other_modes:

                # Skip not-hybrid modes
                if "hybrid" not in mode: continue

                # Get the number of thread per process
                nthreads = nthreads_per_process_for_hybrid_mode(mode)

                # Check whether there is multithreading data for this number of threads
                if nthreads not in processor_counts: continue

                # Get the index for the multithreading entry for this number of threads
                index = processor_counts.index(nthreads)

                # Loop over all phases
                for phase in self.simulation_phases_and_subphases_timing:

                    # Get the serial time and error from the multithreading mode
                    time = self.timing_data[phase][parameter_set]["multithreading"].times[index]
                    error = self.timing_data[phase][parameter_set]["multithreading"].errors[index]

                    # Fill in the processor count of 1, the time and the error on the time
                    self.timing_data[phase][parameter_set][mode].processor_counts.append(nthreads)
                    self.timing_data[phase][parameter_set][mode].times.append(time)
                    self.timing_data[phase][parameter_set][mode].errors.append(error)

    # -----------------------------------------------------------------

    def set_equivalent_timing_data_hybrid(self, parameters_modes_processor_counts_dict_timing):

        """
        This function ...
        :return:
        """

        # Loop over all parameter sets of the timing data
        for parameter_set in parameters_modes_processor_counts_dict_timing:

            # Get all parallelization modes
            modes = parameters_modes_processor_counts_dict_timing[parameter_set]

            # Does a task based hybrid mode exist in the modes?
            reference_modes = dict()
            reference_indices = dict()
            for mode in modes:

                # Get processor counts
                processor_counts = list(parameters_modes_processor_counts_dict_timing[parameter_set][mode])

                # Search for single-processing
                if 1 not in processor_counts: continue
                index = processor_counts.index(1)

                if "(task)" in mode:

                    #time = self.timing_data[phase][parameter_set]["multithreading"].times[index]
                    #error = self.timing_data[phase][parameter_set]["multithreading"].errors[index]

                    ncores = ncores_from_hybridisation_mode(mode)
                    reference_modes[ncores] = mode
                    reference_indices[ncores] = index

            # No reference task modes
            if len(reference_modes) == 0: continue

            # Loop over the reference modes
            for ncores in reference_modes:

                # Get mode
                reference_mode = reference_modes[ncores]

                # Get index for single-process entry
                index = reference_indices[ncores]

                # Search for task+data mode with the same number of cores
                for mode in modes:

                    # Check number of cores and that it is a task+data mode
                    if ncores == ncores_from_hybridisation_mode(mode) and "(task+data)" in mode:

                        # Loop over all phases
                        for phase in self.simulation_phases_and_subphases_timing:

                            # Get time and error
                            time = self.timing_data[phase][parameter_set][reference_mode].times[index]
                            error = self.timing_data[phase][parameter_set][reference_mode].errors[index]

                            # Fill in the process count of 1, the time and the error on the time
                            self.timing_data[phase][parameter_set][mode].processor_counts.append(1)
                            self.timing_data[phase][parameter_set][mode].times.append(time)
                            self.timing_data[phase][parameter_set][mode].errors.append(error)

    # -----------------------------------------------------------------

    def set_equivalent_memory_data(self, parameters_modes_processor_counts_dict_memory):

        """
        This function ...
        :param parameters_modes_processor_counts_dict_memory:
        :return:
        """

        # Loop over all parameter sets of the memory data
        for parameter_set in parameters_modes_processor_counts_dict_memory:

            # Loop over all encountered parallelization modes
            modes = parameters_modes_processor_counts_dict_memory[parameter_set]

            # Does multithreading mode exist in the modes?
            if not "multithreading" in modes: continue

            # Get a list of the other modes
            other_modes = [mode for mode in modes if mode != "multithreading"]

            # If there are no other modes, there is no point in this function
            if len(other_modes) == 0: continue

            # Get the different processor counts for multithreading mode
            processor_counts = list(modes["multithreading"])

            # SEARCH FOR SINGLE-THREADING
            # Check whether serial (singlethreading) is in the multithreading data
            if 1 in processor_counts:

                # Get the index of the serial run in the memory usage and errors lists
                index = processor_counts.index(1)
                # Should be the same as:
                # self.timing_data["total"][parameter_set][mode].processor_counts.index(1)
                # with "total" replaced by any phase
                # Asser this:
                assert index == self.memory_data["total"][parameter_set]["multithreading"].processor_counts.index(1)

                # Loop over all phases
                for phase in self.simulation_phases_and_subphases_memory:

                    # Get the serial memory usage and error from the multithreading mode
                    memory = self.memory_data[phase][parameter_set]["multithreading"].memory[index]
                    error = self.memory_data[phase][parameter_set]["multithreading"].errors[index]

                    # Set this memory usage and error for the pure multiprocessing modes
                    for mode in other_modes:

                        # Skip modes that are not purily multiprocessing
                        if mode.startswith("multiprocessing"): continue

                        # Fill in the processor count of 1, the memory data and the error on the time
                        self.memory_data[phase][parameter_set][mode].processor_counts.append(1)
                        self.memory_data[phase][parameter_set][mode].memory.append(memory)
                        self.memory_data[phase][parameter_set][mode].errors.append(error)

            # ALSO FILL IN MULT-THREADING DATA POINTS FOR HYBRID CURVES
            # Loop over the hybrid modes
            for mode in other_modes:

                # Skip not-hybrid modes
                if "hybrid" not in mode: continue

                # Get the number of thread per process
                nthreads = nthreads_per_process_for_hybrid_mode(mode)

                # Check whether there is multithreading data for this number of threads
                if nthreads not in processor_counts: continue

                # Get the index for the multithreading entry for this number of threads
                index = processor_counts.index(nthreads)

                # Loop over all phases
                for phase in self.simulation_phases_and_subphases_memory:

                    # Get the serial memory usage and error from the multithreading mode
                    memory = self.memory_data[phase][parameter_set]["multithreading"].memory[index]
                    error = self.memory_data[phase][parameter_set]["multithreading"].errors[index]

                    # Fill in the processor count of 1, the memory usage and the error on the memory usage
                    self.memory_data[phase][parameter_set][mode].processor_counts.append(nthreads)
                    self.memory_data[phase][parameter_set][mode].memory.append(memory)
                    self.memory_data[phase][parameter_set][mode].errors.append(error)

    # -----------------------------------------------------------------

    def has_serial_timing_for_parameters(self, parameter_set):

        """
        This function ...
        :param parameter_set:
        :return:
        """

        return parameter_set in self.serial_timing

    # -----------------------------------------------------------------

    def set_missing_data(self, parameters_modes_processor_counts_dict_timing, parameters_modes_processor_counts_dict_memory):

        """
        This function ...
        :param parameters_modes_processor_counts_dict_timing:
        :param parameters_modes_processor_counts_dict_memory:
        :return:
        """

        # Inform the user
        log.info("Setting missing data ...")

        # Set missing serial timing data
        if self.needs_timing: self.set_missing_serial_timing(parameters_modes_processor_counts_dict_timing)

        # Set missing serial memory data
        if self.needs_memory: self.set_missing_serial_memory(parameters_modes_processor_counts_dict_memory)

    # -----------------------------------------------------------------

    def set_missing_serial_timing(self, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        # Check if serial data is found for the different parameters set
        for parameter_set in parameters_modes_processor_counts_dict:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_timing: continue

            # Set missing serial timing through extrapolation of npackages
            if not self.has_serial_timing_for_parameters(parameter_set) and self.config.extrapolation.timing.npackages:
                self.set_missing_serial_timing_npackages(parameter_set, parameters_modes_processor_counts_dict)

            # Set missing serial timing through extrapolation of ncores
            if not self.has_serial_timing_for_parameters(parameter_set) and self.config.extrapolation.timing.ncores:
                self.set_missing_serial_timing_ncores(parameter_set, parameters_modes_processor_counts_dict)

            # Set missing serial timing through extrapolation of nwavelengths
            if not self.has_serial_timing_for_parameters(parameter_set) and self.config.extrapolation.timing.nwavelengths:
                self.set_missing_serial_timing_nwavelengths(parameter_set, parameters_modes_processor_counts_dict)

            # Set missing serial timing as timing that is not really serial ...
            if not self.has_serial_timing_for_parameters(parameter_set): self.set_missing_serial_timing_false(parameter_set)

    # -----------------------------------------------------------------

    def set_missing_serial_timing_npackages(self, parameter_set, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameter_set:
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        # Give warning
        log.warning("Serial (one core) timing data not found, extrapolating timing data created with other numbers of photon packages ...")

        the_other_npackages = None
        the_other_parameter_set = None

        # Loop over the other parameter sets
        other_parameter_sets = parameters_modes_processor_counts_dict.keys()
        other_parameter_sets.remove(parameter_set)
        if len(other_parameter_sets) == 0: return
        for other_parameter_set in other_parameter_sets:

            other_npackages = None

            # before parameter sets were tuples instead of dictionaries / Maps
            #for name in other_parameter_set:
            #    if name == "Packages": other_npackages = other_parameter_set["Packages"]
            #    elif other_parameter_set[name] != parameter_set[name]:
            #        break

            for index in range(len(self.different_parameters_timing)):

                name = self.different_parameters_timing[index]
                if name == "Packages":
                    if isinstance(other_parameter_set, tuple): other_npackages = other_parameter_set[index]
                    else: other_npackages = other_parameter_set
                elif other_parameter_set[index] != parameter_set[index]: break

            # If a break is not encountered
            else:

                # End the loop over the other parameter sets because we found a good one
                the_other_npackages = other_npackages
                the_other_parameter_set = other_parameter_set
                break

        # Only if a good match was found
        if the_other_parameter_set is None: return

        # Debugging
        log.debug("Found the parameter set: " + parameter_set_to_string_inline(the_other_parameter_set, self.different_parameters_timing))

        # Multiplication factor
        npackages_index = self.get_parameter_index_timing("Packages")
        if isinstance(parameter_set, tuple): npackages_factor = float(parameter_set[npackages_index]) / float(the_other_npackages)
        else: npackages_factor = float(parameter_set) / float(the_other_npackages) # the parameter set = only just the npackages value

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

        # GET TIMING DATA FROM PARAMETER SET THAT DIFFERS ONLY IN THE NUMBER OF PACKAGES
        # FILL IN FOR THE PARAMETER SET WHICH WE WERE SEARCHING SERIAL DATA FOR

        # Loop over the required phases
        for phase in self.simulation_phases_and_subphases_timing:

            # Emission phases
            if phase == "stellar" or phase == "dust":

                # Scale the runtime
                time = self.timing_data[phase][the_other_parameter_set][the_mode].times[the_serial_index] * npackages_factor

            # No communication for serial runs
            elif "communication" in phase: time = 0.0

            # Total simulation: not now, we first need the other times
            elif phase == "total": continue

            # Other phases
            else: time = self.timing_data[phase][the_other_parameter_set][the_mode].times[the_serial_index]

            # total
            self.serial_timing[parameter_set][phase].time = time
            self.serial_timing[parameter_set][phase].error = 0.0
            self.serial_timing_ncores[parameter_set][phase] = 1

            # Also use the extrapolated value in the runtimes
            if self.config.extrapolation.timing.in_times:
                for mode in self.timing_data[phase][parameter_set]:  # Loop over all modes but not hybrid modes!
                    if "hybrid" in mode and (phase in multiprocessing_phases or not self.config.extrapolation.timing.for_hybrid): continue
                    self.timing_data[phase][parameter_set][mode].processor_counts.append(1)
                    self.timing_data[phase][parameter_set][mode].times.append(time)
                    self.timing_data[phase][parameter_set][mode].errors.append(0.0)

        # Total simulation
        if "total" in self.simulation_phases_and_subphases_timing:

            setup_time = self.serial_timing[parameter_set]["setup"].time
            stellar_time = self.serial_timing[parameter_set]["stellar"].time
            spectra_time = self.serial_timing[parameter_set]["spectra"].time
            dust_time = self.serial_timing[parameter_set]["dust"].time
            writing_time = self.serial_timing[parameter_set]["writing"].time
            waiting_time = self.serial_timing[parameter_set]["waiting"].time
            communication_time = self.serial_timing[parameter_set]["communication"].time
            intermediate_time = self.serial_timing[parameter_set]["intermediate"].time

            # Calculate the total time
            total_time = setup_time + stellar_time + spectra_time + dust_time + writing_time + waiting_time + communication_time + intermediate_time

            # Set the time
            self.serial_timing[parameter_set]["total"].time = total_time
            self.serial_timing[parameter_set]["total"].error = 0.0
            self.serial_timing_ncores[parameter_set]["total"] = 1

    # -----------------------------------------------------------------

    def set_missing_serial_timing_ncores(self, parameter_set, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameter_set:
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_missing_serial_timing_nwavelengths(self, parameter_set, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameter_set:
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        pass

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
        for phase in self.timing_data:

            max_time = 0.0
            max_time_error = None
            max_time_ncores = float("+inf")

            # Loop over the modes
            for mode in self.timing_data[phase][parameter_set]:

                # Get index for lowest processor count
                index = np.argmin(self.timing_data[phase][parameter_set][mode].processor_counts)

                # Get properties
                ncores = self.timing_data[phase][parameter_set][mode].processor_counts[index]
                time = self.timing_data[phase][parameter_set][mode].times[index]
                time_error = self.timing_data[phase][parameter_set][mode].errors[index]

                # Check condition and adapt accordingly
                if ncores < max_time_ncores or (ncores == max_time_ncores and time > max_time):

                    max_time = time
                    max_time_error = time_error
                    max_time_ncores = ncores

            # Set the time and error
            self.serial_timing[parameter_set][phase].time = max_time
            self.serial_timing[parameter_set][phase].error = max_time_error
            self.serial_timing_ncores[parameter_set][phase] = max_time_ncores

    # -----------------------------------------------------------------

    def has_serial_memory_for_parameters(self, parameter_set):

        """
        This function ...
        :param parameter_set:
        :return:
        """

        return parameter_set in self.serial_memory

    # -----------------------------------------------------------------

    def set_missing_serial_memory(self, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        # Check if serial data is found for the different parameters set
        for parameter_set in parameters_modes_processor_counts_dict:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_memory: continue

            # Set missing serial memory usage through extrapolation of the number of dust cells
            if not self.has_serial_memory_for_parameters(parameter_set) and self.config.extrapolation.memory.ncells:
                self.set_missing_serial_memory_dustcells(parameter_set, parameters_modes_processor_counts_dict)

            # Set missing serial memory usage through extrapolation of the number of wavelengths
            if not self.has_serial_memory_for_parameters(parameter_set) and self.config.extrapolation.memory.nwavelengths:
                self.set_missing_serial_memory_nwavelengths(parameter_set, parameters_modes_processor_counts_dict)

            # Set missing serial memory usage through extrapolation of the number of processes
            if not self.has_serial_memory_for_parameters(parameter_set) and self.config.extrapolation.memory.nprocesses:
                self.set_missing_serial_memory_nprocesses(parameter_set, parameters_modes_processor_counts_dict)

            # Set missing serial memory usage as memory usage that is not really serial ...
            if not self.has_serial_memory_for_parameters(parameter_set): self.set_missing_serial_memory_false(parameter_set)

    # -----------------------------------------------------------------

    def set_missing_serial_memory_dustcells(self, parameter_set, parameters_modes_processor_counts_dict):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_missing_serial_memory_nwavelengths(self, parameter_set, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameter_set:
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_missing_serial_memory_nprocesses(self, parameter_set, parameters_modes_processor_counts_dict):

        """
        This function ...
        :param parameter_set:
        :param parameters_modes_processor_counts_dict:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_missing_serial_memory_false(self, parameter_set):

        """
        This function ...
        :param parameter_set:
        :return:
        """

        # Give warning
        log.warning("Serial (one core) memory data not found, searching for highest memory consumption per process with the least amount of cores for each phase (any parallelization mode)")

        # Loop over the phases
        for phase in self.memory_data:

            max_memory = 0.0
            max_memory_error = None
            max_memory_ncores = float("+inf")

            # Loop over the modes
            for mode in self.memory_data[phase][parameter_set]:

                # Get the index for the lowest number of cores
                index = np.argmin(self.memory_data[phase][parameter_set][mode].processor_counts)

                # Get properties
                ncores = self.memory_data[phase][parameter_set][mode].processor_counts[index]
                memory = self.memory_data[phase][parameter_set][mode].memory[index]
                error = self.memory_data[phase][parameter_set][mode].errors[index]

                # Check condition and adapt accordingly
                if ncores < max_memory_ncores or (ncores == max_memory_ncores and memory > max_memory):

                    max_memory = memory
                    max_memory_error = error
                    max_memory_ncores = ncores

            # Set the memory usage and error
            self.serial_memory[parameter_set][phase].memory = max_memory
            self.serial_memory[parameter_set][phase].error = max_memory_error
            self.serial_memory_ncores[parameter_set][phase] = max_memory_ncores

    # -----------------------------------------------------------------

    def check_coverage_modes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the data coverage for the different parallelization modes ...")

        # Timing
        if self.needs_timing: self.check_coverage_modes_timing()

        # Memory
        if self.needs_memory: self.check_coverage_modes_memory()

    # -----------------------------------------------------------------

    def check_coverage_modes_timing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the timing data coverage for the different parallelization modes ...")

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets_timing:

            # Loop over the modes
            for mode in self.modes_for_timing_parameter_set(parameter_set):

                # Loop over the phases
                for phase in self.timing_data:

                    # Check
                    if len(self.timing_data[phase][parameter_set][mode].processor_counts) < 2:

                        # Debugging
                        if self.has_multiparameters_timing: log.debug("Removing timing data for " + phase_names[phase] + " for " + parameter_set_to_string_inline(parameter_set, self.different_parameters_timing) + " in " + mode + " mode ...")
                        else: log.debug("Removing timing data for " + phase_names[phase] + " in " + mode + " mode ...")

                        # Remove
                        del self.timing_data[phase][parameter_set][mode]

                        # If no modes left for the parameter set
                        if len(self.timing_data[phase][parameter_set]) == 0: del self.timing_data[phase][parameter_set]

    # -----------------------------------------------------------------

    def check_coverage_modes_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the memory data coverage for the different parallelization modes ...")

        #print(self.parameter_sets_memory)
        parameter_set = self.parameter_sets_memory[0]
        #print(self.modes_for_memory_parameter_set(parameter_set))

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets_memory:

            # Loop over the modes
            for mode in self.modes_for_memory_parameter_set(parameter_set):

                # Loop over the phases
                for phase in self.memory_data:

                    #print(self.memory_data[phase][parameter_set][mode].processor_counts)

                    # Check
                    if len(self.memory_data[phase][parameter_set][mode].processor_counts) < 2:

                        # Debugging
                        if self.has_multiparameters_memory: log.debug("Removing memory data for " + phase_names[phase] + " for " + parameter_set_to_string_inline(parameter_set, self.different_parameters_memory) + " in " + mode + " mode ...")
                        else: log.debug("Removing memory data for " + phase_names[phase] + " in " + mode + " mode ...")

                        # Remove
                        del self.memory_data[phase][parameter_set][mode]

                        # If no modes left for the parameter set
                        if len(self.memory_data[phase][parameter_set]) == 0: del self.memory_data[phase][parameter_set]

    # -----------------------------------------------------------------

    def add_fixed_multiprocessing_phases_runtimes(self):

        """
        This function ...
        :return:
        """

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets_timing:

            # Loop over the modes
            for mode in self.modes_for_timing_parameter_set(parameter_set):

                # Skip multithreading
                if mode == "multithreading": continue

                # Get threads per process
                if "hybrid" in mode: threadspp = nthreads_per_process_for_hybrid_mode(mode)
                else: threadspp = 1

                # For one process, ncores = threadspp (x 1)
                ncores = threadspp

                # Check whether 1 process is present in MPI phases
                for phase in multiprocessing_phases:

                    if ncores not in self.timing_data[phase][parameter_set][mode].processor_counts:

                        self.timing_data[phase][parameter_set][mode].processor_counts.append(ncores)
                        self.timing_data[phase][parameter_set][mode].times.append(0.0)
                        self.timing_data[phase][parameter_set][mode].errors.append(0.0)

    # -----------------------------------------------------------------

    @property
    def parameter_sets_timing(self):

        """
        This function ...
        :return:
        """

        #phase = "total" # not always present
        phase = self.timing_data.keys()[0]

        return self.timing_data[phase].keys()

    # -----------------------------------------------------------------

    @property
    def parameter_sets_memory(self):

        """
        This function ...
        :return:
        """

        # phase = "total" # not always present
        phase = self.memory_data.keys()[0]

        return self.memory_data[phase].keys()

    # -----------------------------------------------------------------

    def modes_for_timing_parameter_set(self, parameter_set):

        """
        This function ...
        :param parameter_set:
        :return:
        """

        # phase = "total" # not always present
        phase = self.timing_data.keys()[0]

        return self.timing_data[phase][parameter_set].keys()

    # -----------------------------------------------------------------

    def modes_for_memory_parameter_set(self, parameter_set):

        """
        This function ...
        :param parameter_set:
        :return:
        """

        # phase = "total" # not always present
        phase = self.memory_data.keys()[0]

        #for phase in sel.t

        return self.memory_data[phase][parameter_set].keys()

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
        for parameter_set in self.parameter_sets_timing:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_timing: continue

            # Loop over the phases
            for phase in self.simulation_phases_and_subphases_timing:

                # Don't fit for this phase
                if phase in dont_fit_phases: continue

                # Debugging
                log.debug("Fitting timing data for " + phase_names[phase] + " ...")

                # Create a dictionary that stores the fitted parameters for each different mode
                parameters = dict()
                parameter_errors = dict()

                # Generate the fit function
                behaviour = generate_scaling_behaviour(phase, self.config.fitting)
                fit_function, nparameters = generate_fit_function(behaviour)

                # Debugging
                log.debug("Expected scaling behaviour: " + behaviour_to_function_string(behaviour))

                # Loop over the different parallelization modes (the different curves)
                for mode in self.timing_data[phase][parameter_set]:

                    # Get the list of processor counts, runtimes and errors
                    processor_counts = self.timing_data[phase][parameter_set][mode].processor_counts
                    times = self.timing_data[phase][parameter_set][mode].times
                    errors = self.timing_data[phase][parameter_set][mode].errors

                    # Sort the lists
                    processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                    # Set the weights of the different timing points for the fitting procedure
                    weights = errors if not np.any(np.isinf(errors)) else None
                    if np.count_nonzero(errors) == 0: weights = None

                    # Calculate the normalized processor counts (relative to the number of processors used for the serial run)
                    normalized_processor_counts = processor_counts / self.serial_timing_ncores[parameter_set][phase]

                    # Hack so that number of parameters doesn't need to be inferred by Scipy from the fit function prescription (because it uses *args)
                    from scipy.optimize.minpack import _initialize_feasible, prepare_bounds
                    n = nparameters
                    bounds = (-np.inf, np.inf)
                    lb, ub = prepare_bounds(bounds, n)
                    p0 = _initialize_feasible(lb, ub)

                    # Fit parameters for the speedups to Amdahl's law
                    #print(mode, p0, normalized_processor_counts, times, weights)
                    popt, pcov = curve_fit(fit_function, normalized_processor_counts, times, sigma=weights, absolute_sigma=False, p0=p0)
                    perr = np.sqrt(np.diag(pcov))

                    parameters[mode] = popt
                    parameter_errors[mode] = perr

                # Add the parameters
                self.timing_fit_functions[phase][parameter_set] = fit_function
                self.timing_fit_parameters[phase][parameter_set] = parameters
                self.timing_fit_parameter_errors[phase][parameter_set] = parameter_errors

    # -----------------------------------------------------------------

    def fit_timing_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting the timing data ...")

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets_timing:

            # Skip parameter set
            if parameter_set in self.ignore_parameter_sets_timing: continue

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
                serial = Measurement(serial_time, serial_error)

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
                        time = Measurement(times[i], errors[i])

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

        # Loop over the parameter sets
        for parameter_set in self.parameter_sets_memory:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_memory: continue

            # Loop over the phases
            for phase in self.simulation_phases_and_subphases_memory:

                # Get the serial (1 process) memory consumption (and error) for this phase (create a Quantity object)
                #serial_memory = self.serial_memory[phase].memory
                #serial_error = self.serial_memory[phase].error
                #serial = Quantity(serial_memory, serial_error)

                # Create a dictionary that stores the fitted parameters for each different mode
                parameters = dict()
                parameter_errors = dict()

                # Loop over the different parallelization modes (the different curves)
                for mode in self.memory_data[phase][parameter_set]:

                    # Skip modes that are not data parallel
                    if not "task+data" in mode: continue

                    # Get the list of processor counts, memory and errors
                    processor_counts = self.memory_data[phase][parameter_set][mode].processor_counts
                    memories = self.memory_data[phase][parameter_set][mode].memory
                    errors = self.memory_data[phase][parameter_set][mode].errors

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
                    if len(processor_counts) < 6:

                        popt, pcov = curve_fit(memory_scaling, nprocesses, memories, sigma=weights, absolute_sigma=False)
                        perr = np.sqrt(np.diag(pcov))
                        parameters[mode] = Map({"a": popt[0], "a_error": perr[0], "b": popt[1], "b_error": perr[1], "c": 0.0, "c_error": 0.0})

                    else:

                        popt, pcov = curve_fit(modified_memory_scaling, nprocesses, memories, sigma=weights, absolute_sigma=False)
                        perr = np.sqrt(np.diag(pcov))
                        parameters[mode] = Map({"a": popt[0], "a_error": perr[0], "b": popt[1], "b_error": perr[1], "c": popt[2], "c_error": perr[2]})

                # Add the parameters
                self.memory_fit_parameters[phase][parameter_set] = parameters

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
        log.info("Plotting the runtimes for each phase ...")

        # Loop over the phases
        for phase in self.simulation_phases_timing: self.plot_runtimes_phase(phase)

    # -----------------------------------------------------------------

    def plot_speedups(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the speedups for each phase ...")

        # Loop over the phases
        for phase in self.parallel_simulation_phases_timing: self.plot_speedups_phase(phase)

    # -----------------------------------------------------------------

    def plot_efficiencies(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the efficiencies ...")

        # Loop over the phases
        for phase in self.parallel_simulation_phases_timing: self.plot_efficiencies_phase(phase)

    # -----------------------------------------------------------------

    def plot_cpu_times(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the CPU times ...")

        # Loop over the phases
        for phase in self.simulation_phases_timing: self.plot_cpu_times_phase(phase)

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the memory usage for each phase ...")

        # Loop over the phases to be plotted
        for phase in self.simulation_phases_memory: self.plot_memory_phase(phase)

    # -----------------------------------------------------------------

    def plot_memory_gain(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the memory gain for each phase ...")

        # Loop over the phases
        for phase in self.simulation_phases_memory: self.plot_memory_gain_phase(phase)

    # -----------------------------------------------------------------

    def plot_total_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the total memory usage for each phase ...")

        # Loop over the phases
        for phase in self.simulation_phases_memory: self.plot_total_memory_phase(phase)

    # -----------------------------------------------------------------

    @property
    def figsize(self):

        """
        This function ...
        :return:
        """

        return (self.config.plot.xsize, self.config.plot.ysize)

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
        plt.figure(figsize=self.figsize)
        plt.clf()

        # Set background
        set_background(plt.gcf(), plt.gca(), self.config.plot)

        # Create a set that stores the tick labels for the plot
        ticks = set()
        ticks_ncores = set()

        # Phases for this plot
        if phase == "communication" and self.config.split_communication: plot_phases = communication_phases
        else: plot_phases = [phase]
        multi_phases = len(plot_phases) > 1

        # Check whether we have multiple parameter sets
        multi_parameter_sets = len(self.timing_data[plot_phases[0]]) > 1

        # The plot handles for the different parameter sets (or phases if there is only one parameter set)
        handles = defaultdict(list)

        # Set unit for runtimes
        if self.config.normalize_runtimes: ylabel = "Normalized runtime"
        else: ylabel = "Runtime (s)"

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Loop over the different parameter sets (different ski files)
            for parameter_set in self.timing_data[plot_phase]:

                # Ignore
                if parameter_set in self.ignore_parameter_sets_timing: continue

                # Debugging
                log.debug("Plotting the runtimes for the " + phase_names[plot_phase] + " for the simulations with parameters: " + parameter_set_to_string_inline(parameter_set, self.different_parameters_timing))

                # Get the serial runtime (and error) for this phase (create a Quantity object)
                if self.config.normalize_runtimes:
                    serial_time = self.serial_timing[parameter_set][plot_phase].time
                    serial_error = self.serial_timing[parameter_set][plot_phase].error
                    serial = Measurement(serial_time, serial_error)
                else: serial = None

                # Loop over the different parallelization modes (the different curves)
                for mode in sorted(self.timing_data[plot_phase][parameter_set].keys()): # sort alphabetically

                    # Get the list of processor counts, runtimes and errors
                    processor_counts = self.timing_data[plot_phase][parameter_set][mode].processor_counts
                    times = self.timing_data[plot_phase][parameter_set][mode].times
                    errors = self.timing_data[plot_phase][parameter_set][mode].errors

                    # Sort the lists
                    processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                    # Get the array of x values
                    if self.config.hybridisation: x_values = processor_counts
                    # actually, in hybridisation mode,
                    # processor_counts contains actually the process counts
                    else:
                        x_values = get_x_values(processor_counts, mode, self.config.x_quantity)
                        if x_values is None: continue

                    # Normalize runtimes
                    if self.config.normalize_runtimes:

                        new_times = []
                        new_errors = []

                        for i in range(len(processor_counts)):

                            # Create a quantity for the current runtime
                            time = Measurement(times[i], errors[i])

                            # Calculate the normalized runtim
                            speedup = time / serial

                            # Adjust the speedup range
                            #speedup_range.adjust(speedup.value)

                            # Add the value and the propagated error of the normalized runtime
                            new_times.append(speedup.value)
                            new_errors.append(speedup.error)

                        # Convert to arrays
                        times = np.array(new_times)
                        errors = np.array(new_errors)

                    # Determine the label
                    label = mode

                    # Add phase to label
                    if multi_parameter_sets and multi_phases: label += " (" + plot_phase + ")"

                    # Plot the data points for this mode
                    fmt = '' if self.config.plot.connect_points else 'o'
                    handle = plt.errorbar(x_values, times, errors, marker='.', label=label,
                                          linewidth=self.config.plot.linewidth, fmt=fmt, markersize=self.config.plot.markersize)

                    # Determine key for handle dictionary
                    if not multi_parameter_sets and multi_phases: key = plot_phase
                    else: key = parameter_set

                    # Add the handle to the dictionary
                    handles[key].append(handle)

                    # Add the appropriate ticks
                    ticks |= set(x_values)
                    ticks_ncores |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads) and the y axis (time)
        if self.config.xlog: plt.xscale('log')
        if self.config.ylog and phase not in phases_not_logaritmic_runtimes: plt.yscale('log')

        # Check if there less than 2 data points
        if len(ticks) < 2:
            log.warning("Not enough data for plotting")
            return

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks_ncores = sorted(ticks_ncores)
        #ticks.append(ticks[-1] * 2)

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Plot a line that denotes linear scaling (runtime = serial runtime / ncores), for parallel phases (including the total simulation)
            if not self.config.hybridisation and plot_phase in parallel_phases:

                # Loop over the parameter sets
                for parameter_set in self.timing_data[plot_phase]:

                    # Ignore
                    if parameter_set in self.ignore_parameter_sets_timing: continue

                    # Get the serial runtime (and error) for this phase (create a Quantity object)
                    serial_time = self.serial_timing[parameter_set][plot_phase].time
                    serial_error = self.serial_timing[parameter_set][plot_phase].error
                    serial = Measurement(serial_time, serial_error)

                    # Calculate the ideal runtimes
                    runtimes = [serial.value / ncores for ncores in ticks_ncores]

                    # Plot the line
                    plt.plot(ticks, runtimes, linestyle='--', linewidth=self.config.plot.linewidth)

            # Plot the fit
            if not self.config.hybridisation and self.config.fit and self.config.fitting.plot_fit:

                # Loop over the parameter sets
                for parameter_set in self.timing_fit_functions[plot_phase]:

                    # Get the fit function
                    fit_function = self.timing_fit_functions[plot_phase][parameter_set]

                    # Get the number of processers taken as the reference for normalization, and thus calculation of the speedups and as reference for the fit
                    #reference_ncores = self.serial_timing_ncores[phase]

                    fit_x_values = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                    fit_ncores = np.logspace(np.log10(ticks_ncores[0]), np.log10(ticks_ncores[-1]), 50)
                    for mode in self.timing_fit_parameters[plot_phase][parameter_set]:

                        # Get the parameter values
                        parameters = self.timing_fit_parameters[plot_phase][parameter_set][mode]

                        # Get the parameter errors
                        parameter_errors = self.timing_fit_parameter_errors[plot_phase][parameter_set][mode]

                        # Calculate the fitted times
                        fit_times = [fit_function(n, *parameters) for n in fit_ncores]

                        # Add the plot
                        plt.plot(fit_x_values, fit_times, color="grey", linewidth=self.config.plot.linewidth)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.yaxis.set_major_formatter(LogFormatter())
        #ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))

        plt.xlim(ticks[0], ticks[-1])

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set borders
        set_borders(ax, self.config.plot)

        # Add axis labels
        if self.config.hybridisation: plt.xlabel("Number of processes", fontsize=self.config.plot.label_fontsize)
        else: plt.xlabel("Number of " + self.config.x_quantity, fontsize=self.config.plot.label_fontsize)
        plt.ylabel(ylabel, fontsize=self.config.plot.label_fontsize)

        # Add the legends
        # if self.config.hybridisation: plt.legend(title="Number of cores")
        # else: plt.legend(title="Parallelization modes")
        pa_or_ph = "phases" if not multi_parameter_sets and multi_phases else "parameters"
        add_legends(ax, handles, self.different_parameters_timing, self.config.plot, "runtime", pa_or_ph)

        # Set the plot title
        if self.config.plot.add_titles:
            title = "Scaling of the " + phase_labels_timing[phase].lower()
            title = split_in_lines(title)
            plt.suptitle(title, fontsize=self.config.plot.title_fontsize)

        # Set file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "runtimes_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path, transparent=True)
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
        plt.figure(figsize=self.figsize)
        plt.clf()

        # Set background
        set_background(plt.gcf(), plt.gca(), self.config.plot)

        # Create a set that stores the tick labels for the plot
        ticks = set()
        ticks_ncores = set()

        # Keep track of the minimal and maximal speedup
        speedup_range = RealRange.zero()

        # Phases for this plot
        if phase == "communication" and self.config.split_communication: plot_phases = communication_phases
        else: plot_phases = [phase]
        multi_phases = len(plot_phases) > 1

        # Check whether we have multiple parameter sets
        multi_parameter_sets = len(self.timing_data[plot_phases[0]]) > 1

        # The plot handles for the different parameter sets (or phases if there is only one parameter set)
        handles = defaultdict(list)

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Loop over the different parameter sets (different ski files)
            for parameter_set in self.timing_data[plot_phase]:

                # Ignore
                if parameter_set in self.ignore_parameter_sets_timing: continue

                # Debugging
                log.debug("Plotting the speedups for the " + phase_names[plot_phase] + " for the simulations with parameters: " + parameter_set_to_string_inline(parameter_set, self.different_parameters_timing))

                # Get the serial runtime (and error) for this phase (create a Quantity object)
                serial_time = self.serial_timing[parameter_set][plot_phase].time
                serial_error = self.serial_timing[parameter_set][plot_phase].error
                serial = Measurement(serial_time, serial_error)

                # Loop over the different parallelization modes (the different curves)
                for mode in sorted(self.timing_data[plot_phase][parameter_set].keys()): # sort alphabetically

                    # Get the list of processor counts, runtimes and errors
                    processor_counts = self.timing_data[plot_phase][parameter_set][mode].processor_counts
                    times = self.timing_data[plot_phase][parameter_set][mode].times
                    errors = self.timing_data[plot_phase][parameter_set][mode].errors

                    # Sort the lists
                    processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                    # Get the array of x values
                    if self.config.hybridisation: x_values = processor_counts
                    # actually, in hybridisation mode,
                    # processor_counts contains actually the process counts
                    else:
                        x_values = get_x_values(processor_counts, mode, self.config.x_quantity)
                        if x_values is None: continue

                    # Calculate the speedups and the errors on the speedups
                    speedups = []
                    speedup_errors = []
                    for i in range(len(processor_counts)):

                        # Create a quantity for the current runtime
                        time = Measurement(times[i], errors[i])

                        # Calculate the speedup based on the current runtime and the serial runtime
                        speedup = serial / time

                        # Adjust the speedup range
                        speedup_range.adjust(speedup.value)

                        # Add the value and the propagated error of the speedup to the appropriate lists
                        speedups.append(speedup.value)
                        speedup_errors.append(speedup.error)

                    # Convert to arrays
                    speedups = np.array(speedups)
                    speedup_errors = np.array(speedup_errors)

                    # Determine the label
                    label = mode

                    # Add phase to label
                    if multi_parameter_sets and multi_phases: label += " (" + plot_phase + ")"

                    # Eliminate Nans
                    not_nan = np.logical_not(np.isnan(speedups))
                    processor_counts = processor_counts[not_nan]
                    speedups = speedups[not_nan]
                    speedup_errors = speedup_errors[not_nan]

                    # If there are NaNs in the speedup errors
                    if np.any(np.isnan(speedup_errors)): speedup_errors = np.zeros_like(speedup_errors)

                    # If all speedup points are zero for some reason, skip
                    if np.all(speedups == 0): continue

                    # Plot the data points for this curve
                    fmt = '' if self.config.plot.connect_points else 'o'
                    handle = plt.errorbar(x_values, speedups, speedup_errors, marker='.', label=label,
                                          linewidth=self.config.plot.linewidth, fmt=fmt, markersize=self.config.plot.markersize)

                    # Determine key for handle dictionary
                    if not multi_parameter_sets and multi_phases: key = plot_phase
                    else: key = parameter_set

                    # Add the plot handle
                    handles[key].append(handle)

                    # Add the appropriate ticks
                    ticks |= set(x_values)
                    ticks_ncores |= set(processor_counts)

        # Use a logarithmic scale for both axes
        if self.config.xlog: plt.xscale('log')
        if self.config.ylog: plt.yscale('log')

        # No data points
        if len(ticks) == 0:
            log.warning("Could not make a plot for the speedups of the " + phase_names[phase] + " because there was no valid data")
            return

        # Check if there less than 2 data points
        if len(ticks) < 2:
            log.warning("Not enough data for plotting")
            return

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks_ncores = sorted(ticks_ncores)
        #ticks.append(ticks[-1] * 2)

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Plot the fit
            if not self.config.hybridisation and self.config.fit and self.config.fitting.plot_fit:

                # Loop over the parameter sets
                for parameter_set in self.timing_fit_functions[plot_phase]:

                    # Get the fit function
                    fit_function = self.timing_fit_functions[plot_phase][parameter_set]

                    fit_x_values = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                    fit_ncores = np.logspace(np.log10(ticks_ncores[0]), np.log10(ticks_ncores[-1]), 50)
                    for mode in self.timing_fit_parameters[plot_phase][parameter_set]:

                        # Get the parameter values
                        parameters = self.timing_fit_parameters[plot_phase][parameter_set][mode]

                        # Get the parameter errors
                        parameter_errors = self.timing_fit_parameter_errors[plot_phase][parameter_set][mode]

                        # Calculate the fitted speedups
                        fit_speedups = [serial.value / fit_function(n, *parameters) for n in fit_ncores]

                        # Add the plot
                        plt.plot(fit_x_values, fit_speedups, color="grey", linewidth=self.config.plot.linewidth)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        if not self.config.hybridisation: ax.set_yticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlim(ticks[0], ticks[-1])
        #if not self.config.hybridisation: plt.ylim(ticks[0], ticks[-1])
        plt.ylim(speedup_range.min, speedup_range.max)

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set borders
        set_borders(ax, self.config.plot)

        # Add axis labels
        if self.config.hybridisation: plt.xlabel("Number of processes", fontsize=self.config.plot.label_fontsize)
        else: plt.xlabel("Number of " + self.config.x_quantity, fontsize=self.config.plot.label_fontsize)
        plt.ylabel("Speedup", fontsize=self.config.plot.label_fontsize)

        # Plot a line that denotes linear scaling (speedup = nthreads)
        if not self.config.hybridisation and phase in parallel_phases: plt.plot(ticks, ticks, linestyle='--', linewidth=self.config.plot.linewidth)

        # Add the legends
        pa_or_ph = "phases" if not multi_parameter_sets and multi_phases else "parameters"
        add_legends(ax, handles, self.different_parameters_timing, self.config.plot, "speedup", pa_or_ph)

        # Set the plot title
        if self.config.plot.add_titles:
            title = "Speedup of the " + phase_labels_timing[phase].lower()
            title = split_in_lines(title)
            plt.suptitle(title, fontsize=self.config.plot.title_fontsize)

        # Set file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "speedups_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path, transparent=True)
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
        plt.figure(figsize=self.figsize)
        plt.clf()

        # Set background
        set_background(plt.gcf(), plt.gca(), self.config.plot)

        # Create a set that stores the tick labels for the plot
        ticks = set()
        ticks_ncores = set()

        # Phases for this plot
        if phase == "communication" and self.config.split_communication: plot_phases = communication_phases
        else: plot_phases = [phase]
        multi_phases = len(plot_phases) > 1

        # Check whether we have multiple parameter sets
        multi_parameter_sets = len(self.timing_data[plot_phases[0]]) > 1

        # The plot handles for the different parameter sets (or phases if there is only one parameter set)
        handles = defaultdict(list)

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Loop over the different parameter sets (different ski files)
            for parameter_set in self.timing_data[plot_phase]:

                # Ignore
                if parameter_set in self.ignore_parameter_sets_timing: continue

                # Debugging
                log.debug("Plotting the efficiencies for the " + phase_names[plot_phase] + " for the simulations with parameters: " + parameter_set_to_string_inline(parameter_set, self.different_parameters_timing))

                # Get the serial runtime (and error) for this phase (create a Quantity object)
                serial_time = self.serial_timing[parameter_set][plot_phase].time
                serial_error = self.serial_timing[parameter_set][plot_phase].error
                serial = Measurement(serial_time, serial_error)
                serial_ncores = self.serial_timing_ncores[parameter_set][plot_phase] if plot_phase in self.serial_timing_ncores[parameter_set] else 1

                # Loop over the different parallelization modes (the different curves)
                for mode in sorted(self.timing_data[plot_phase][parameter_set].keys()): # sort alphabetically

                    # Get the list of processor counts, runtimes and errors
                    processor_counts = self.timing_data[plot_phase][parameter_set][mode].processor_counts
                    times = self.timing_data[plot_phase][parameter_set][mode].times
                    errors = self.timing_data[plot_phase][parameter_set][mode].errors

                    # Sort the lists
                    processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                    # Get array of number of used cores
                    if self.config.hybridisation: ncores = np.ones(len(processor_counts)) * int(mode.split(" cores")[0])
                    else: ncores = processor_counts

                    # Get the array of x values
                    if self.config.hybridisation: x_values = processor_counts
                    # actually, in hybridisation mode,
                    # processor_counts contains actually the process counts
                    else:
                        x_values = get_x_values(processor_counts, mode, self.config.x_quantity)
                        if x_values is None: continue

                    # Calculate the efficiencies and the errors on the efficiencies
                    efficiencies = []
                    efficiency_errors = []
                    for i in range(len(processor_counts)):

                        # Create a quantity for the current runtime
                        time = Measurement(times[i], errors[i])

                        # Calculate the efficiency based on the current runtime and the serial runtime
                        speedup = serial / time
                        efficiency = speedup.value / (ncores[i]/serial_ncores)
                        efficiency_error = speedup.error / (ncores[i]/serial_ncores)

                        # Add the value and the propagated error of the efficiency to the appropriate lists
                        efficiencies.append(efficiency)
                        efficiency_errors.append(efficiency_error)

                    # Determine the label
                    label = mode

                    # Add phase to label
                    if multi_parameter_sets and multi_phases: label += " (" + plot_phase + ")"

                    # Plot the data points for this curve
                    fmt = '' if self.config.plot.connect_points else 'o'
                    handle = plt.errorbar(x_values, efficiencies, efficiency_errors, marker='.', label=label,
                                          linewidth=self.config.plot.linewidth, fmt=fmt, markersize=self.config.plot.markersize)

                    # Determine key for handle dictionary
                    if not multi_parameter_sets and multi_phases: key = plot_phase
                    else: key = parameter_set

                    # Add the handle
                    handles[key].append(handle)

                    # Add the appropriate ticks
                    ticks |= set(x_values)
                    ticks_ncores |= set(processor_counts)

        # Use a logaritmic scale for the x axis (nthreads)
        if self.config.xlog: plt.xscale('log')
        if self.config.ylog: log.warning("Not using y log scale for efficieny plots")

        # Check if there less than 2 data points
        if len(ticks) < 2:
            log.warning("Not enough data for plotting")
            return

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks_ncores = sorted(ticks_ncores)
        #ticks.append(ticks[-1] * 2)

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Plot fit
            if not self.config.hybridisation and self.config.fit and self.config.fitting.plot_fit:

                # Loop over the parameter sets
                for parameter_set in self.timing_fit_functions[plot_phase]:

                    # Get the fit function
                    fit_function = self.timing_fit_functions[plot_phase][parameter_set]

                    fit_x_values = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                    fit_ncores = np.logspace(np.log10(ticks_ncores[0]), np.log10(ticks_ncores[-1]), 50)
                    for mode in self.timing_fit_parameters[plot_phase][parameter_set]:

                        # Get the parameter values
                        parameters = self.timing_fit_parameters[plot_phase][parameter_set][mode]

                        # Get the parameter errors
                        parameter_errors = self.timing_fit_parameter_errors[plot_phase][parameter_set][mode]

                        # Calculate the fitted efficiencies
                        fit_efficiencies = [serial.value / fit_function(n, *parameters) / n for n in fit_ncores]

                        # Add the plot
                        plt.plot(fit_x_values, fit_efficiencies, color="grey", linewidth=self.config.plot.linewidth)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=False))
        #ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlim(ticks[0], ticks[-1])
        #plt.ylim(0, 1.1)

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set borders
        set_borders(ax, self.config.plot)

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes", fontsize=self.config.plot.label_fontsize)
        else: plt.xlabel("Number of " + self.config.x_quantity, fontsize=self.config.plot.label_fontsize)
        plt.ylabel("Efficiency", fontsize=self.config.plot.label_fontsize)

        # Add the legends
        # if self.config.hybridisation: plt.legend(title="Number of cores")
        # else: plt.legend(title="Parallelization modes")
        pa_or_ph = "phases" if not multi_parameter_sets and multi_phases else "parameters"
        add_legends(ax, handles, self.different_parameters_timing, self.config.plot, "efficiency", pa_or_ph)

        # Set the plot title
        if self.config.plot.add_titles:
            title = "Efficiency of the " + phase_labels_timing[phase].lower()
            title = split_in_lines(title)
            plt.suptitle(title, fontsize=self.config.plot.title_fontsize)

        # Determine the path
        if self.config.output is not None: file_path = fs.join(self.config.output, "efficiencies_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path, transparent=True)
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
        plt.figure(figsize=self.figsize)
        plt.clf()

        # Set background
        set_background(plt.gcf(), plt.gca(), self.config.plot)

        # Create a set that stores the tick labels for the plot
        ticks = set()
        ticks_ncores = set()

        # Phases for this plot
        if phase == "communication" and self.config.split_communication: plot_phases = communication_phases
        else: plot_phases = [phase]
        multi_phases = len(plot_phases) > 1

        # Check whether we have multiple parameter sets
        multi_parameter_sets = len(self.timing_data[plot_phases[0]]) > 1

        # The plot handles for the different parameter sets (or phases if there is only one parameter set)
        handles = defaultdict(list)

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Loop over the different parameter sets (different ski files)
            for parameter_set in self.timing_data[plot_phase]:

                # Ignore
                if parameter_set in self.ignore_parameter_sets_timing: continue

                # Debugging
                log.debug("Plotting the CPU times for the " + phase_names[plot_phase] + " for the simulations with parameters: " + parameter_set_to_string_inline(parameter_set, self.different_parameters_timing))

                # Loop over the different parallelization modes (the different curves)
                for mode in sorted(self.timing_data[plot_phase][parameter_set].keys()): # sort alphabetically

                    # Get the list of processor counts, runtimes and errors
                    processor_counts = self.timing_data[plot_phase][parameter_set][mode].processor_counts
                    times = self.timing_data[plot_phase][parameter_set][mode].times
                    errors = self.timing_data[plot_phase][parameter_set][mode].errors

                    # Sort the lists
                    processor_counts, times, errors = sort_lists(processor_counts, times, errors, to_arrays=True)

                    # Get array of number of used cores
                    if self.config.hybridisation: ncores = np.ones(len(processor_counts)) * int(mode.split(" cores")[0])
                    else: ncores = processor_counts

                    # Get the array of x values
                    if self.config.hybridisation: x_values = processor_counts
                    # actually, in hybridisation mode,
                    # processor_counts contains actually the process counts
                    else:
                        x_values = get_x_values(processor_counts, mode, self.config.x_quantity)
                        if x_values is None: continue

                    # Get list of process count
                    #processes = nprocesses_from_mode(mode, processor_counts)

                    # Multiply to get total
                    times *= ncores
                    errors *= ncores

                    # Determine the label
                    label = mode

                    # Add phase to label
                    if multi_parameter_sets and multi_phases: label += " (" + plot_phase + ")"

                    # Plot the data points for this mode
                    fmt = '' if self.config.plot.connect_points else 'o'
                    handle = plt.errorbar(x_values, times, errors, marker='.', label=label,
                                          fmt=fmt, linewidth=self.config.plot.linewidth, markersize=self.config.plot.markersize)

                    # Determine key for handle dictionary
                    if not multi_parameter_sets and multi_phases: key = plot_phase
                    else: key = parameter_set

                    # Add the plot handle
                    handles[key].append(handle)

                    # Add the appropriate ticks
                    ticks |= set(x_values)
                    ticks_ncores |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        if self.config.xlog: plt.xscale("log")
        if self.config.ylog and phase not in phases_not_logaritmic_runtimes: plt.yscale("log")

        # Check if there less than 2 data points
        if len(ticks) < 2:
            log.warning("Not enough data for plotting")
            return

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks_ncores = sorted(ticks_ncores)
        #ticks.append(ticks[-1] * 2)

        # Loop over the plot phases
        for plot_phase in plot_phases:

            # Plot a line that denotes linear scaling (CPU time = serial runtime / ncores)
            if not self.config.hybridisation and plot_phase in parallel_phases:

                # Loop over the parameter sets
                for parameter_set in self.timing_data[plot_phase]:

                    # Ignore
                    if parameter_set in self.ignore_parameter_sets_timing: continue

                    # Get the serial runtime (and error) for this phase (create a Quantity object)
                    serial_time = self.serial_timing[parameter_set][plot_phase].time
                    serial_error = self.serial_timing[parameter_set][plot_phase].error
                    serial = Measurement(serial_time, serial_error)

                    # Calculate the ideal runtimes
                    runtimes = [serial.value] * len(ticks_ncores)

                    # Plot the line
                    plt.plot(ticks, runtimes, linestyle='--', linewidth=self.config.plot.linewidth)

            # Plot the fit
            if not self.config.hybridisation and self.config.fit and self.config.fitting.plot_fit:

                # Loop over the parameter sets
                for parameter_set in self.timing_fit_functions[plot_phase]:

                    # Get the fit function
                    fit_function = self.timing_fit_functions[plot_phase][parameter_set]

                    # Get the number of processers taken as the reference for normalization, and thus calculation of the speedups and as reference for the fit
                    # reference_ncores = self.serial_timing_ncores[phase]

                    fit_x_values = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                    fit_ncores = np.logspace(np.log10(ticks_ncores[0]), np.log10(ticks_ncores[-1]), 50)
                    for mode in self.timing_fit_parameters[plot_phase][parameter_set]:

                        # Get the parameter values
                        parameters = self.timing_fit_parameters[plot_phase][parameter_set][mode]

                        # Get the parameter errors
                        parameter_errors = self.timing_fit_parameter_errors[plot_phase][parameter_set][mode]

                        # Calculate the fitted times
                        fit_times = [fit_function(n, *parameters) * n for n in fit_ncores]

                        # Add the plot
                        plt.plot(fit_x_values, fit_times, color="grey", linewidth=self.config.plot.linewidth)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        #ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlim(ticks[0], ticks[-1])

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set borders
        set_borders(ax, self.config.plot)

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes", fontsize=self.config.plot.label_fontsize)
        else: plt.xlabel("Number of " + self.config.x_quantity, fontsize=self.config.plot.label_fontsize)
        plt.ylabel("CPU time (s)", fontsize=self.config.plot.label_fontsize)

        # Add the legends
        # if self.config.hybridisation: plt.legend(title="Number of cores")
        # else: plt.legend(title="Parallelization modes")
        pa_or_ph = "phases" if not multi_parameter_sets and multi_phases else "parameters"
        add_legends(ax, handles, self.different_parameters_timing, self.config.plot, "CPU-time", pa_or_ph)

        # Set the plot title
        if self.config.add_titles:
            title = "Scaling of the total CPU time of the " + phase_labels_timing[phase].lower()
            title = split_in_lines(title)
            plt.title(title, fontsize=self.config.plot.title_fontsize)

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "cpu_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path, transparent=True)
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

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.figsize)
        plt.clf()

        # Set background
        set_background(plt.gcf(), plt.gca(), self.config.plot)

        # Create a set that stores the tick labels for the plot
        ticks = set()
        ticks_ncores = set()

        # Plot handles
        handles = defaultdict(list)

        # Check whether we have multiple parameter sets
        multi_parameter_sets = len(self.memory_data[phase]) > 1

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.memory_data[phase]:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_memory: continue

            # Get the serial memory usage (and error) for this phase (create a Quantity object)
            if self.config.normalize_memory:
                serial_memory = self.serial_memory[parameter_set][phase].memory
                serial_error = self.serial_memory[parameter_set][phase].error
                serial = Measurement(serial_memory, serial_error)
            else: serial = None

            # Loop over the different parallelization modes (the different curves)
            for mode in sorted(self.memory_data[phase][parameter_set].keys()): # sort alphabetically

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.memory_data[phase][parameter_set][mode].processor_counts
                memories = self.memory_data[phase][parameter_set][mode].memory
                errors = self.memory_data[phase][parameter_set][mode].errors

                # Sort the lists
                processor_counts, memories, errors = sort_lists(processor_counts, memories, errors, to_arrays=True)

                # Get the array of x values
                if self.config.hybridisation: x_values = processor_counts
                # actually, in hybridisation mode,
                # processor_counts contains actually the process counts
                else:
                    x_values = get_x_values(processor_counts, mode, self.config.x_quantity)
                    if x_values is None: continue

                # Calculate the gains and the errors on the gains
                if self.config.normalize_memory:

                    memories_new = []
                    errors_new = []
                    for i in range(len(processor_counts)):

                        # Create a quantity for the current memory usage
                        memory = Measurement(memories[i], errors[i])

                        # Calculate the efficiency based on the current memory usage and the serial memory usage
                        memory_new = memory / serial

                        # Add the value and the propagated error of the normalized memory
                        memories_new.append(memory_new.value)
                        errors_new.append(memory_new.error)

                    # Convert to arrays
                    memories = np.array(memories_new)
                    errors = np.array(errors_new)

                # Determine the label
                if len(self.parameters_memory_left_after_ignored) > 0: label = mode + parameter_set_to_string_for_label(parameter_set, self.different_parameters_timing)
                else: label = mode

                # Plot the data points for this mode
                fmt = '' if self.config.plot.connect_points else 'o'
                handle = plt.errorbar(x_values, memories, errors, marker='.', label=label,
                                      linewidth=self.config.plot.linewidth, fmt=fmt, markersize=12)

                # Add the handle
                handles[parameter_set].append(handle)

                # Add the appropriate ticks
                ticks |= set(x_values)
                ticks_ncores |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Check if there less than 2 data points
        if len(ticks) < 2:
            log.warning("Not enough data for plotting")
            return

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks_ncores = sorted(ticks_ncores)
        #ticks.append(ticks[-1] * 2)

        # data model: self.memory_fit_parameters[phase][parameter_set][mode]

        # Plot the fit
        if not self.config.hybridisation and self.config.fit and self.config.fitting.plot_fit:

            # Loop over the parameter sets
            for parameter_set in self.memory_fit_parameters[phase]:

                # Plot the fitted curves
                fit_x_values = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
                fit_ncores = np.logspace(np.log10(ticks_ncores[0]), np.log10(ticks_ncores[-1]), 50)
                for mode in self.memory_fit_parameters[phase][parameter_set]:

                    # Get the parameter values
                    a = self.memory_fit_parameters[phase][parameter_set][mode].a
                    b = self.memory_fit_parameters[phase][parameter_set][mode].b
                    c = self.memory_fit_parameters[phase][parameter_set][mode].c

                    # Calculate the fitted memory usages
                    fit_memories = [modified_memory_scaling(nprocesses_from_mode_single(mode, ncores), a, b, c) for ncores in fit_ncores]

                    # Add the plot
                    plt.plot(fit_x_values, fit_memories, color="grey")

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        #ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlim(ticks[0], ticks[-1])

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set borders
        set_borders(ax, self.config.plot)

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes", fontsize=self.config.plot.label_fontsize)
        else: plt.xlabel("Number of " + self.config.x_quantity, fontsize=self.config.plot.label_fontsize)
        plt.ylabel("Memory usage per process (GB)", fontsize=self.config.plot.label_fontsize)

        # Add the legends
        # if self.config.hybridisation: plt.legend(title="Number of cores")
        # else: plt.legend(title="Parallelization modes")
        add_legends(ax, handles, self.different_parameters_memory, self.config.plot, "memory", "parameters")

        # Set the plot title
        if self.config.plot.add_titles:
            title = "Scaling of the memory usage of the " + phase_names[phase].lower()
            title = split_in_lines(title)
            plt.suptitle(title, fontsize=self.config.plot.title_fontsize)

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "memory_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path, transparent=True)
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

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.figsize)
        plt.clf()

        # Set background
        set_background(plt.gcf(), plt.gca(), self.config.plot)

        # Create a set that stores the tick labels for the plot
        ticks = set()
        ticks_ncores = set()

        # The plot handles
        handles = defaultdict(list)

        # Check whether we have multiple parameter sets
        multi_parameter_sets = len(self.memory_data[phase]) > 1

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.memory_data[phase]:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_memory: continue

            # Get the serial memory (and error) for this phase (create a Quantity object)
            serial_memory = self.serial_memory[parameter_set][phase].memory
            serial_error = self.serial_memory[parameter_set][phase].error
            serial = Measurement(serial_memory, serial_error)

            # Loop over the different parallelization modes (the different curves)
            for mode in sorted(self.memory_data[phase][parameter_set].keys()): # sort alphabetically

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.memory_data[phase][parameter_set][mode].processor_counts
                memories = self.memory_data[phase][parameter_set][mode].memory
                errors = self.memory_data[phase][parameter_set][mode].errors

                # Sort the lists
                processor_counts, memories, errors = sort_lists(processor_counts, memories, errors, to_arrays=True)

                # Get the array of x values
                if self.config.hybridisation: x_values = processor_counts
                # actually, in hybridisation mode,
                # processor_counts contains actually the process counts
                else:
                    x_values = get_x_values(processor_counts, mode, self.config.x_quantity)
                    if x_values is None: continue

                # Calculate the gains and the errors on the gains
                gains = []
                gain_errors = []
                for i in range(len(processor_counts)):

                    # Create a quantity for the current memory usage
                    memory = Measurement(memories[i], errors[i])

                    # Calculate the efficiency based on the current memory usage and the serial memory usage
                    gain = serial / memory

                    # Add the value and the propagated error of the gain to the appropriate lists
                    gains.append(gain.value)
                    gain_errors.append(gain.error)

                # Determine the label
                if len(self.parameters_memory_left_after_ignored) > 0: label = mode + parameter_set_to_string_for_label(parameter_set, self.different_parameters_timing)
                else: label = mode

                # Plot the data points for this mode
                handle = plt.errorbar(x_values, gains, gain_errors, marker='.', label=label, linewidth=self.config.plot.linewidth)

                # Add the handle
                handles[parameter_set].append(handle)

                # Add the appropriate ticks
                ticks |= set(x_values)
                ticks_ncores |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale("log")
        plt.yscale("log")

        # Check if there less than 2 data points
        if len(ticks) < 2:
            log.warning("Not enough data for plotting")
            return

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks_ncores = sorted(ticks_ncores)
        #ticks.append(ticks[-1] * 2)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        #ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlim(ticks[0], ticks[-1])

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set borders
        set_borders(ax, self.config.plot)

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes", fontsize=self.config.plot.label_fontsize)
        else: plt.xlabel("Number of " + self.config.x_quantity, fontsize=self.config.plot.label_fontsize)
        plt.ylabel("Memory gain", fontsize=self.config.plot.label_fontsize)

        # Add the legends
        # if self.config.hybridisation: plt.legend(title="Number of cores")
        # else: plt.legend(title="Parallelization modes")
        add_legends(ax, handles, self.different_parameters_memory, self.config.plot, "memory-gain", "parameters")

        # Set the plot title
        if self.config.plot.add_titles:
            title = "Scaling of the memory gain of the " + phase_names[phase].lower()
            title = split_in_lines(title)
            plt.suptitle(title, fontsize=self.config.plot.title_fontsize)

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "memorygain_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path, transparent=True)
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

        # Initialize figure with the appropriate size
        plt.figure(figsize=self.figsize)
        plt.clf()

        # Set background
        set_background(plt.gcf(), plt.gca(), self.config.plot)

        # Create a set that stores the tick labels for the plot
        ticks = set()
        ticks_ncores = set()

        # The plot handles
        handles = defaultdict(list)

        # Check whether we have multiple parameter sets
        multi_parameter_sets = len(self.memory_data[phase]) > 1

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.memory_data[phase]:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_memory: continue

            # Loop over the different parallelization modes (the different curves)
            for mode in sorted(self.memory_data[phase][parameter_set].keys()): # sort alphabetically

                # Get the list of processor counts, runtimes and errors
                processor_counts = self.memory_data[phase][parameter_set][mode].processor_counts
                memories = self.memory_data[phase][parameter_set][mode].memory
                errors = self.memory_data[phase][parameter_set][mode].errors

                # Sort the lists
                processor_counts, memories, errors = sort_lists(processor_counts, memories, errors, to_arrays=True)

                # Get list of process count
                if self.config.hybridisation: processes = processor_counts
                else: processes = nprocesses_from_mode(mode, processor_counts)

                # Get the array of x values
                if self.config.hybridisation: x_values = processor_counts
                # actually, in hybridisation mode,
                # processor_counts contains actually the process counts
                else:
                    x_values = get_x_values(processor_counts, mode, self.config.x_quantity)
                    if x_values is None: continue

                # Multiply the memory usage for each processor count with the corresponding number of processes (to get the total)
                memories *= processes
                errors *= processes

                # Determine the label
                if len(self.parameters_memory_left_after_ignored) > 0: label = mode + parameter_set_to_string_for_label(parameter_set, self.different_parameters_timing)
                else: label = mode

                # Plot the data points for this mode
                handle = plt.errorbar(x_values, memories, errors, marker='.', label=label, linewidth=self.config.plot.linewidth)

                # Add the plot handle
                handles[parameter_set].append(handle)

                # Add the appropriate ticks
                ticks |= set(x_values)
                ticks_ncores |= set(processor_counts)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale("log")
        plt.yscale("log")

        # Check if there less than 2 data points
        if len(ticks) < 2:
            log.warning("Not enough data for plotting")
            return

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        #ticks.append(ticks[-1] * 2)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        #ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        plt.xlim(ticks[0], ticks[-1])

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set borders
        set_borders(ax, self.config.plot)

        # Add axis labels and a legend
        if self.config.hybridisation: plt.xlabel("Number of processes", fontsize=self.config.plot.label_fontsize)
        else: plt.xlabel("Number of " + self.config.x_quantity, fontsize=self.config.plot.label_fontsize)
        plt.ylabel("Total memory usage (all processes) (GB)", fontsize=self.config.plot.label_fontsize)

        # Add the legends
        # if self.config.hybridisation: plt.legend(title="Number of cores")
        # else: plt.legend(title="Parallelization modes")
        add_legends(ax, handles, self.different_parameters_memory, self.config.plot, "total-memory", "parameters")

        # Set the plot title
        if self.config.plot.add_titles:
            title = "Memory scaling for " + phase_names[phase].lower()
            title = split_in_lines(title)
            plt.suptitle(title, fontsize=self.config.plot.title_fontsize)

        # Determine file path
        if self.config.output is not None: file_path = fs.join(self.config.output, "totalmemory_" + phase + ".pdf")
        else: file_path = None

        # Save the figure
        if file_path is not None: plt.savefig(file_path, transparent=True)
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

        # Loop over the different parameter sets (different ski files)
        for parameter_set in self.timing_data["total"]:

            # Ignore
            if parameter_set in self.ignore_parameter_sets_timing: continue

            # Debugging
            log.debug("Plotting the scaling timeline for " + parameter_set_to_string_inline(parameter_set, self.different_parameters_timing) + " ...")

            # Loop over the different parallelization modes
            for mode in sorted(self.timing_data["total"][parameter_set].keys()): # sort alphabetically

                # Determine plot path
                if self.config.output is not None: plot_file_path = fs.join(self.config.output, "timeline" + parameter_set_to_string_for_filename(parameter_set, self.different_parameters_timing) + "_" + mode + ".pdf")
                else: plot_file_path = None

                # Debugging
                if plot_file_path is not None: log.debug("Saving timeline plot file for the " + mode + " mode to '" + plot_file_path + "' ...")

                # Initialize a data structure to contain the start times and endtimes of the different simulation phases,
                # for the different processor counts (data is indexed on the simulation phase)
                data = []
                ncores_list = []
                nprocs_list = []

                # Initialize
                data.append(["setup", [], []])
                data.append(["stellar", [], []])
                data.append(["spectra", [], []])
                data.append(["dust", [], []])
                data.append(["write", [], []])
                data.append(["wait", [], []])
                data.append(["comm", [], []])

                # Get serial number of processors
                serial_nproc = self.serial_timing_ncores[parameter_set]["total"]

                # Add serial CPU times
                if self.config.timelines.add_serial and serial_nproc not in self.timing_data["total"][parameter_set][mode].processor_counts:

                    # Get the runtimes
                    setup_time = self.serial_timing[parameter_set]["setup"].time * serial_nproc
                    stellar_time = self.serial_timing[parameter_set]["stellar"].time * serial_nproc
                    spectra_time = self.serial_timing[parameter_set]["spectra"].time * serial_nproc
                    dust_time = self.serial_timing[parameter_set]["dust"].time * serial_nproc
                    writing_time = self.serial_timing[parameter_set]["writing"].time * serial_nproc
                    waiting_time = self.serial_timing[parameter_set]["waiting"].time * serial_nproc
                    communication_time = self.serial_timing[parameter_set]["communication"].time * serial_nproc

                    # Add number of serial cores
                    ncores_list.append(serial_nproc)

                    # Add row
                    add_timeline_row(data, setup_time, stellar_time, spectra_time, dust_time, writing_time, waiting_time, communication_time)

                # Loop over the different processor counts
                for j in range(len(self.timing_data["total"][parameter_set][mode].processor_counts)):

                    # Get the processor count
                    if self.config.hybridisation:
                        processors = int(mode.split(" cores")[0])
                        processes = self.timing_data["total"][parameter_set][mode].processor_counts[j]
                    else:
                        processors = self.timing_data["total"][parameter_set][mode].processor_counts[j]
                        processes = nprocesses_from_mode_single(mode, processors)

                    #print(self.timing_data["total"][parameter_set][mode].processor_counts)
                    #print(self.timing_data["setup"][parameter_set][mode].processor_counts)

                    # Get the average runtimes for the different phases corresponding to the current processor count
                    setup_time = self.timing_data["setup"][parameter_set][mode].times[j] * processors
                    stellar_time = self.timing_data["stellar"][parameter_set][mode].times[j] * processors
                    spectra_time = self.timing_data["spectra"][parameter_set][mode].times[j] * processors
                    dust_time = self.timing_data["dust"][parameter_set][mode].times[j] * processors
                    writing_time = self.timing_data["writing"][parameter_set][mode].times[j] * processors
                    waiting_time = self.timing_data["waiting"][parameter_set][mode].times[j] * processors
                    communication_time = self.timing_data["communication"][parameter_set][mode].times[j] * processors

                    # Add the processor and process count
                    ncores_list.append(processors)
                    nprocs_list.append(processes)

                    # Add row
                    add_timeline_row(data, setup_time, stellar_time, spectra_time, dust_time, writing_time, waiting_time, communication_time)

                # Set the plot title
                if self.config.plot.add_titles: title = "Scaling timeline"
                else: title = None

                # Create the plot
                if self.config.hybridisation: create_timeline_plot(data, nprocs_list, plot_file_path,
                                                                   percentages=self.config.timelines.percentages,
                                                                   totals=True, unordered=True, cpu=True, title=title,
                                                                   rpc='p', add_border=self.config.plot.add_border,
                                                                   label_fontsize=self.config.plot.label_fontsize,
                                                                   figsize=self.figsize,
                                                                   title_fontsize=self.config.plot.title_fontsize,
                                                                   ticks_fontsize=self.config.plot.ticks_fontsize)
                else: create_timeline_plot(data, ncores_list, plot_file_path,
                                           percentages=self.config.timelines.percentages, totals=True, unordered=True,
                                           cpu=True, title=title, rpc='c', add_border=self.config.plot.add_border,
                                           label_fontsize=self.config.plot.label_fontsize, figsize=self.figsize,
                                           title_fontsize=self.config.plot.title_fontsize, ticks_fontsize=self.config.plot.ticks_fontsize)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the timing table
        if self.needs_timing and self.timing is not None: self.write_timing()

        # Write the memory table
        if self.needs_memory and self.memory is not None: self.write_memory()

        # Write the timing data
        if self.needs_timing: self.write_timing_data()

        # Write the memory data
        if self.needs_memory: self.write_memory_data()

        # Write timing fit data
        if self.needs_timing: self.write_timing_fits()

        # Write memory fit data
        if self.needs_memory: self.write_memory_fits()

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

    def write_timing_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the timing data ...")

        # Write
        path = fs.join(self.config.output, "timing_data.dat")
        write_dict(self.timing_data, path)

        # Write serial
        path = fs.join(self.config.output, "serial_timing_data.dat")
        write_dict(self.serial_timing, path)

        # Write serial ncores
        path = fs.join(self.config.output, "serial_timing_ncores.dat")
        write_dict(self.serial_timing_ncores, path)

        # Write different parameters
        path = fs.join(self.config.output, "different_parameters_timing.dat")
        write_list(self.different_parameters_timing, path)

        # Write ignore
        path = fs.join(self.config.output, "ignore_parameter_sets_timing.dat")
        write_list(list(self.ignore_parameter_sets_timing), path)

    # -----------------------------------------------------------------

    def write_memory_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the memory data ...")

        # Write
        path = fs.join(self.config.output, "memory_data.dat")
        write_dict(self.memory_data, path)

        # Write serial
        path = fs.join(self.config.output, "serial_memory_data.dat")
        write_dict(self.serial_memory, path)

        # Write serial ncores
        path = fs.join(self.config.output, "serial_memory_ncores.dat")
        write_dict(self.serial_memory_ncores, path)

        # Write different parameters
        path = fs.join(self.config.output, "different_parameters_memory.dat")
        write_list(self.different_parameters_memory, path)

        # Write ignore
        path = fs.join(self.config.output, "ignore_parameter_sets_memory.dat")
        write_list(list(self.ignore_parameter_sets_memory), path)

    # -----------------------------------------------------------------

    def write_timing_fits(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the timing fit results ...")

        # Loop over the phases
        for phase in self.timing_fit_parameters:

            # Loop over the parameter sets
            for parameter_set in self.timing_fit_parameters[phase]:

                # Get the modes
                modes = self.timing_fit_parameters[phase][parameter_set].keys()

                mode_list = []
                other_columns = defaultdict(list)

                # Fill columns
                for mode in modes:

                    mode_list.append(mode)

                    fit_parameters = self.timing_fit_parameters[phase][parameter_set][mode]
                    fit_parameter_errors = self.timing_fit_parameter_errors[phase][parameter_set][mode]

                    for i in range(len(fit_parameters)):

                        other_columns[alphabet[i]].append(fit_parameters[i])
                        other_columns[alphabet[i] + "_error"].append(fit_parameter_errors[i])

                # Create a data file to contain the fitted parameters
                directory = self.config.output
                parameter_file_path = fs.join(directory, "parameters_timing_" + parameter_set_to_string_for_filename(parameter_set, self.different_parameters_timing) + "_" + phase + ".dat")

                # Create the parameters table and write to file
                data = [mode_list]
                names = ["Parallelization mode"]

                # Set data and names
                for column_name in other_columns:
                    data.append(other_columns[column_name])
                    names.append(column_name)

                # Write the table
                table = Table(data=data, names=names)
                tables.write(table, parameter_file_path)

    # -----------------------------------------------------------------

    def write_memory_fits(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the memory fit results ...")

        # Loop over the phases
        for phase in self.memory_fit_parameters:

            # Loop over the parameter sets
            for parameter_set in self.memory_fit_parameters[phase]:

                # Get the modes
                modes = self.memory_fit_parameters[phase][parameter_set].keys()

                mode_list = []
                a_list = []
                a_error_list = []
                b_list = []
                b_error_list = []
                c_list = []
                c_error_list = []

                for mode in modes:

                    parameters = self.memory_fit_parameters[phase][parameter_set][mode]

                    mode_list.append(mode)
                    a_list.append(parameters.a)
                    a_error_list.append(parameters.a_error)
                    b_list.append(parameters.b)
                    b_error_list.append(parameters.b_error)
                    c_list.append(parameters.c)
                    c_error_list.append(parameters.c_error)

                # Create a data file to contain the fitted parameters
                directory = self.config.output
                parameter_file_path = fs.join(directory, "parameters_memory_" + parameter_set_to_string_for_filename(parameter_set, self.different_parameters_memory) + "_" + phase + ".dat")

                # Create the parameters table and write to file
                data = [mode_list, a_list, a_error_list, b_list, b_error_list, c_list, c_error_list]
                names = ["Parallelization mode", "Parameter a", "Error on a", "Parameter b", "Error on b", "Parameter c",
                         "Error on c"]
                table = Table(data=data, names=names)
                tables.write(table, parameter_file_path)

# -----------------------------------------------------------------

def nthreads_per_process_for_hybrid_mode(mode):

    """
    This function ...
    :param mode:
    :return:
    """

    threads = int(mode.split("(")[1].split(" threads)")[0])
    return threads

# -----------------------------------------------------------------

def ncores_from_mode_single(mode, nprocesses):

    """
    This function ...
    :param mode:
    :param nprocesses:
    :return:
    """

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

def ncores_from_hybridisation_mode(mode):

    """
    This function ...
    :param mode:
    :return:
    """

    return int(mode.split(" cores")[0])

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

def generate_scaling_behaviour(phase, config):

    """
    This function ...
    :param phase:
    :param config:
    :return:
    """

    behaviour = set()

    # Pure scaling term
    if phase in config.pure_scaling_behaviour:
        for component in config.pure_scaling_behaviour[phase]: behaviour.add(component)

    # Composite scaling term
    elif phase in config.composite_scaling_behaviour:
        for composite_phase in config.composite_scaling_behaviour[phase]:
            for component in generate_scaling_behaviour(composite_phase, config): behaviour.add(component)

    # Derived (all specific communication phases derive the same function from the total communication 'phase')
    elif phase in derived_scaling_behaviour:
        for component in config.pure_scaling_behaviour[derived_scaling_behaviour[phase]]: behaviour.add(component)

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

def get_hybridisation_mode(ncores, data_parallel):

    """
    This function ...
    :return:
    """

    if data_parallel: return str(ncores) + " cores (task+data)"
    else: return str(ncores) + " cores (task)"

# -----------------------------------------------------------------

def parameter_set_to_string_inline(parameter_set, parameter_names):

    """
    This function ...
    :param parameter_set:
    :param parameter_names:
    :return:
    """

    # string that says 'empty' (no different parameters for the simulations)
    if parameter_set == "empty": string = "none"

    # tuple
    elif isinstance(parameter_set, tuple): string = ", ".join([parameter_names[index] + ": " + stringify.stringify_not_list(parameter_set[index], scientific=True)[1] for index in range(len(parameter_names))])

    # only one value
    else: string = parameter_names[0] + ": " + stringify.stringify_not_list(parameter_set, scientific=True)[1]

    # return
    return string

# -----------------------------------------------------------------

def parameter_set_to_string_for_label(parameter_set, parameter_names):

    """
    This function ...
    :param parameter_set:
    :param parameter_names:
    :return:
    """

    # string that says 'empty' (no different parameters for the simulations)
    if parameter_set == "empty": string = ""

    # tuple
    elif isinstance(parameter_set, tuple): string = " [" + ", ".join([parameter_names[index] + ": " + stringify.stringify_not_list(parameter_set[index], scientific=True)[1] for index in range(len(parameter_names))]) + "]"

    # only one value
    else: string = " [" + parameter_names[0] + ": " + stringify.stringify_not_list(parameter_set, scientific=True)[1] + "]"

    # return
    return string

# -----------------------------------------------------------------

def parameter_set_to_string_for_legend(parameter_set, parameter_names):

    """
    This function ...
    :param parameter_set:
    :param parameter_names:
    :return:
    """

    # string that says 'empty'
    if parameter_set == "empty": string = None

    # tuple
    elif isinstance(parameter_set, tuple): string = "\n".join([parameter_names[index] + ": " + stringify.stringify_not_list(parameter_set[index], scientific=True)[1] for index in range(len(parameter_names))])

    # only one value
    else: string = parameter_names[0] + ": " + stringify.stringify_not_list(parameter_set, scientific=True)[1]

    # Return
    return string

# -----------------------------------------------------------------

def parameter_set_to_string_for_filename(parameter_set, parameter_names):

    """
    This function ...
    :param parameter_set:
    :param parameter_names:
    :return:
    """

    if parameter_set == "empty": string = ""

    elif isinstance(parameter_set, tuple):

        string = "_" + "_".join([parameter_names[index] + "=" + stringify.stringify_not_list(parameter_set[index], scientific=True)[1] for index in range(len(parameter_names))])

    # Only one value
    else: string = "_" + parameter_names[0] + "=" + stringify.stringify_not_list(parameter_set, scientific=True)[1]

    # Return
    return string

# -----------------------------------------------------------------

def all_equal(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    first = lst[0]

    for index in range(1,len(lst)):
        if lst[index] != first: return False

    return True

# -----------------------------------------------------------------

def add_timeline_row(data, setup_time, stellar_time, spectra_time, dust_time, writing_time, waiting_time,
                     communication_time):

    """
    This function ...
    :param data:
    :param setup_time:
    :param stellar_time:
    :param spectra_time:
    :param dust_time:
    :param writing_time:
    :param waiting_time:
    :param communication_time:
    :return:
    """

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

# -----------------------------------------------------------------

def behaviour_to_function_string(behaviour):

    """
    This function ...
    :param behaviour:
    :return:
    """

    string = "T(N_c) ="

    alphabet_iterator = iter(alphabet)

    terms = []

    for power in behaviour:

        if power == "log": terms.append(alphabet_iterator.next() + " x log(N_c)")
        elif power == -1: terms.append(alphabet_iterator.next() + " / N_c")
        elif power == 0: terms.append(alphabet_iterator.next())
        elif power == 1: terms.append(alphabet_iterator.next() + " x N_c")
        else: terms.append("N_c^" + str(power))

    string += " " + " + ".join(terms)

    return string

# -----------------------------------------------------------------

def add_legends(ax, handles, different_parameters, config, property, pa_or_ph):

    """
    This function ...
    :param ax:
    :param handles:
    :param different_parameters:
    :param config:
    :param property:
    :param pa_or_ph:
    :return:
    """

    nlegends = len(handles)

    if config.legend_below:

        if nlegends > 1: percentage = 25.
        else: percentage = 10.

        # Shrink current axis's height by a certain percentage on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * percentage/100., box.width, box.height * (100-percentage)/100.])

    # Define inverse transform, transforms display coordinates (pixels) to axes coordinates
    inv = ax.transAxes.inverted()

    # Find how may pixels there are on the x-axis
    x_pixels = ax.transAxes.transform((1, 0)) - ax.transAxes.transform((0, 0))

    # Find how many pixels there are on the y-axis
    y_pixels = ax.transAxes.transform((0, 1)) - ax.transAxes.transform((0, 0))

    axbox = ax.get_position()

    x_value = 0.0
    y_value = 0.0

    plot_width = axbox.width
    plot_height = axbox.height


    # Generate one legend now to get the dimensions of the legend
    #legend = plt.legend(handles=handles[handles.keys()[0]], title="dummy title", loc="upper center")
    #legend_ax = ax.add_artist(legend)
    #plt.draw()

    #legend.draw()
    #p = legend.get_window_extent()
    #print(p)

    #legbounds = legend.get_frame().get_bbox().bounds
    #print(legbounds)

    #legbox = legend.get_bbox_to_anchor()
    #legheight = legbox.height
    #legwidth = legbox.width

    #print()
    #print(plot_width, plot_height)
    #print(legwidth, legheight)
    #print(x_pixels, y_pixels)
    #print()

    # Loop over the parameter sets
    height = 0.0
    width = 0.0
    width_per_legend = 1. / len(handles) * 1.1
    for index in range(len(handles)):

        # Get the parameter set (or phase)
        parameter_set = handles.keys()[index]

        # Determine the legend title
        if pa_or_ph == "parameters": legend_title = parameter_set_to_string_for_legend(parameter_set, different_parameters)
        else: legend_title = phase_names_for_legend[parameter_set].capitalize()
        if legend_title is not None: legend_title =  r"\underline{" + legend_title + "}"

        # Create a legend for this parameter set
        if config.legend_below:

            if nlegends > 1:
                legend = plt.legend(handles=handles[parameter_set], loc="upper center", title=legend_title,
                                    bbox_to_anchor=(0.2 + width_per_legend * index, -0.15), fancybox=False, shadow=False,
                                    ncol=1)
            else:
                legend = plt.legend(handles=handles[parameter_set], loc="upper center", title=legend_title,
                                    bbox_to_anchor=(0.5, -0.11), fancybox=False, shadow=False, ncol=2)

        elif nlegends > 1:

            if properties_behaviour[property] == "ascending":
                #legend = plt.legend(handles=handles[parameter_set], title=legend_title, loc=(axbox.x0 + x_value + 0.65*plot_width, 0.05 + y_value + index * 0.25 * height))
                legend = plt.legend(handles=handles[parameter_set], title=legend_title, loc=(axbox.x0 + x_value + plot_width - 0.5*legwidth, 0.5*legheight + y_value + index * 0.25 * height))
            elif properties_behaviour[property] == "descending":
                #legend = plt.legend(handles=handles[parameter_set], title=legend_title, loc=(axbox.x0 + x_value + 0.65*plot_width - index * 0.5 * width, 0.05 + y_value + 0.95 * plot_height))
                legend = plt.legend(handles=handles[parameter_set], title=legend_title, loc=(axbox.x0 + x_value + plot_width - 0.5*legwidth - index * 0.5 * width, 0.05 + y_value + plot_height - 0.5*legheight))

            if index == 0:
                first_box = legend.axes.get_position()
                height = first_box.height
                width = first_box.width

        else: legend = plt.legend(handles=handles[parameter_set], title=legend_title, loc='best')

        frame = legend.get_frame()

        # Set legend frame color and line width
        if config.add_legend_border: frame.set_linewidth(config.legend_borderwidth)
        else: frame.set_linewidth(0)
        frame.set_edgecolor(config.legend_bordercolor)

        # Set background color
        frame.set_facecolor('0.85')
        legend.legendPatch.set_alpha(0.75)

        # Move to foreground
        legend.set_zorder(100+index)

        #  Add the legend manually to the current axes (except when it is the last)
        if index != len(handles) - 1:
            legend_ax = ax.add_artist(legend)
            #legend_ax.set_zorder(100+index)

        # Set fontsize
        plt.setp(legend.get_title(), fontsize=str(config.legend_title_fontsize))

# -----------------------------------------------------------------

def set_borders(ax, config):

    """
    This function ...
    :param ax:
    :param config:
    :return:
    """

    # Set border width
    if config.add_borders: [i.set_linewidth(config.borderwidth) for i in ax.spines.itervalues()]

    # Remove borders
    else:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis=u'both', which=u'both', length=0)

# -----------------------------------------------------------------

def set_grid(config):

    """
    This function ...
    :param config:
    :return:
    """

    if config.add_grid: plt.grid(linewidth=config.grid_linewidth, linestyle=config.grid_linestyle, color=config.grid_color)

# -----------------------------------------------------------------

def set_background(fig, ax, config):

    """
    This function ...
    :param fig:
    :param ax:
    :param config:
    :return:
    """

    if config.transparent_background:

        # Set transparent background
        for item in [fig, ax]:
            item.patch.set_visible(False)

# -----------------------------------------------------------------

def get_x_values(processor_counts, mode, x_quantity):

    """
    This function ...
    :param processor_counts:
    :param mode:
    :param x_quantity:
    :return:
    """

    # Convert x quantity if necessary
    # Cores on x axis
    if x_quantity == "cores": x_values = processor_counts

    # Processes on x axis
    elif x_quantity == "processes":

        # Skip multithreading modes
        if mode == "multithreading":
            log.warning("Cannot plot the runtimes for " + mode + " mode when the number of processes is to be used for the x axis")
            return None

        # Get list of nprocesses
        x_values = nprocesses_from_mode(mode, processor_counts)

    # Threads on x axis
    elif x_quantity == "threads":

        # Skip multiprocessing and hybrid modes
        if mode.startswith("multiprocessing"):
            log.warning("Cannot plot the runtimes for " + mode + " mode when the number of threads is to be used for the x axis")
            return None

        elif mode.startswith("hybrid"):
            log.warning("Cannot plot the runtimes for " + mode + " mode when the number of threads is to be used for the x axis")
            return None

        # Get list of nthreads = ncores in multithreading mode
        x_values = processor_counts

    # Invalid
    else: raise ValueError("Invalid x quantity")

    # Return
    return x_values

# -----------------------------------------------------------------
