#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.plot.scaling import scaling_properties, simulation_phases, communication_phases
from pts.core.plot.scaling import pure_scaling_behaviour, composite_scaling_behaviour, phase_names
from pts.core.config.plot import definition as plot_definition

# -----------------------------------------------------------------

parallelization_modes = ["multithreading", "multiprocessing", "hybrid"]

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Flags
definition.add_flag("recursive", "look for simulation in directories recursively", True)

# Add optional
definition.add_positional_optional("properties", "string_list", "plot the scaling of these properties", choices=scaling_properties, default=scaling_properties)
definition.add_positional_optional("phases", "string_list", "simulation phases for which to do the plotting", choices=simulation_phases, default=["total"])

# Timing and memory table paths
definition.add_optional("timing_table", "file_path", "path to the timing table file")
definition.add_optional("memory_table", "file_path", "path to the memory table file")

# Data input directory
definition.add_optional("data_input", "directory_path", "path to the directory where plotting data has been written out, for re-use")

definition.add_flag("hybridisation", "plot as a function of number of processes for constant number of cores")

definition.add_optional("modes", "string_list", "parallelization modes to plot (not for hybridization mode)", choices=parallelization_modes, default=parallelization_modes)
definition.add_flag("use_task_parallel", "use data from task parallelized simulations in the plots", True)
definition.add_flag("use_task_data_parallel", "use data from task+data parallelized simulations in the plots", True)

# Input
definition.add_optional("input", "directorypath_list", "input director(y)(ies)", letter='i')

# Output
definition.add_optional("output", "directory_path", "output directory", letter="o")

# Report which simulations are used for the plotting
definition.add_flag("report_simulations", "report which simulations are used for the data")

# Plotting options
definition.import_section("plot", "plotting options", plot_definition)

definition.add_flag("xlog", "log scale for x axis", True)
definition.add_flag("ylog", "log scale for y axis", True)

# The quantity to be used for the x axis
definition.add_optional("x_quantity", "string", "quantity to be used for the x axis (not for hybridization mode)", choices=["cores", "processes", "threads"], default="cores")

# Sigma level for plotting error bars
definition.add_optional("sigma_level", "real", "sigma level for plotting error bars", 1.0)

# Make fit and plot fit
definition.add_flag("fit", "fit theoretical curves to timing and memory data", True)

# Fitting options
definition.add_section("fitting", "fitting options")
definition.sections["fitting"].add_flag("plot_fit", "plot the fitted relations", True)

# Pure scaling behaviour
definition.sections["fitting"].add_section("pure_scaling_behaviour", "pure scaling behaviour")
for phase in pure_scaling_behaviour:
    definition.sections["fitting"].sections["pure_scaling_behaviour"].add_optional(phase, "mixed_list", "scaling behaviour for the " + phase_names[phase], default=pure_scaling_behaviour[phase])

# Composite scaling behaviour
definition.sections["fitting"].add_section("composite_scaling_behaviour", "composite scaling behaviour")
for phase in composite_scaling_behaviour:
    definition.sections["fitting"].sections["composite_scaling_behaviour"].add_optional(phase, "string_list", "scaling behaviour for the " + phase_names[phase], default=composite_scaling_behaviour[phase])

# Normalize runtimes and memory usages
definition.add_flag("normalize_runtimes", "normalize runtimes for plotting")
definition.add_flag("normalize_memory", "normalize memory usage for plotting")

# Split communication into subphases
definition.add_flag("split_communication", "split the different communication steps")
definition.add_optional("communication_subphases", "string_list", "communication subphases to plot", communication_phases)

# Enable all properties and phases
definition.add_flag("all", "plot everything (enable all properties and phases)")
definition.add_flag("all_timing", "plot everything related to timing")
definition.add_flag("all_memory", "plot everything related to memory usage")

# Enable all phases
definition.add_flag("all_phases", "enable all phases")

# Enable all properties
definition.add_flag("all_properties", "enable all properties")

# Timelines
definition.add_section("timelines", "options for plotting timelines")
definition.sections["timelines"].add_flag("add_serial", "add CPU times of serial run for comparison")
definition.sections["timelines"].add_flag("percentages", "add percentages for the different phases compared to the total simulation", True)

# FROM HERE: ADVANCED STUFF: USE WITH CARE
definition.add_flag("hetero", "not necessarily a single ski")

# EXTRAPOLATION
definition.add_section("extrapolation", "extrapolate ...")

# TIMING
definition.sections["extrapolation"].add_section("timing", "extrapolation of timing data")
definition.sections["extrapolation"].sections["timing"].add_flag("ncores", "extrapolate the data to a number of cores of one to get a serial timing")
definition.sections["extrapolation"].sections["timing"].add_flag("npackages", "extrapolate the number of photon packages of a serial run to obtain a serial run for a series of simulations with a higher number of packages (requires 'hetero' to be enabled)")
definition.sections["extrapolation"].sections["timing"].add_flag("nwavelengths", "extrapolate the number of wavelengths of a serial run to obtain a serial run for a series of simulations with a higher number of wavelengths (requires 'hetero' to be enabled) [THIS OPTION IS VERY TRICKY: LOAD BALANCING CAN VARY!]")
definition.sections["extrapolation"].sections["timing"].add_flag("in_times", "not only extrapolate for normalizing, but also plot the serial time as if it were a genuine data point")
definition.sections["extrapolation"].sections["timing"].add_flag("for_hybrid", "also extrapolate to serial for hybrid runs for visualisation purposes (which can per definition not be run on 1 core)")

# MEMORY
definition.sections["extrapolation"].add_section("memory", "extrapolation of memory data")
definition.sections["extrapolation"].sections["memory"].add_flag("nprocesses", "extrapolate the data to a number of processes of one to get serial memory data points")
definition.sections["extrapolation"].sections["memory"].add_flag("nwavelengths", "extrapolate the number of wavelengths of a serial run to obtain a serial run for a series of simulations with a higher number of wavelengths (required 'hetero' to be enabled)")
definition.sections["extrapolation"].sections["memory"].add_flag("extrapolate_ncells", "extrapolate the number of dust cells")
definition.sections["extrapolation"].sections["memory"].add_flag("in_memory", "not only extrapolate for normalizing, but also plot the serial time as if it were a genuine data point")

# -----------------------------------------------------------------

definition.add_flag("tolerant", "don't compare timing and memory tables in a strict way", True)

# -----------------------------------------------------------------
