#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.core.launch.analyser import all_steps
from pts.core.config.analyse_simulation import definition as analysis_definition

# -----------------------------------------------------------------

host_ids = find_host_ids()

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Assignment file
definition.add_positional_optional("assignment", "file_path", "path of assignment file")

# Status file
definition.add_optional("status", "file_path", "path of status file")

# Remotes for which to find corresponding simulations
definition.add_optional("remotes", "string_list", "remote hosts for which to look for matching simulations (not necessary when assignment is specified)", default=host_ids, choices=host_ids)

# To specify the simulations
definition.add_optional("simulation_names", "string_list", "names of the simulations to look for")
definition.add_optional("simulation_ids", "integer_list", "IDs of the simulations (only if one remote host is specified)")
definition.add_flag("from_directories", "use directory names as simulation names")

# Timing and memory table
definition.add_optional("timing", "file_path", "timing table path")
definition.add_optional("memory", "file_path", "memory table path")

# Commands to be run
definition.add_optional("commands", "string_list", "commands to be run in interactive mode")

# Interactive mode
definition.add_flag("interactive", "use interactive mode", default=None)

# Offline?
definition.add_flag("offline", "offline mode")

# Dry: don't actually launch any simulation
definition.add_flag("dry", "dry mode")

# Fix success flags
definition.add_flag("fix_success", "fix success flags in assignment table if necessary")

# Backup
definition.add_flag("backup_simulations", "backup simulation files")
definition.add_flag("backup_assignment", "backup assignment table")
definition.add_optional("backup_path", "directory_path", "backup directory")
definition.add_optional("backup_dir_path", "directory_path", "directory to contain the backup directory")
definition.add_optional("backup_dirname", "string", "name of the backup directory")

# -----------------------------------------------------------------

# Flags
definition.add_flag("local", "treat simulations without a match as local simulations", False)
definition.add_flag("warn_local", "give a warning for each simulation that didn't have a match for any host ID")
definition.add_flag("success", "success flag to fill in into the assignment table for all simulations", True)

# -----------------------------------------------------------------

# Move simulations
definition.add_flag("move", "move simulations", False)
definition.add_flag("move_running", "move running simulations", False)
definition.add_optional("move_simulations", "string_list", "simulation names for moving")
definition.add_optional("move_remotes", "string_list", "host IDs for moving simulations from")
definition.add_flag("prompt_simulations_move", "prompt before moving a particular simulation", None)
definition.add_optional("move_to_host", "host", "move simulations to this remote host")

# -----------------------------------------------------------------

# Showing
definition.add_flag("show", "showing", True)
definition.add_flag("show_assignment", "show the assignment scheme")
definition.add_flag("show_status", "show the simulation status")
definition.add_flag("show_runtimes", "show runtimes")
definition.add_flag("show_memory", "show memory")

# -----------------------------------------------------------------

runtimes_phases = ["total", "setup", "stellar", "spectra", "dust", "writing", "waiting", "communication", "intermediate"]
memory_phases = ["total", "setup", "stellar", "spectra", "dust", "writing"]

# Plotting
definition.add_flag("plot", "plotting", True)
definition.add_flag("plot_runtimes", "plot runtimes")
definition.add_optional("plot_runtimes_phases", "string_list", "simulation phases for which to plot the runtimes", runtimes_phases, choices=runtimes_phases)
definition.add_flag("plot_memory", "plot memory usage")
definition.add_optional("plot_memory_phases", "string_list", "simulation phases for which to plot the memory usage", memory_phases, choices=memory_phases)

# Observed SEDs as reference for plotting simulation SEDs
definition.add_optional("reference_seds", "string_filepath_dictionary", "file paths of SEDs to use as reference for plotting simulated SEDs")

# -----------------------------------------------------------------

# Caching simulation output
definition.add_optional("cache_path", "directory_path", "path to be used for caching")
definition.add_optional("cache_root", "directory_path", "path of directory up in the hierarchy with respect to the simulation directories that should be set equivalent to the 'cache_path' and the intermediate directory structure should be created if necessary")
definition.add_flag("cache_output", "cache the output of retrieved simulations")
definition.add_flag("cache_datacubes", "cache the datacubes of retrieved simulations")
definition.add_flag("cache_misc", "cache the misc output of analysed simulations")
definition.add_flag("cache_images", "cache the images output of analysed simulations")
definition.add_flag("cache_after_analysis", "cache immediately after a simulation has been analysed")

# -----------------------------------------------------------------

# Writing
definition.add_flag("write", "writing", True)
definition.add_flag("write_assignment", "write the assignent scheme", None)
definition.add_flag("write_status", "write the status table", False)
definition.add_flag("write_moved", "write the moved simulations table", False)
definition.add_flag("write_relaunched", "write the relaunched simulations table", False)
definition.add_flag("write_commands", "write the commands", False)

# -----------------------------------------------------------------

# Retrieve?
definition.add_flag("retrieve", "retrieve finished simulations", False)

# -----------------------------------------------------------------

# Analyse?
definition.add_flag("analyse", "run analysis", False)

# -----------------------------------------------------------------

# (Re)analyse only certain simulations
definition.add_optional("analyse_simulations", "string_list", "simulation names for analysis")
definition.add_optional("reanalyse_simulations", "string_list", "simulation names for re-analysis")

# (Re)analyse only for certain host IDs
definition.add_optional("analyse_remotes", "string_list", "host IDs for analysis")
definition.add_optional("reanalyse_remotes", "string_list", "host IDs for re-analysis")

# Prompt
definition.add_flag("prompt_simulations_analysis", "prompt before analysing a particular simulation (by default, all retrieved simulations are analysed without prompt)", False)
definition.add_flag("prompt_simulations_reanalysis", "prompt before re-analysing a particular simulation (by default, will be promted if renalaysis_simulations are not defined)", None)

# Add section for analysis options (also for re-analysis)
definition.import_section("analysis", "analyser options", analysis_definition)

# Re-analysis options
# steps: extraction, plotting, misc, batch, scaling
# features:
#  - extraction: progress, timeline, memory
#  - plotting: = progress, timeline, memory, seds, grids
#  - misc: = rgb, animations, fluxes, fluxes_from_images, images
definition.add_optional("reanalyse", "string_list", "apply re-analysis of these steps", choices=all_steps)
definition.add_optional("features_reanalysis", "string_list", "re-analyse only certain features (if a single re-analysis step is defined)")

# -----------------------------------------------------------------

definition.add_optional("info_tables", "filepath_list", "tables with additional information for the simulations")

# -----------------------------------------------------------------

definition.add_flag("info_scientific", "use scientific notation for formatting info values (default is automatic)", None)
definition.add_optional("info_ndecimal_places", "positive_integer", "number of decimal places for formatting info values", 3)

# -----------------------------------------------------------------

# Shared input?
definition.add_flag("shared_input", "whether the different simulations share their input files", None)

# -----------------------------------------------------------------

# List of local screen script paths
definition.add_optional("screen_scripts", "filepath_list", "filepaths of local screen scripts associated with the simulations")

# -----------------------------------------------------------------
