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

# Timing and memory table
definition.add_optional("timing", "file_path", "timing table path")
definition.add_optional("memory", "file_path", "memory table path")

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

# -----------------------------------------------------------------

# Showing
definition.add_flag("show", "showing", True)
definition.add_flag("show_assignment", "show the assignment scheme")
definition.add_flag("show_status", "show the simulation status")
definition.add_flag("show_runtimes", "show runtimes")
definition.add_flag("show_memory", "show memory")

# -----------------------------------------------------------------

# Plotting
definition.add_flag("plot", "plotting", True)
definition.add_flag("plot_runtimes", "plot runtimes")
definition.add_flag("plot_memory", "plot memory usage")

# -----------------------------------------------------------------

# Writing
definition.add_flag("write", "writing", True)
definition.add_flag("write_assignment", "write the assignent scheme", None)
definition.add_flag("write_status", False)

# -----------------------------------------------------------------

# Analyse?
definition.add_flag("analyse", "analysis", False)

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
