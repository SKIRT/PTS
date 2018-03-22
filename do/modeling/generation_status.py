#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_status View the status of the simulations of a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_yn
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.launch.manager import extra_columns
from pts.core.remote.host import find_host_ids
from pts.modeling.fitting.manager import GenerationManager

# -----------------------------------------------------------------

# Load the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation name
definition.add_required("generation", "string", "generation name")

# Show parameters
definition.add_optional("extra", "string_list", "show extra info", choices=extra_columns)

# Options for the manager
definition.add_flag("offline", "offline mode")
definition.add_flag("lazy", "lazy mode")
definition.add_flag("find_simulations", "find missing simulations by searching on simulation name", True)
definition.add_optional("find_remotes", "string_list", "find missing simulations in these remote hosts", default=all_host_ids, choices=all_host_ids)
definition.add_flag("produce_missing", "produce missing simulation files", False)
definition.add_flag("check_paths", "check simulation paths", False)
definition.add_flag("correct_paths", "correct simulation paths instead of raising errors", False)
definition.add_flag("confirm_correction", "confirm before correcting paths", False)
definition.add_flag("fix_success", "check success flags in assignment table")
definition.add_flag("check_analysis", "check analysis output", False)

# Write the status table
definition.add_flag("write_status", "write the status", False)

# Correct analysis status
definition.add_flag("correct_status", "correct the status", True)

# Get configuration
config = parse_arguments("generation_status", definition, "View the status of the simulations of a certain generation")

# -----------------------------------------------------------------

# Create simulation manager
manager = GenerationManager()

# Set fitting run and generation name
manager.config.run = config.run
manager.config.generation = config.generation

# Set options
manager.config.extra = config.extra
manager.config.offline = config.offline
manager.config.lazy = config.lazy
manager.config.find_simulations = config.find_simulations
manager.config.find_remotes = config.find_remotes
manager.config.produce_missing = config.produce_missing
manager.config.check_paths = config.check_paths
manager.config.correct_paths = config.correct_paths
manager.config.confirm_correction = config.confirm_correction
manager.config.fix_success = config.fix_success
manager.config.check_analysis = config.check_analysis
manager.config.write_status = config.write_status
manager.config.correct_status = config.correct_status

# Not interactive
manager.config.interactive = False

# Set the status command
if config.extra is not None: status_command = "status " + ",".join(config.extra)
else: status_command = "status"

# Set status command
manager.config.commands = [status_command]

# Run the generation manager
manager.run()

# -----------------------------------------------------------------
