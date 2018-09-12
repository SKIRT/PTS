#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.remote.host import find_host_ids
from pts.core.config.manage_simulations import definition

# -----------------------------------------------------------------

# Load the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# -----------------------------------------------------------------

# Generation name
definition.add_required("generation", "string", "generation name")

# -----------------------------------------------------------------

# Plotting
definition.add_flag("plot_chisquared", "plot chi squared")

# Flags
definition.add_flag("lazy", "lazy mode")
definition.add_flag("find_simulations", "find missing simulations by searching on simulation name", True)
definition.add_optional("find_remotes", "string_list", "find missing simulations in these remote hosts", default=all_host_ids, choices=all_host_ids)
definition.add_flag("produce_missing", "produce missing simulation files", False)
definition.add_flag("check_paths", "check simulation paths", False)
definition.add_flag("correct_paths", "correct simulation paths instead of raising errors", False)
definition.add_flag("confirm_correction", "confirm before correcting paths", False)
definition.add_flag("check_analysis", "check analysis output", False)
definition.add_flag("check_status", "check the simulation status", True)

# Caching
definition.add_optional("cache_volume", "string", "name of the volume to be used for caching")

# Correct analysis status
definition.add_flag("correct_status", "correct the status", True)

# -----------------------------------------------------------------

# Interactive
definition.flags["interactive"].default = True

# -----------------------------------------------------------------
