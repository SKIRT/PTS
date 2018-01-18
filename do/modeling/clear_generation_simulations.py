#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.clear_generation_simulations Remove the simulation files for a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_proceed
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.core.tools.stringify import tostr
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Get the modeling commands
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation for which to remove the simulations (none means all)")
definition.add_positional_optional("remote", "string", "remote host for which to remove the simulations (none means all")

# Get configuration
config = parse_arguments("clear_generation_output", definition)

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Check
if not fitting_run.is_generation(config.generation): raise ValueError("Generation doesn't exist")

# Get generation path
generation_path = fitting_run.get_generation_path(config.generation)

# -----------------------------------------------------------------

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Loop over the remote hosts
for host_id in generation.host_ids:

    # Clear for this host?
    if config.remote is not None and host_id != config.remote: continue

    # Get paths of the simulation files
    filepaths = generation.get_simulation_paths_for_host(host_id, id_or_name="name")
    nfilepaths = len(filepaths)

    # Proceed?
    if not prompt_proceed("proceed removing the simulation files for simulations '" + tostr(filepaths.keys())): continue

    # Remove
    paths = filepaths.values()
    fs.remove_files(paths)

# -----------------------------------------------------------------
