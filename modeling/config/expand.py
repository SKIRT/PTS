#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Copy the definition
definition = definition.copy()

# -----------------------------------------------------------------

# Set the modeling path
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Remote hosts
all_host_ids = find_host_ids()
has_remotes = len(all_host_ids) > 0

# -----------------------------------------------------------------

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation name
definition.add_required("generation", "string", "generation name")

# -----------------------------------------------------------------

# Parameters in which to expand
definition.add_required("parameters", "string_list", "parameters for which to expand the range")  # choices are supposed to be the free parameters of a fitting run
definition.add_required("direction", "string_or_string_string_dictionary", "direction in which to expand") #choices=directions)
definition.add_required("npoints", "integer_or_string_integer_dictionary", "number of grid points to add")

# -----------------------------------------------------------------

# Remote or local execution
#if has_remotes: definition.add_positional_optional("remotes", "string_list", "remote hosts to use", default=environment.modeling_configuration.fitting_host_ids, choices=all_host_ids)
#else: definition.add_fixed("remotes", [])
definition.add_positional_optional("host", "host", "remote host to use")
definition.add_flag("local", "run everything locally")

# -----------------------------------------------------------------

# Options
definition.add_flag("attached", "run remote simulations in attached mode")
definition.add_flag("group", "group simulations in larger jobs")
definition.add_optional("walltime", "real", "the preferred walltime per job (for schedulers)")

# -----------------------------------------------------------------

# Update flags
definition.add_flag("update_individuals", "update the individuals table", True)
definition.add_flag("update_parameters", "update the parameters table", True)
definition.add_flag("update_info", "update the generation info", True)
definition.add_flag("update_generations", "update the generations table", True)

# -----------------------------------------------------------------

# Parallelization
definition.add_optional("parallelization", "parallelization", "parallelization scheme for the simulations")
definition.add_optional("nnodes", "positive_integer", "number of computation nodes to use for the simulations")

# -----------------------------------------------------------------
