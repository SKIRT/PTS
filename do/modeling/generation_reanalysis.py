#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_reanalysis Reanalyse simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.launch.analyser import reanalyse_simulation, steps, batch
from pts.core.basics.log import log
from pts.core.config.analyse_simulation import definition as analysis_definition

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

all_steps = steps + [batch]

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Simulations to reanalyse
definition.add_optional("simulations", "string_list", "simulation names")

# Reanalyse which steps?
definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=all_steps, default=all_steps)
definition.add_positional_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")
definition.add_optional("not_steps", "string_list", "don't analyse these steps", choices=all_steps)
definition.add_optional("not_features", "string_list", "don't analyse these features (if a single not_step is defined)")

# Add section for analysis options
definition.import_section("analysis", "analyser options", analysis_definition)

# Create the configuration
config = parse_arguments("generation_reanalysis", definition)

# -----------------------------------------------------------------

# Load the fitting run and the generation
fitting_run = runs.load(config.name)
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Loop over the simulations
for simulation in generation.simulations:

    # Check
    if config.simulations is not None and simulation.name not in config.simulations: continue

    # Inform the user
    log.info("Re-analysing simulation '" + simulation.name + "' ...")

    # Reanalyse the simulation
    reanalyse_simulation(simulation, config.steps, config.features, not_steps=config.not_steps, not_features=config.not_features, config=config.analysis)

# -----------------------------------------------------------------
