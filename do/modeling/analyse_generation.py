#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.analyse_generation Analyse the simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_proceed
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.launch.analyser import reanalyse_simulation, all_steps, analyse_simulation
from pts.core.config.analyse_simulation import definition as analysis_definition
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

# Load fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# The fitting run name
if runs.empty: raise ValueError("No fitting runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation
definition.add_required("generation", "string", "generation name")

# Simulation names
definition.add_positional_optional("simulations", "string_list", "simulation names")
definition.add_flag("prompt_simulations", "prompt before analysing a particular simulation")

# Re-analyse?
definition.add_optional("reanalyse", "string_list", "apply re-analysis of these steps", choices=all_steps)
definition.add_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")

# Add section for analysis options
definition.import_section("analysis", "analyser options", analysis_definition)

# -----------------------------------------------------------------

# Create the configuration
config = parse_arguments("analyse_generation", definition, description="Analyse the simulations of a generation")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.run)

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Check
if not generation.has_assignment_table: raise RuntimeError("No assignment for this generation")

# -----------------------------------------------------------------

# Get the simulations
if config.simulations is not None: simulations = [generation.get_simulation(name) for name in config.simulations]
else: simulations = generation.retrieved_simulations

# -----------------------------------------------------------------

# Loop over the retrieved simulations
for simulation in simulations:

    # Get parameter values
    parameter_values = generation.get_parameter_values_for_simulation(simulation.name)
    parameter_string = ", ".join(label + ": " + tostr(value) for label, value in parameter_values.items())

    # Analyse?
    if config.prompt_simulations and not prompt_proceed("analyse simulation '" + simulation.name + "' [" + parameter_string + "] ?"): continue

    # Re-analysis
    if config.reanalyse is not None:

        # Inform the user
        log.info("Re-analysing simulation '" + simulation.name + "' (" + simulation.host_id + " " + str(simulation.id) + ") ...")

        # Reanalyse simulation
        reanalyse_simulation(simulation, config.reanalyse, config.features)

    # Normal analysis
    else:

        # Inform the user
        log.info("Analysing simulation '" + simulation.name + "' (" + simulation.host_id + " " + str(simulation.id) + ") ...")

        # Analyse the simulation
        analyse_simulation(simulation, config=config.analysis)

# -----------------------------------------------------------------
