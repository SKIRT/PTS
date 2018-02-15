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
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.launch.analyser import all_steps
from pts.core.config.analyse_simulation import definition as analysis_definition
from pts.core.launch.manager import SimulationManager

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

# For each simulation:
# Get parameter values
# parameter_values = generation.get_parameter_values_for_simulation(simulation.name)
# parameter_string = ", ".join(label + ": " + tostr(value) for label, value in parameter_values.items())

# -----------------------------------------------------------------

# Create simulation manager
manager = SimulationManager()

# Set options
manager.config.analyse = True
manager.config.analysis = config.analysis
manager.config.reanalyse = config.reanalyse
manager.config.features_reanalysis = config.features
manager.config.analyse_simulations = config.simulations
manager.config.prompt_simulations_reanalysis = config.prompt_simulations

# Run the manager
manager.run(assignment=generation.assignment_table)

# -----------------------------------------------------------------
