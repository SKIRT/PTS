#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.compare_generations Compare simulation parameters between generations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd

# -----------------------------------------------------------------

# Load the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generation names
definition.add_positional_optional("generations", "string_pair", "names of generations to compare")

# Create the configuration
config = parse_arguments("compare_generations", definition, "Compare simulation parameters between generations")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# -----------------------------------------------------------------

# Get names of generations
if config.generations is None:

    # Test number of generations
    if fitting_run.ngenerations != 2: raise RuntimeError("Number of generations is different from 2")

    # Get the two generation names
    generation_name_a, generation_name_b = fitting_run.generation_names

# Get the specified generation names
else: generation_name_a, generation_name_b = config.generations

# -----------------------------------------------------------------

# Load the generations
generation_a = fitting_run.get_generation(generation_name_a)
generation_b = fitting_run.get_generation(generation_name_b)

# Load the generation parameters tables
parameters_a = generation_a.parameters_table
parameters_b = generation_b.parameters_table

# -----------------------------------------------------------------

closest = dict()

# -----------------------------------------------------------------

# Loop over the simulations of one generation
for index in range(len(parameters_a)):

    # Get the simulation name
    simulation_name = parameters_a.get_simulation_name(index)

    # Get parameter values
    values = parameters_a.get_parameter_values(index)

    # Find the simulation with the closest parameter values for the other generation
    closest_index = parameters_b.closest_simulation_for_parameter_values(values, return_index=True)
    closest_simulation_name = parameters_b.get_simulation_name(closest_index)

    # Get the closest values
    closest_values = parameters_b.get_parameter_values(closest_index)

    # Loop over the labels
    diff_keys = []
    for label in fitting_run.free_parameter_labels:

        # Get both values
        value = values[label]
        closest_value = closest_values[label]

        # Calculate the difference
        diff = np.exp(abs(np.log(value / closest_value)))

        # Add the difference
        diff_keys.append(diff)

    # Calculate difference
    diff = np.sqrt(np.sum(np.power(diff_keys, 2)))
    #print(diff)

    if closest_simulation_name in closest:
        current_diff = closest[closest_simulation_name][1]
        if diff < current_diff: closest[closest_simulation_name] = (simulation_name, diff)
    else: closest[closest_simulation_name] = (simulation_name, diff)

# -----------------------------------------------------------------

matches = []
for closest_simulation_name in sorted(closest, key=lambda name: closest[name][1]):
    simulation_name = closest[closest_simulation_name][0]
    diff = closest[closest_simulation_name][1]
    matches.append((simulation_name, closest_simulation_name, diff))

# -----------------------------------------------------------------

for match in matches: print(match)

# -----------------------------------------------------------------
