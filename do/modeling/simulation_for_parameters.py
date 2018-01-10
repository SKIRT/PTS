#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.simulation_for_parameters Show the simulation from a generation corresponding to specific parameter values.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.basics.configuration import prompt_choices
from pts.core.tools import formatting as fmt
from pts.core.tools import filesystem as fs
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Flags
definition.add_flag("sed", "show the SED")

# Specifiy the indices of the different parameters
definition.add_optional("indices", "integer_list", "indices of the chosen parameter values")

# Create the configuration
config = parse_arguments("simulation_for_parameters", definition)

# -----------------------------------------------------------------

# Get the generation
fitting_run = runs.load(config.name)
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Initialize dictionary for the chosen parameter values
values = dict()

# Loop over the free parameters
for j, label in enumerate(generation.parameter_labels):

    # Get unique values
    unique_values = generation.unique_parameter_values[label]
    unique_values = list(sorted(unique_values))

    # Get value from index
    if config.indices is not None: value = unique_values[config.indices[j]]

    # Prompt for value
    else:

        # Get parameter description
        description = fitting_run.parameter_descriptions[label]

        # Prompt for the value of this parameter
        value = prompt_choices(label, description, unique_values)

    # Set the chosen value
    values[label] = value

# -----------------------------------------------------------------

# Get the simulation name
simulation_name = generation.get_simulation_name_for_parameter_values(values)

# -----------------------------------------------------------------

print("")
print("simulation name: " + fmt.underlined + fmt.green + simulation_name + fmt.reset)
print("")

# -----------------------------------------------------------------

# Show parameter values
for label in generation.parameter_labels: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(values[label]))

# -----------------------------------------------------------------

# Is analysed?
if generation.is_analysed(simulation_name):

    # Get chi squared value
    chisqr = generation.get_chi_squared(simulation_name)

    # Show
    print(" - chi squared: " + tostr(chisqr))

# -----------------------------------------------------------------

# Load the simulation
simulation = generation.get_simulation_or_basic(simulation_name)

# -----------------------------------------------------------------

# Show SED
if config.sed:
    if not generation.has_sed_plot(simulation_name): raise IOError("No SED plot")
    else: fs.open_file(generation.get_simulation_sed_plot_path(simulation_name))

# -----------------------------------------------------------------
