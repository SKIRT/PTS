#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_best_simulations Show the best fitting simulation from a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr
from pts.core.tools import filesystem as fs
from pts.core.tools import nr
from pts.core.tools import numbers

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

# Rank
definition.add_positional_optional("nsimulations", "positive_integer", "number of simulations to show", 5)

# Options
definition.add_flag("sed", "show the SED plot")

# Create the configuration
config = parse_arguments("show_best_simulations", definition)

# -----------------------------------------------------------------

# Load the fitting run and the generation
fitting_run = runs.load(config.name)
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Get parameter ranges
parameter_ranges = fitting_run.free_parameter_ranges
parameter_labels = parameter_ranges.keys()

# Get parameter units
parameter_units = fitting_run.parameter_units

# Show ranges
print("")
print("Parameter ranges:")
print("")
for label in parameter_labels: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(parameter_ranges[label]))

# -----------------------------------------------------------------

# Get the scales
grid_settings = fitting_run.grid_settings
parameter_scales = dict()
for label in parameter_labels:
    key = label + "_scale"
    parameter_scales[label]= grid_settings[key]

# -----------------------------------------------------------------

# Get the chi squared table
chi_squared = generation.chi_squared_table

# -----------------------------------------------------------------

# Load parameters table
parameters = generation.parameters_table

# -----------------------------------------------------------------

unique_values = parameters.unique_parameter_values
unique_values_scalar = dict()
for label in unique_values:
    values = list(sorted([value.to(parameter_units[label]).value for value in unique_values[label]]))
    unique_values_scalar[label] = values
#print(unique_values)
#print(unique_values_scalar)

# -----------------------------------------------------------------

initial_parameter_values = fitting_run.first_guess_parameter_values
#print(initial_parameter_values)

# Show initial parameter values
print("")
print("Initial parameter values:")
print("")
for label in parameter_labels: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(initial_parameter_values[label]))

# -----------------------------------------------------------------

# Get best simulation name
simulation_names = chi_squared.get_best_simulation_names(config.nsimulations)

# Loop over the simulations
for index, simulation_name in enumerate(simulation_names):

    # Get the chi squared value
    chisq = chi_squared.chi_squared_for(simulation_name)

    # Get parameter values
    parameter_values = parameters.parameter_values_for_simulation(simulation_name)

    # Show
    print("")
    print(" " + str(index+1) + ". " + fmt.green + fmt.underlined + simulation_name + fmt.reset)
    print("")

    print("  - chi squared: " + str(chisq))
    print("  - parameters:")
    for label in parameter_values:

        unit = parameter_units[label]
        initial_value = initial_parameter_values[label]
        initial_value_scalar = initial_value.to(unit).value
        value = parameter_values[label]
        value_scalar = value.to(unit).value
        nunique_values = len(unique_values_scalar[label])

        # Get index of value and index of initial value
        i = unique_values_scalar[label].index(value_scalar)
        j = nr.locate_continuous(unique_values_scalar[label], initial_value_scalar, scale=parameter_scales[label])
        if not numbers.is_integer(j, absolute=False): raise NotImplementedError("Not yet implemented")
        j = int(round(j))

        # Show
        #indicator = "| "
        indicator = "[ "
        for _ in range(nunique_values):
            if _ == j: character = "+"
            else: character = "o"
            if _ == i: indicator += fmt.red + character + fmt.reset + " "
            else: indicator += character + " "
        #indicator += "|"
        indicator += "]"

        print("     * " + fmt.bold + label + fmt.reset + ": " + tostr(value) + "   " + indicator)

    # Show SED?
    if config.sed:

        # Determine simulation path
        simulation_path = generation.get_simulation_path(simulation_name)

        # Simulation paths
        output_path = fs.join(simulation_path, "out")
        extract_path = fs.join(simulation_path, "extr")
        plot_path = fs.join(simulation_path, "plot")
        misc_path = fs.join(simulation_path, "misc")

        # Determine SED path
        sed_plot_path = fs.join(plot_path, "sed.pdf")

        # Show the plot
        fs.open_file(sed_plot_path)

print("")

# -----------------------------------------------------------------
