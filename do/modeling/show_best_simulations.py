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
from pts.core.basics.distribution import Distribution
from pts.core.plot.distribution import plot_distributions
from pts.modeling.fitting.refitter import get_and_show_best_simulations
from pts.core.tools import filesystem as fs

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
definition.add_flag("counts", "show the counts")
definition.add_flag("plot_counts", "plot the counts")
definition.add_flag("statistics", "show statistics", True)

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
    parameter_scales[label] = grid_settings[key]

# -----------------------------------------------------------------

# Get the chi squared table
chi_squared = generation.chi_squared_table

# -----------------------------------------------------------------

# Load parameters table
parameters = generation.parameters_table

# -----------------------------------------------------------------

initial_parameter_values = fitting_run.first_guess_parameter_values

# Show initial parameter values
print("")
print("Initial parameter values:")
print("")
for label in parameter_labels: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(initial_parameter_values[label]))

# -----------------------------------------------------------------

# Show
simulation_names, counts = get_and_show_best_simulations(config.nsimulations, parameter_labels, chi_squared, parameters, parameter_units, parameter_scales, initial_parameter_values)
print("")

# -----------------------------------------------------------------

# Show SED?
if config.sed:

    # Loop over the simulations
    for simulation_name in simulation_names:

        # Determine the plot path
        sed_plot_path = generation.get_simulation_sed_plot_path(simulation_name)

        # Show the plot
        fs.open_file(sed_plot_path)

# -----------------------------------------------------------------

# Make counts distributions
counts_distributions = dict()
for label in parameter_labels: counts_distributions[label] = Distribution.from_counts(label, counts[label].values(), counts[label].keys(), sort=True)

# -----------------------------------------------------------------

# Show counts
if config.counts:
    print("Counts in best simulations:")
    for label in parameter_labels:
        print("")
        print(" - " + fmt.bold + label + fmt.reset + ":")
        print("")
        for value in sorted(counts[label].keys()):
            count = counts[label][value]
            relcount = float(count) / config.nsimulations
            print("    * " + tostr(value) + ": " + str(counts[label][value]) + " (" + tostr(relcount*100) + "%)")

# -----------------------------------------------------------------

# Plot the distributions
if config.plot_counts: plot_distributions(counts_distributions, logscale=True, panels=True, frequencies=True)

# -----------------------------------------------------------------

# Show statistics?
if config.statistics:

    # Get best simulation name and chi squared
    best_simulation_name, best_chi_squared = chi_squared.best_simulation_name_and_chi_squared

    # Get best simulation parameter values
    best_parameter_values = parameters.parameter_values_for_simulation(best_simulation_name)

    # Most probable model: should be same as simulation with lowest chi squared
    most_probable_simulation_name = generation.most_probable_model
    #print(most_probable_simulation_name)
    assert most_probable_simulation_name == best_simulation_name

    print("")
    print("Statistics:")

    # Loop over the free parameter labels
    for label in parameter_labels:

        print("")
        print(" - " + fmt.bold + label + fmt.reset + ":")
        print("")

        # Get most probable parameter value
        most_probable_value = generation.get_most_probable_parameter_value(label)

        print("    * Initial guess value: " + tostr(initial_parameter_values[label], ndigits=3))
        print("    * Best simulation value: " + tostr(best_parameter_values[label], ndigits=3))
        print("    * Most probable value: " + tostr(most_probable_value, ndigits=3) + " " + tostr(parameter_units[label]))
        if config.nsimulations > 1: print("    * Most counted in " + str(config.nsimulations) + " best simulations: " + tostr(counts_distributions[label].most_frequent, ndigits=3) + " " + tostr(parameter_units[label]))

# -----------------------------------------------------------------
