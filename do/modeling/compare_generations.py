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
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr, yes_or_no
from pts.core.plot.sed import plot_seds

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

# Number of best maching simulations
definition.add_optional("nsimulations", "positive_integer", "number of best matching simulations", 5)

# Plot SEDs or fluxes
definition.add_flag("plot_seds", "compare simulated SEDs of matching simulations")
definition.add_flag("plot_fluxes", "compare the mock observed SEDs of matching simulations")
definition.add_flag("plot_datacube_residuals", "plot the residuals of the matching simulation datacubes")
definition.add_flag("plot_image_residuals", "plot the residuals of the matching simulation mock observed images")

# -----------------------------------------------------------------

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

    # Add to the dictionary
    if closest_simulation_name in closest:
        current_diff = closest[closest_simulation_name][1]
        if diff < current_diff: closest[closest_simulation_name] = (simulation_name, diff)
    else: closest[closest_simulation_name] = (simulation_name, diff)

# -----------------------------------------------------------------

# Keep only the closest matches
matches = []

# Sort from best match to worse match
for closest_simulation_name in sorted(closest, key=lambda name: closest[name][1]):

    simulation_name = closest[closest_simulation_name][0]
    diff = closest[closest_simulation_name][1]
    matches.append((simulation_name, closest_simulation_name, diff))
    if len(matches) == config.nsimulations: break

# -----------------------------------------------------------------

# Determine common parameter units
parameter_units = dict()
for label in fitting_run.free_parameter_labels:
    unit = parameters_a.column_unit(label)
    parameter_units[label] = unit

# -----------------------------------------------------------------

print("")
print(fmt.blue + "GENERATION INFO:" + fmt.reset)
print("")

# -----------------------------------------------------------------

def show_generation_info(generation):

    """
    This fnuction ...
    :param generation:
    :return:
    """

    print(" - " + fmt.bold + "Method: " + fmt.reset + generation.method)
    print(" - " + fmt.bold + "Wavelength grid: " + fmt.reset + generation.wavelength_grid_name)
    print(" - " + fmt.bold + "Representation: " + fmt.reset + generation.model_representation_name)
    print(" - " + fmt.bold + "Number of photon packages: " + fmt.reset + yes_or_no(generation.npackages))
    print(" - " + fmt.bold + "Dust self-absorption: " + fmt.reset + yes_or_no(generation.selfabsorption))
    print(" - " + fmt.bold + "Transient heating: " + fmt.reset + yes_or_no(generation.transient_heating))
    print(" - " + fmt.bold + "Spectral convolution: " + fmt.reset + yes_or_no(generation.spectral_convolution))
    print(" - " + fmt.bold + "Use images: " + fmt.reset + yes_or_no(generation.use_images))

# -----------------------------------------------------------------

# Show info of generation a
print(fmt.underlined + fmt.yellow + generation_name_a + fmt.reset)
print("")
show_generation_info(generation_a)
print("")

# Show info of generation b
print(fmt.underlined + fmt.cyan + generation_name_b + fmt.reset)
print("")
show_generation_info(generation_b)
print("")

# -----------------------------------------------------------------

reference_seds = OrderedDict()
reference_seds["Observed clipped fluxes"] = environment.observed_sed
reference_seds["Observed truncated fluxes"] = environment.truncated_sed

# -----------------------------------------------------------------

print(fmt.blue + "CLOSEST SIMULATIONS:" + fmt.reset)
print("")

# Loop over the matches
for index, match in enumerate(matches):

    # Get simulation names
    name_a, name_b = match[0], match[1]

    # Show simulation names
    print(str(index+1) + ". " + fmt.underlined + fmt.green + name_a + fmt.reset + ", " + fmt.underlined + fmt.green + name_b + fmt.reset + ":")
    print("")

    # Get parameter values
    values_a = parameters_a.parameter_values_for_simulation(name_a)
    values_b = parameters_b.parameter_values_for_simulation(name_b)

    # Loop over the parameter values
    for label in fitting_run.free_parameter_labels:

        unit = parameter_units[label]
        value_a = values_a[label].to(unit).value
        value_b = values_b[label].to(unit).value
        reldiff = np.exp(abs(np.log(value_a/value_b))) - 1
        print(" - " + fmt.bold + label + fmt.reset + ": " + fmt.yellow + tostr(value_a, decimal_places=4) + fmt.reset + ", " + fmt.cyan + tostr(value_b, decimal_places=4) + fmt.reset + " " + tostr(unit) + " (" + tostr(reldiff*100, decimal_places=3) + "%)")

    print("")

    # Plot SEDs
    if config.plot_seds:

        # Create SEDs
        seds = OrderedDict()
        #seds.update(reference_seds)
        seds[generation_name_a] = generation_a.get_simulation_sed(name_a)
        seds[generation_name_b] = generation_b.get_simulation_sed(name_b)

        # Plot
        plot_seds(seds)

    # Plot fluxes
    if config.plot_fluxes:

        # Create SEDs
        seds = OrderedDict()
        seds.update(reference_seds)
        seds[generation_name_a] = generation_a.get_simulation_mock_sed(name_a)
        seds[generation_name_b] = generation_b.get_simulation_mock_sed(name_b)

        # Plot
        plot_seds(seds)

# -----------------------------------------------------------------
