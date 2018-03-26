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

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import formatting as fmt
from pts.modeling.fitting.statistics import FittingStatistics

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

# Create the fitting statistics
statistics = FittingStatistics()

# Setup
statistics.setup()

# -----------------------------------------------------------------

# Compare generation info
print("")
print(fmt.blue + "GENERATION INFO:" + fmt.reset)
print("")

# Show generation a
statistics.show_generation_info(generation_name_a)

# Show generation b
statistics.show_generation_info(generation_name_b)

# -----------------------------------------------------------------

# Compare simulations
print(fmt.blue + "CLOSEST SIMULATIONS:" + fmt.reset)
print("")

# Show parameters comparison
statistics.show_closest_parameters(generation_name_a, generation_name_b, nsimulations=config.nsimulations)

# -----------------------------------------------------------------

# Plot SEDs
if config.plot_seds: statistics.plot_closest_seds(generation_name_a, generation_name_b, nsimulations=config.nsimulations)

# Plot fluxes
if config.plot_fluxes: statistics.plot_closest_fluxes(generation_name_a, generation_name_b, nsimulations=config.nsimulations)

# -----------------------------------------------------------------
