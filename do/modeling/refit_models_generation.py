#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.refit_models_generation Refit the models (calculate the chis squared) for the simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.fitting.modelanalyser import SEDFitModelAnalyser
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Load the modeling environment and the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Fitting run
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation
definition.add_required("generation", "string", "name of the generation")

# Simulation names
definition.add_positional_optional("simulations", "string_list", "names of the simulations to refit")

# Create the configuration
config = parse_arguments("refit_models_generation", definition, "Refit the models (calculate the chis squared) for the simulations of a generation")

# -----------------------------------------------------------------

# Load the fitting run and the generation
fitting_run = runs.load(config.run)
generation = fitting_run.get_generation(config.generation)

# Set the simulation names
if config.simulations is not None: simulation_names = config.simulations
else: simulation_names = generation.simulation_names

# -----------------------------------------------------------------

# Open the chi squared table
chi_squared = generation.chi_squared_table

# -----------------------------------------------------------------

# Loop over the simulations, remove the chi squared entry
for simulation_name in simulation_names:

    # Does simulation exist in table?
    if not chi_squared.has_simulation(simulation_name): continue

    # Remove entry from table
    #chi_squared.remove_simulation(simulation_name)

# -----------------------------------------------------------------

# Save the chi squared table
chi_squared.save()

# -----------------------------------------------------------------

# Loop over the simulations, run the SED fit model analyser
for simulation_name in simulation_names:

    # Debugging
    log.debug("Refitting the '" + simulation_name + "' simulation ...")

    # Get the simulation object
    simulation = generation.get_simulation_or_basic(simulation_name)

    # Get the mock sed
    mock_sed = generation.get_simulation_mock_sed(simulation_name)

    # Create the fit model analyser
    analyser = SEDFitModelAnalyser()

    # Set options
    analyser.config.update_generation = False
    analyser.config.write = False

    # Run the analyser
    analyser.run(simulation=simulation, fitting_run=fitting_run, mock_sed=mock_sed)

    # Get the chi squared
    chi_squared = analyser.chi_squared

    # Get the differences table
    differences = analyser.differences

    # Show chi squared
    print(chi_squared)

# -----------------------------------------------------------------
