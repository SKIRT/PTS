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

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.fitting.modelanalyser import SEDFitModelAnalyser
from pts.core.basics.log import log
from pts.core.tools import numbers
from pts.core.basics.distribution import Distribution
from pts.core.plot.distribution import plot_distribution
from pts.core.tools import filesystem as fs

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
nsimulations = len(simulation_names)

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

nzeros = 0

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
    chisq = analyser.chi_squared

    # Get the differences table
    differences = analyser.differences
    differences.sort()

    distribution = Distribution.from_columns("Weighed squared difference", differences.wavelengths(add_unit=False), differences.chi_squared_terms(), unit="micron")
    plot_distribution(distribution, logscale=True, statistics=False)
    exit()

    # Calculate the probability
    probability = np.exp(-0.5 * chisq)
    if probability == 0: nzeros += 1

    # Show chi squared
    log.debug("The chi squared value for simulation '" + simulation_name + "' is " + str(chisq) + " and the probability is " + str(probability))

    # Debugging
    log.debug("Adding to the chi squared table ...")

    # Set the chi squared value
    chi_squared.add_entry(simulation_name, chisq)

    # Save the chi squared table
    chi_squared.save()

    # Debugging
    log.debug("Writing the differences table ...")

    # Set the differences path
    differences_path = generation.get_simulation_misc_differences_path(simulation_name)

    # Save the differences table
    differences.saveto(differences_path)

# -----------------------------------------------------------------

# Show number of zeros
if nzeros > 5: log.warning(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probabilities of zero")
else: log.debug(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probability of zero")

# -----------------------------------------------------------------
