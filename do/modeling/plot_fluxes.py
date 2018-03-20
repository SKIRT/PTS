#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_fluxes Plot the mock observed fluxes from the simulations of a generation together with the observed fluxes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.plot.sed import SEDPlotter
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Set the modeling path
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to plot the fluxes
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Create the configuration
config = parse_arguments("plot_fluxes", definition, "Plot the mock observed fluxes from the simulations of a generation together with the observed fluxes")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Get the clipped and truncated observed SED
clipped = environment.observed_sed
truncated = environment.truncated_sed

# -----------------------------------------------------------------

# Loop over the simulations
for simulation_name in generation.simulation_names:

    # No fluxes?
    #if not generation.has_misc_fluxes(simulation_name): continue
    if not generation.has_mock_sed(simulation_name): continue

    # Debugging
    log.debug("Plotting fluxes for the '" + simulation_name + "' simulation ...")

    # Get the misc/fluxes directory path
    fluxes_path = generation.get_simulation_misc_fluxes_path(simulation_name)

    # Get the mock sed
    sed = generation.get_mock_sed(simulation_name)

    # Determine the path to the plot file
    path = fs.join(fluxes_path, "earth_fluxes.pdf")

    # Initialize the plotter
    plotter = SEDPlotter()

    # Add the SED
    plotter.add_sed(sed, "Mock observation")

    # Get the SED
    simulated_sed = generation.get_simulation_sed(simulation_name)

    # Add
    plotter.add_sed(simulated_sed, "Simulation")

    # Add references
    plotter.add_sed(clipped, "Observed clipped fluxes")
    plotter.add_sed(truncated, "Observed truncated fluxes")

    # Run the plotter
    plotter.run(output=path)

# -----------------------------------------------------------------
