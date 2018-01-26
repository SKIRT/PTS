#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_generation_seds Plot the seds of all simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.plot.sed import SEDPlotter
from pts.core.basics.log import log

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run name
if runs.empty: raise ValueError("No fitting runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run for which to (re)make the wavelength grids", choices=runs.names)

# Generation
definition.add_required("generation", "string", "generation name")

# Get configuration
config = parse_arguments("plot_generation_seds", definition, "Plot the seds of all simulations of a generation")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.run)

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Create the SED plotter
plotter = SEDPlotter()

# -----------------------------------------------------------------

# Inform the user
log.info("Adding the observed SEDs ...")

# Add observed SEDs
plotter.add_sed(environment.observed_sed, "Clipped observed fluxes")
plotter.add_sed(environment.truncated_sed, "Truncated observed fluxes")

# -----------------------------------------------------------------

# Add simulation SEDs
log.info("Adding the simulated SEDs ...")

# Loop over the simulation names
nsimulations = generation.nsimulations
for index, simulation_name in enumerate(generation.simulation_names):

    # Debugging
    log.debug("Adding SED of the '" + simulation_name + "' simulation (" + str(index+1) + " of " + str(nsimulations) + ") ...")

    # Get the simulated SED
    sed = generation.get_simulation_sed(simulation_name)

    # Add the SED
    plotter.add_sed(sed, simulation_name, ghost=True)

# -----------------------------------------------------------------

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
