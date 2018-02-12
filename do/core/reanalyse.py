#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.reanalyse Re-analyse certain simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import all_host_ids
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.launch.analyser import reanalyse_simulation, all_steps
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=all_host_ids())
definition.add_required("ids", "integer_list", "simulation IDs")
definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=all_steps, default=all_steps)
definition.add_positional_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")
definition.add_optional("not_steps", "string_list", "don't analyse these steps", choices=all_steps)
definition.add_optional("not_features", "string_list", "don't analyse these features (if a single not_step is defined)")

# Read the command line arguments
config = parse_arguments("reanalyse", definition, description="Re-analyse a certain simulation")

# -----------------------------------------------------------------

# Loop over the simulations
for simulation_id in config.ids:

    # Load the simulation
    simulation = get_simulation_for_host(config.remote, config.id)

    # Check whether retrieved (specific to remote simulations)
    if not simulation.retrieved: raise ValueError("The simulation has not been retrieved yet")

    # Inform the user
    log.info("Re-analysing simulation '" + simulation.name + "' (" + config.remote + " " + str(config.id) + ") ...")

    # Reanalyse simulation
    reanalyse_simulation(simulation, config.steps, config.features, not_steps=config.not_steps, not_features=config.not_features)

# -----------------------------------------------------------------
