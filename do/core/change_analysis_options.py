#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.change_analysis_options Change certain analysis options for a single or multiple simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.tools import introspection
from pts.do.core.change_simulation_settings import save_simulation, save_simulations

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_positional_optional("id", "positive_integer", "simulation ID")
definition.add_positional_optional("matching", "string", "only adapt settings with a name matching this string", suggestions=["remote"])

# -----------------------------------------------------------------

# Select certain properties
definition.add_optional("contains", "string", "only adapt properties containing this string in their name")
definition.add_optional("not_contains", "string", "don't adapt properties containing this string in their name")
definition.add_optional("exact_name", "string", "only adapt properties with this exact string as their name")
definition.add_optional("exact_not_name", "string", "don't adapt properties with this exact string as their name")
definition.add_optional("startswith", "string", "only adapt properties whose name starts with this string")
definition.add_optional("endswith", "string", "only adapt properties whose name starts with this string")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("change_analysis_options", definition, description="Change certain analysis options for a simulation")

# -----------------------------------------------------------------

# Check settings
if config.matching is not None:
    if config.contains is not None: raise ValueError("Cannot specify both matching string and containing string")
    config.contains = config.matching

# -----------------------------------------------------------------

# No IDs specified?
if config.ids is None: config.ids = introspection.simulation_ids_for_host(config.remote)

# -----------------------------------------------------------------

nids = len(config.ids)
has_single_id = nids == 1

# -----------------------------------------------------------------

# Adapt a single simulation
if has_single_id:

    # Open the simulation object
    single_id = config.ids[0]
    simulation = get_simulation_for_host(config.remote, single_id)

    # Update
    simulation.update_analysis_options()

    # Check whether analysis options are defined
    simulation.analysis.prompt_properties(contains=config.contains, not_contains=config.not_contains, exact_name=config.exact_name, exact_not_name=config.exact_not_name, startswith=config.startswith, endswith=config.endswith)

    # Save the simulation
    save_simulation(simulation)

# -----------------------------------------------------------------

# Adapt multiple simulations
else:

    # Load the simulations and put them in a dictionary
    simulations = OrderedDict()
    for simulation_id in config.ids: simulations[simulation_id] = get_simulation_for_host(config.remote, simulation_id)

    # Update analysis options in each simulation
    for simulation_id in config.ids: simulations[simulation_id].update_analysis_options()

    # Create a dictionary to contain a flag for each simulation that tells whether it has changed
    changed = dict()
    for simulation_id in simulations: changed[simulation_id] = False

    # Loop over the properties
    #for name in properties:
        #pass

    raise NotImplementedError("Not yet implemented")

    # Save the simulations
    save_simulations(simulations, changed=changed)

# -----------------------------------------------------------------
