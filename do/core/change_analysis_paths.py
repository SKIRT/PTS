#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.change_analysis_paths Change certain analysis paths of a single or multiple simulations.

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
from pts.do.core.change_simulation_settings import save_simulations
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_positional_optional("ids", "integer_list", "simulation IDs")

# -----------------------------------------------------------------

definition.add_optional("extraction", "string", "extraction path")
definition.add_optional("plotting", "string", "plotting path")
definition.add_optional("misc", "string", "miscellaneous output path")

# -----------------------------------------------------------------

definition.add_flag("relative", "treat the specified paths as relative to the corresponding simulation directories")
definition.add_flag("create", "create directories that don't yet exist", True)

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("change_analysis_paths", definition, description="Change certain analysis paths")

# -----------------------------------------------------------------

# No IDs specified?
if config.ids is None: config.ids = introspection.simulation_ids_for_host(config.remote)

# -----------------------------------------------------------------

# Load the simulations and put them in a dictionary
simulations = OrderedDict()
for simulation_id in config.ids: simulations[simulation_id] = get_simulation_for_host(config.remote, simulation_id)

# Update analysis options in each simulation
for simulation_id in config.ids: simulations[simulation_id].update_analysis_options()

# Create a dictionary to contain a flag for each simulation that tells whether it has changed
changed = dict()
for simulation_id in simulations: changed[simulation_id] = False

# -----------------------------------------------------------------

# Extraction path
if config.extraction is not None:

    # Debugging
    log.debug("Changing extraction paths ...")

    # Loop over the simulations
    for simulation_id in simulations:

        # Debugging
        log.debug("Changing extraction path of simulation with ID '" + str(simulation_id) + "' ...")

        # Get the simulation
        simulation = simulations[simulation_id]
        name = simulation.name
        base_path = simulation.base_path

        # Get extraction path
        path = config.extraction.replace("$name", name)

        # Determine full path
        if config.relative: path = fs.join(base_path, path)

        # Create if necessary
        if not fs.is_directory(path) and config.create: fs.create_directory(path, recursive=True)

        # Set path
        simulation.analysis.extraction.path = path

        # Changed
        changed[simulation_id] = True

# -----------------------------------------------------------------

# Plotting path
if config.plotting is not None:

    # Debugging
    log.debug("Changing plotting paths ...")

    # Loop over the simulations
    for simulation_id in simulations:

        # Debugging
        log.debug("Changing plotting path of simulation with ID '" + str(simulation_id) + "' ...")

        # Get the simulation
        simulation = simulations[simulation_id]
        name = simulation.name
        base_path = simulation.base_path

        # Get plotting path
        path = config.plotting.replace("$name", name)

        # Determine full path
        if config.relative: path = fs.join(base_path, path)

        # Create if necessary
        if not fs.is_directory(path) and config.create: fs.create_directory(path, recursive=True)

        # Set path
        simulation.analysis.plotting.path = path

        # Changed
        changed[simulation_id] = True

# -----------------------------------------------------------------

# Misc path
if config.misc is not None:

    # Debugging
    log.debug("Changing misc paths ...")

    # Loop over the simulations
    for simulation_id in simulations:

        # Debugging
        log.debug("Changing misc path of simulation with ID '" + str(simulation_id) + "' ...")

        # Get the simulation
        simulation = simulations[simulation_id]
        name = simulation.name
        base_path = simulation.base_path

        # Get plotting path
        path = config.misc.replace("$name", name)

        # Determine full path
        if config.relative: path = fs.join(base_path, path)

        # Create if necessary
        if not fs.is_directory(path) and config.create: fs.create_directory(path, recursive=True)

        # Set path
        simulation.analysis.misc.path = path

        # Changed
        changed[simulation_id] = True

# -----------------------------------------------------------------

# Save the simulations
save_simulations(simulations, changed=changed)

# -----------------------------------------------------------------
