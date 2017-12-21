#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.create_assignment Create simulation assignment file from simulation names.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulations_for_host
from pts.core.launch.batchlauncher import SimulationAssignmentTable

# -----------------------------------------------------------------

host_ids = find_host_ids()

# -----------------------------------------------------------------

# Set the number of allowed open file handles
fs.set_nallowed_open_files(1024)

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_positional_optional("names", "string_list", "simulation names")
definition.add_optional("remotes", "string_list", "remote hosts for which to look for matching simulations", default=host_ids, choices=host_ids)
definition.add_flag("from_directories", "use directory names as simulation names")
definition.add_flag("local", "treat simulations without a match as local simulations", False)
definition.add_flag("success", "success flag to fill in into the assignment table for all simulations", True)
definition.add_flag("show", "show the assignment", False)
config = parse_arguments("create_assignment", definition, "Create simulation assignment file from simulation names")

# -----------------------------------------------------------------

# Set simulation names
if config.from_directories:
    if config.names is not None: raise ValueError("Cannot specify names with 'from_directories' enabled")
    config.names = fs.directories_in_path(returns="name")

# -----------------------------------------------------------------

# Get simulations for each remote host
simulations = dict()
for host_id in config.remotes:

    # Get the simulations
    simulations_host = get_simulations_for_host(host_id, as_dict=True)

    # Set the simulations
    simulations[host_id] = simulations_host

# -----------------------------------------------------------------

# Create assignment
assignment = SimulationAssignmentTable()

# -----------------------------------------------------------------

# Loop over the simulation names and look for matches
for simulation_name in config.names:

    the_host_id = None

    # Loop over the remotes
    for host_id in simulations:
        if simulation_name in simulations[host_id]:
            the_host_id = host_id
            break

    # No remote host for this simulation
    if the_host_id is None:
        if config.local:
            assignment.add_local_simulation(simulation_name, success=config.success)
        else: raise ValueError("Cannot determine the host for simulation '" + simulation_name + "'")

    # Remote host was found
    else:

        # Get the simulation
        simulation = simulations[host_id][simulation_name]

        # Get simulation properties
        simulation_id = simulation.id
        cluster_name = simulation.cluster_name

        # Add to assignment
        assignment.add_remote_simulation(simulation_name, host_id, cluster_name=cluster_name, simulation_id=simulation_id, success=config.success)

# -----------------------------------------------------------------

if config.show: print(assignment)

# -----------------------------------------------------------------

# Write the assignment
assignment.saveto("assignment.dat")

# -----------------------------------------------------------------
