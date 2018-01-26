#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.move_simulations Move the simulation objects from the SKIRT run directory to another directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulation_paths_for_host
from pts.core.basics.log import log

# -----------------------------------------------------------------

host_ids = find_host_ids()

# -----------------------------------------------------------------

# Set the number of allowed open file handles
#fs.set_nallowed_open_files(1024)

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("names", "string_list", "simulation names")
definition.add_optional("remotes", "string_list", "remote hosts for which to look for matching simulations", default=host_ids, choices=host_ids)
definition.add_flag("from_directories", "use directory names as simulation names")
definition.add_optional("not_contains", "string_list", "ignore directories containing this string in their name")
definition.add_optional("exact_not_name", "string_list", "ignore directories with these names")
definition.add_optional("output", "directory_path", "output directory")
definition.add_flag("per_host", "put the simulations in a separate subdirectory for each host")
definition.add_flag("rename", "rename the simulation files to have the name of the simulation")
definition.add_flag("ignore_missing", "ignore missing simulations (otherwise an error is thrown)")

# Create the configuration
config = parse_arguments("move_simulations", definition, "Move the simulation objects from the SKIRT run directory to another directory")

# -----------------------------------------------------------------

# Set simulation names
if config.from_directories:
    if config.names is not None: raise ValueError("Cannot specify names with 'from_directories' enabled")
    config.names = fs.directories_in_path(returns="name", not_contains=config.not_contains, exact_not_name=config.exact_not_name)

# -----------------------------------------------------------------

# Get simulations for each remote host
simulation_paths = dict()
for host_id in config.remotes:

    # Get the simulation paths
    paths_host = get_simulation_paths_for_host(host_id, as_dict=True)

    # Set the simulation paths
    simulation_paths[host_id] = paths_host

# -----------------------------------------------------------------

# Determine the output directory
output_path = config.output if config.output is not None else fs.cwd()

# -----------------------------------------------------------------

# Loop over the simulation names and look for matches
for simulation_name in config.names:

    the_host_id = None

    # Loop over the remotes, look for match
    for host_id in simulation_paths:
        if simulation_name in simulation_paths[host_id]:
            the_host_id = host_id
            break

    # Simulation file not found
    if the_host_id is None:
        if config.ignore_missing: continue
        else: raise ValueError("Cannot find simulation file for simulation '" + simulation_name + "'")

    # Determine the output path
    if config.per_host:
        new_path = fs.join(output_path, the_host_id)
        if not fs.is_directory(new_path): fs.create_directory(new_path)
    else: new_path = output_path

    # Get the original simulation file path
    filepath = simulation_paths[the_host_id][simulation_name]

    # Debugging
    log.debug("Moving the '" + simulation_name + "' simulation ...")

    # Move the file
    if config.rename: new_name = simulation_name + ".sim"
    else: new_name = None
    fs.move_file(filepath, new_path, new_name=new_name)

# -----------------------------------------------------------------
