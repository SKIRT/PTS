#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.simulations Show the simulations for a certain remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulation_for_host, get_simulation_ids
from pts.core.tools import introspection
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr
from pts.core.tools import filesystem as fs
from pts.core.remote.mounter import mount_remote
from pts.core.remote.remote import get_home_path

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("host_ids", "string_list", "names of the remote hosts", choices=find_host_ids())
definition.add_positional_optional("simulation_ids", "integer_list", "simulations to show")
definition.add_optional("names", "string_list", "simulation names to show")

# Add flag
definition.add_flag("analysis", "show analysis options", False)

# Add flags
definition.add_flag("open_output", "open the output directory (only for one simulation)")
definition.add_flag("open_remote_output", "open the remote output directory (only for one simulation)")
definition.add_flag("open_input", "open the input directory (only for one simulation)")
definition.add_flag("open_remote_input", "open the remote input directory (only for one simulation)")
definition.add_flag("open_base", "open the simulation directory (only for one simulation)")
definition.add_flag("open_remote_base", "open the remote simulation directory (only for one simulation)")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("simulations", definition, description="Show the simulations for a certain remote host")

# -----------------------------------------------------------------

# Check
if config.simulation_ids is not None and len(config.host_ids) > 1: raise ValueError("Cannot specify simulation IDs for more than one host")

# -----------------------------------------------------------------

# More than one simulation
if len(config.host_ids) > 1: multiple_simulations = True
else:
    if config.names is not None: multiple_simulations = len(config.names) > 1
    elif config.simulation_ids is not None: multiple_simulations = len(config.simulation_ids) > 1
    else: multiple_simulations = True

# -----------------------------------------------------------------

# Check flags
if multiple_simulations:
    if config.open_output: raise ValueError("Multiple simulations")
    if config.open_remote_output: raise ValueError("Multiple simulations")
    if config.open_input: raise ValueError("Multiple simulations")
    if config.open_remote_input: raise ValueError("Multiple simulations")
    if config.open_base: raise ValueError("Multiple simulations")
    if config.open_remote_base: raise ValueError("Multiple simulations")

# -----------------------------------------------------------------

# Loop over the hosts
print("")
for host_id in config.host_ids:

    # Show the local analysis runs
    print(fmt.yellow + host_id.upper() + ":" + fmt.reset)
    print("")

    # Get simulation IDS
    if config.names is not None:
        if config.simulation_ids is not None: raise ValueError("Cannot specify both simulation names and IDs")
        simulation_ids = get_simulation_ids(host_id, config.names)
    else: simulation_ids = config.simulation_ids

    # Loop over the simulations
    for simulation_id in introspection.simulation_ids_for_host(host_id):

        # Check simulation ID
        if simulation_ids is not None and simulation_id not in simulation_ids: continue

        # Open the simulation object
        simulation = get_simulation_for_host(host_id, simulation_id)

        # Show the name
        print(" - " + fmt.underlined + fmt.blue + str(simulation_id) + fmt.reset)
        print("")

        # Show basic info
        print(fmt.green + "   BASIC INFO" + fmt.reset)
        print("")
        print("   - " + fmt.bold + "name: " + fmt.reset + simulation.name)
        print("   - " + fmt.bold + "prefix: " + fmt.reset + simulation.prefix())
        print("   - " + fmt.bold + "ski path: " + fmt.reset + simulation.ski_path)
        if simulation.has_input:
            print("   - " + fmt.bold + "input files: " + fmt.reset)
            input = simulation.input
            print("")
            for name in input.names:
                print("      * " + fmt.bold + name + ": " + fmt.reset + input.get_filepath(name))
            print("")
        print("   - " + fmt.bold + "output path: " + fmt.reset + simulation.output_path)
        print("   - " + fmt.bold + "submitted at: " + fmt.reset + simulation.submitted_at)
        print("")

        # More info
        print("   - " + fmt.bold + "remote simulation path: " + fmt.reset + simulation.remote_simulation_path)
        print("   - " + fmt.bold + "remote input path(s): " + fmt.reset + tostr(simulation.remote_input_path))
        print("   - " + fmt.bold + "remote output path: " + fmt.reset + simulation.remote_output_path)
        print("")

        # More info
        print(fmt.green + "   RETRIEVAL" + fmt.reset)
        print("")
        print("   - " + fmt.bold + "remove_remote_input: " + fmt.reset + str(simulation.remove_remote_input))
        print("   - " + fmt.bold + "remove_remote_output: " + fmt.reset + str(simulation.remove_remote_output))
        print("   - " + fmt.bold + "remove_remote_simulation_directory: " + fmt.reset + str(simulation.remove_remote_simulation_directory))
        print("   - " + fmt.bold + "remove_local_output: " + fmt.reset + str(simulation.remove_local_output))
        if simulation.retrieve_types is not None: print("   - " + fmt.bold + "retrieve_types: " + fmt.reset + tostr(simulation.retrieve_types, delimiter=", "))
        print("")

        if simulation.handle is not None:
            print(fmt.green + "   EXECUTION" + fmt.reset)
            print("")
            for line in simulation.handle.to_lines(line_prefix="   - ", fancy=True): print(line)
            print("")

        # Show analysis info
        if config.analysis:
            print(fmt.green + "   ANALYSIS" + fmt.reset)
            print("")
            for line in simulation.analysis.to_lines(line_prefix="  "): print(line)
            print("")

        # Open
        mount_path = None
        remote_home_path = None
        if config.open_output: fs.open_directory(simulation.output_path)
        if config.open_remote_output:
            if mount_path is None: mount_path = mount_remote(host_id)
            if remote_home_path is None: remote_home_path = get_home_path(host_id)
            relative_output_path = fs.relative_to(simulation.remote_output_path, remote_home_path)
            output_path = fs.join(mount_path, relative_output_path)
            fs.open_directory(output_path)
        if config.open_input: fs.open_directory(simulation.input_path)
        if config.open_remote_input:
            if mount_path is None: mount_path = mount_remote(host_id)
            if remote_home_path is None: remote_home_path = get_home_path(host_id)
            relative_input_path = fs.relative_to(simulation.remote_input_path, remote_home_path)
            input_path = fs.join(mount_path, relative_input_path)
            fs.open_directory(input_path)
        if config.open_base: fs.open_directory(simulation.base_path)
        if config.open_remote_base:
            if mount_path is None: mount_path = mount_remote(host_id)
            if remote_home_path is None: remote_home_path = get_home_path(host_id)
            relative_base_path = fs.relative_to(simulation.remote_simulation_path, remote_home_path)
            base_path = fs.join(mount_path, relative_base_path)
            fs.open_directory(base_path)

# -----------------------------------------------------------------
