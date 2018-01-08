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
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.tools import introspection
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("host_ids", "string_list", "names of the remote hosts", choices=find_host_ids())
definition.add_positional_optional("simulation_ids", "integer_list", "simulations to show")
definition.add_optional("names", "string_list", "simulation names to show")

# Add flag
definition.add_flag("analysis", "show analysis options", False)

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("simulations", definition, description="Show the simulations for a certain remote host")

# -----------------------------------------------------------------

# Loop over the hosts
print("")
for host_id in config.host_ids:

    # Show the local analysis runs
    print(fmt.yellow + host_id.upper() + ":" + fmt.reset)
    print("")

    # Loop over the simulations
    for simulation_id in introspection.simulation_ids_for_host(host_id):

        # Check simulation ID
        if config.simulation_ids is not None and simulation_id not in config.simulation_ids: continue

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

# -----------------------------------------------------------------
