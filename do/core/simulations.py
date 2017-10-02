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
        print("   - " + fmt.bold + "name: " + fmt.reset + simulation.name)
        print("   - " + fmt.bold + "prefix: " + fmt.reset + simulation.prefix())
        print("   - " + fmt.bold + "ski path: " + fmt.reset + simulation.ski_path)
        print("   - " + fmt.bold + "input path(s): " + fmt.reset + tostr(simulation.input_path))
        print("   - " + fmt.bold + "output path: " + fmt.reset + simulation.output_path)
        print("")

        # More info
        print("   - " + fmt.bold + "remote simulation path: " + fmt.reset + simulation.remote_simulation_path)
        print("   - " + fmt.bold + "remote input path(s): " + fmt.reset + tostr(simulation.remote_input_path))
        print("   - " + fmt.bold + "remote output path: " + fmt.reset + simulation.remote_output_path)
        print("")

        # Show analysis info
        if config.analysis:
            print(simulation.analysis)
            print("")

# -----------------------------------------------------------------
