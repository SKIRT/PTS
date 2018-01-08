#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_simulation_log Show the log output of a remote SKIRT simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.tools import filesystem as fs
from pts.core.remote.remote import Remote
from pts.core.simulation.remote import get_simulation_for_host, get_simulation_id

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_positional_optional("id", "positive_integer", "simulation ID")
definition.add_optional("name", "string", "simulation name")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("show_simulation_log", definition, description="Show the log output of a remote SKIRT simulation")

# -----------------------------------------------------------------

# Determine simulation ID
if config.name is not None:
    if config.id is not None: raise ValueError("Cannot specifiy both name and simulation ID")
    simulation_id = get_simulation_id(config.remote, config.name)
else: simulation_id = config.id

# -----------------------------------------------------------------

# Open the simulation
simulation = get_simulation_for_host(config.remote, simulation_id)

# The name of the ski file (the simulation prefix)
ski_name = simulation.prefix()

# Simulation is retrieved
if simulation.retrieved:

    # Determine the path to the simulation log file
    local_log_file_path = simulation.log_file_path

    # Read the log file
    lines = fs.read_lines(local_log_file_path)

# Not yet retrieved
else:

    # The path to the simulation log file
    remote_log_file_path = simulation.remote_log_file_path

    # Create and setup the remote
    remote = Remote()
    remote.setup(config.remote)

    # Check whether the log file exists
    if not remote.is_file(remote_log_file_path): raise RuntimeError("The log file does not (yet) exist remotely")

    # Read the log file
    lines = remote.read_lines(remote_log_file_path)

# Print the lines of the log file
for line in lines: print(line)

# -----------------------------------------------------------------
