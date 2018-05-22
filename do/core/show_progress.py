#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_progress Show the progress of a remotely running simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import all_host_ids
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.simulation.status import LogSimulationStatus
from pts.core.remote.remote import Remote
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=all_host_ids())
definition.add_required("id", "positive_integer", "simulation ID")
definition.add_flag("debug_output", "show all simulation output in debug mode")

# Read the command line arguments
config = parse_arguments("show_progress", definition, description="Show the progress of a remotely running simulation")

# -----------------------------------------------------------------

# Load the simulation
simulation = get_simulation_for_host(config.remote, config.id)

# -----------------------------------------------------------------

# Load the remote
remote = Remote()
if not remote.setup(host_id=config.remote): raise RuntimeError("Could not connect to the remote host")

# -----------------------------------------------------------------

# Create the status object
status = LogSimulationStatus(simulation.remote_log_file_path, remote=remote, debug_output=config.debug_output)

# Show the simulation progress
with log.no_debugging(): success = status.show_progress(simulation.handle)

# Check whether not crashed
if not success: raise RuntimeError("The simulation crashed")

# -----------------------------------------------------------------
