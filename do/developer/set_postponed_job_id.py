#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.set_postponed_job_id Does what it is named

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.simulation.simulation import RemoteSimulation
from pts.core.basics.handle import ExecutionHandle

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "ID of the remote host", choices=find_host_ids())
definition.add_required("simulation_id", "integer", "ID of the simulation")
definition.add_required("job_id", "integer", "job ID")

# Create config
config = parse_arguments("set_postponed_job_id", definition, add_cwd=False)

# -----------------------------------------------------------------

# Determine the path to the run directory for the specified remote host
host_run_path = fs.join(introspection.skirt_run_dir, config.remote)

# -----------------------------------------------------------------

# Determine the path to the simulation file
simulation_file_path = fs.join(host_run_path, config.simulation_id + ".sim")
if not fs.is_file(simulation_file_path): raise ValueError("Simulation file '" + simulation_file_path + "' does not exist")

# -----------------------------------------------------------------

# Open the simulation file
simulation = RemoteSimulation.from_file(simulation_file_path)

# Create execution handle
handle = ExecutionHandle.job(config.job_id, config.host_id)

# Set the handle
simulation.handle = handle

# Save the simulation
simulation.save()

# -----------------------------------------------------------------
