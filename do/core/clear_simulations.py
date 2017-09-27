#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.clear Clear SKIRT simulations for a certain remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.simulation.simulation import RemoteSimulation
from pts.core.remote.remote import Remote
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add settings
definition.add_positional_optional("remotes", "string_list", "the IDs of the remote hosts for which to clear the simulations", choices=find_host_ids(), default=find_host_ids())
definition.add_positional_optional("ids", "integer_list", "the IDs of the simulations to clear")
definition.add_flag("full", "fully clear the simulations, also remove remote simulation directories")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("clear_tasks", definition, description="Clear PTS tasks for a certain remote host")

# -----------------------------------------------------------------

# Loop over the remote hosts
for host_id in config.remotes:

    # Check whether the remote is available
    if config.full:
        remote = Remote()
        if not remote.setup(host_id):
            log.warning("The remote host '" + host_id + "' is not available: skipping ...")
            continue
    else: remote = None

    # Determine the path to the run directory for the specified remote host
    host_run_path = fs.join(introspection.skirt_run_dir, host_id)

    # Check if there are simulations
    if not fs.is_directory(host_run_path):
        log.debug("No run directory for host '" + host_id + "'")
        continue
    if fs.is_empty(host_run_path): log.debug("No simulations for host '" + host_id + "'")

    # Loop over the simulation files in the run directory
    for path, name in fs.files_in_path(host_run_path, extension="sim", returns=["path", "name"], sort=int):

        # Skip
        if config.ids is not None and int(name) not in config.ids: continue

        # Inform the user
        log.info("Removing simulation " + name + " ...")

        # Fully clear
        if config.full:

            # Debugging
            log.debug("Removing the remote files and directories ...")

            # Load the simulation
            simulation = RemoteSimulation.from_file(path)

            # Remove the simulation from the remote
            simulation.remove_from_remote(remote, full=True)

        # Debugging
        log.debug("Removing the simulation file ...")

        # Remove the file
        fs.remove_file(path)

    # Success
    log.success("All cleared for host '" + host_id + "'")

# -----------------------------------------------------------------
