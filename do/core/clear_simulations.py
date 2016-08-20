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
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.host import find_host_ids
from pts.core.tools import logging, introspection
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "the name of the remote host for which to clear the simulations", choices=find_host_ids())

# Add optional
definition.add_positional_optional("ids", "integer_list", "the IDs of the simulations to clear")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
setter = ArgumentConfigurationSetter("clear_tasks", "Clear PTS tasks for a certain remote host")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level)
log.start("Starting clear_simulations ...")

# -----------------------------------------------------------------

# Determine the path to the run directory for the specified remote host
host_run_path = fs.join(introspection.skirt_run_dir, config.remote)

# Loop over the simulation files in the run directory
for path, name in fs.files_in_path(host_run_path, extension="sim", returns=["path", "name"]):

    # Skip
    if config.ids is not None and int(name) not in config.ids: continue

    # Inform the user
    log.info("Removing simulation " + name + " ...")

    # Remove the file
    fs.remove_file(path)

# -----------------------------------------------------------------
