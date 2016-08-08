#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.clean_tasks Clean PTS tasks for a certain remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.host import find_host_ids
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "the name of the remote host for which to clean the tasks", choices=find_host_ids())

# -----------------------------------------------------------------

# Parse the arguments into a configuration
setter = ArgumentConfigurationSetter("clean_tasks", "Clean PTS tasks for a certain remote host")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(config.path, time.unique_name("status") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting clean_tasks ...")

# -----------------------------------------------------------------

# Determine the path to the run directory for the specified remote host
host_run_path = fs.join(introspection.pts_run_dir, config.remote)

# Remove all files in the run directory
for path, name in fs.files_in_path(host_run_path, extension="fits", returns=["path", "name"]):

    # Inform the user
    log.info("Removing task " + name + " ...")

    # Remove the file
    fs.remove_file(path)

# -----------------------------------------------------------------
