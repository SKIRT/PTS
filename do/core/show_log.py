#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_log Show the log output of a remote PTS task.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, ArgumentConfigurationSetter
from pts.core.basics.host import find_host_ids
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.remote import Remote
from pts.core.basics.task import Task
from pts.core.tools import introspection

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "the name of the remote host", choices=find_host_ids())
definition.add_required("id", "positive_integer", "the name of the remote host for which to show/retrieve simulations and tasks")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
#setter = InteractiveConfigurationSetter("show_log", "Show the log output of a remote PTS task")
setter = ArgumentConfigurationSetter("show_log", "Show the log output of a remote PTS task")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), time.unique_name("status") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting show_log ...")

# -----------------------------------------------------------------

# Create and setup the remote
remote = Remote()
remote.setup(config.remote)

# Open the task
task_path = fs.join(introspection.pts_run_dir, config.remote, str(config.id) + ".task")
task = Task.from_file(task_path)

# Check whether the log file is present
log_output_path = task.remote_log_path
log_path = None
for filename in remote.files_in_path(log_output_path):
    if "log" in filename:
        log_path = fs.join(log_output_path, filename)
        break

if log_path is None: raise RuntimeError("Log does not exist remotely")

# Read the text file
lines = remote.read_text_file(log_path)

# Print the lines of the log file
for line in lines: print(line)

# -----------------------------------------------------------------
