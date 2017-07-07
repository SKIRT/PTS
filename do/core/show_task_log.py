#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_task_log Show the log output of a remote PTS task.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.tools import filesystem as fs
from pts.core.remote.remote import Remote
from pts.core.basics.task import Task
from pts.core.tools import introspection

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_required("id", "positive_integer", "task ID")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
#setter = InteractiveConfigurationSetter("show_log", "Show the log output of a remote PTS task")
config = parse_arguments("show_task_log", definition, description="Show the log output of a remote PTS task")

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
lines = remote.read_lines(log_path)

# Print the lines of the log file
for line in lines: print(line)

# -----------------------------------------------------------------
