#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_job Show the standard or error output for a particular job script.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.remote.jobscript import JobScript

# -----------------------------------------------------------------

host_ids = find_host_ids(schedulers=True)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Jobscript filepath
definition.add_required("jobscript", "file_path", "name of the jobscript file")

# Host
definition.add_required("host", "host", "remote host", choices=host_ids)

# Output or error
definition.add_required("output_or_error", "string", "show output or error")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("show_job", definition, description="Show the log output of a remote SKIRT simulation")

# -----------------------------------------------------------------

# Load the job script
jobscript = JobScript.from_file(config.jobscript)

# -----------------------------------------------------------------

output_path = jobscript.pbs_options["o"]
error_path = jobscript.pbs_options["e"]

# -----------------------------------------------------------------

# Create and setup the remote
remote = Remote(host_id=config.host)

# -----------------------------------------------------------------

# Get filepath
if config.output_or_error == "error": filepath = error_path
elif config.output_or_error == "output": filepath = output_path
else: raise ValueError("Invalid value for 'output_or_error'")

# -----------------------------------------------------------------

# Check whether the log file exists
if not remote.is_file(filepath): raise RuntimeError("The output file does not (yet) exist remotely")

# Read the log file
lines = remote.read_lines(filepath)

# Print the lines of the log file
for line in lines: print(line)

# -----------------------------------------------------------------
