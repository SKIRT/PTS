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
from pts.core.tools import filesystem as fs
from pts.core.simulation.remote import get_simulation_for_host, get_simulation_id

# -----------------------------------------------------------------

host_ids = find_host_ids(schedulers=True)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Output or error
definition.add_required("show", "string", "show output or error or logfile")

# Host
definition.add_required("host", "host", "remote host", choices=host_ids)

# Jobscript filepath
definition.add_positional_optional("jobscript", "file_path", "name of the jobscript file")
definition.add_optional("from_directory", "directory_path", "check the jobscript files in directory")
definition.add_optional("job_id", "positive_integer", "select the jobscript with this job ID")
definition.add_optional("simulation_name", "string", "select the jobscript with this simulation name")
definition.add_optional("simulation_id", "positive_integer", "select the jobscript with this simulation ID")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("show_job", definition, description="Show the log output of a remote SKIRT simulation")

# -----------------------------------------------------------------

# From directory
if config.from_directory:

    jobscripts = dict()

    # Open the files
    for path, name in fs.files_in_cwd(extension="sh", returns=["path", "name"]):
        jobscript = JobScript.from_file(path)
        jobscripts[name] = jobscript

    # From simulation name
    if config.simulation_name is not None: jobscript = jobscripts[config.simulation_name]

    # From job ID
    elif config.job_id is not None:

        for simulation_name in jobscripts:
            simulation_id = get_simulation_id(config.host.id, simulation_name)
            simulation = get_simulation_for_host(config.host.id, simulation_id)
            if simulation.handle.value == config.job_id:
                jobscript = jobscripts[simulation_name]
                break

    # From simulation ID
    elif config.simulation_id is not None:

        for simulation_name in jobscripts:
            simulation_id = get_simulation_id(config.host.id, simulation_name)
            if simulation_id == config.simulation_id:
                jobscript = jobscripts[simulation_name]
                break

    # From nothing?
    else: raise ValueError("Not enough input")

# Load
else: jobscript = JobScript.from_file(config.jobscript)

# -----------------------------------------------------------------

output_path = jobscript.pbs_options["o"]
error_path = jobscript.pbs_options["e"]
skirt_command = jobscript.commands[0][0]
simulation_output_path = skirt_command.split("-o ")[1].split(" ")[0]
simulation_prefix = fs.strip_extension(fs.name(skirt_command.split(".ski")[0].split(" ")[-1]))
logfilepath = fs.join(simulation_output_path, simulation_prefix + "_log.txt")

# -----------------------------------------------------------------

# Create and setup the remote
remote = Remote(host_id=config.host)

# -----------------------------------------------------------------

# Get filepath
if config.show == "error": filepath = error_path
elif config.show == "output": filepath = output_path
elif config.show == "log": filepath = logfilepath
else: raise ValueError("Invalid value for 'output_or_error'")

# -----------------------------------------------------------------

# Check whether the log file exists
if not remote.is_file(filepath): raise RuntimeError("The output file does not (yet) exist remotely")

# Read the log file
lines = remote.read_lines(filepath)

# Print the lines of the log file
for line in lines: print(line)

# -----------------------------------------------------------------
