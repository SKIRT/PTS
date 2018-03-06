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
from pts.core.simulation.jobscript import SKIRTJobScript
from pts.core.tools import filesystem as fs
from pts.core.simulation.remote import get_simulation_for_host, get_simulation_id
from pts.core.basics.log import log
from pts.core.basics import containers
from pts.core.tools import formatting as fmt

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
definition.add_flag("from_cwd", "from current working directory")
definition.add_optional("job_ids", "integer_list", "select the jobscripts with these job IDs")
definition.add_optional("names", "string_list", "select the jobscripts with these simulation names")
definition.add_optional("ids", "integer_list", "select the jobscripts with these simulation IDs")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("show_job", definition, description="Show the log output of a remote SKIRT simulation")

# -----------------------------------------------------------------

if config.from_cwd: config.from_directory = fs.cwd()

# -----------------------------------------------------------------

# From directory
if config.from_directory:

    all_jobscripts = dict()

    # Open the files
    for path, name in fs.files_in_cwd(extension="sh", returns=["path", "name"]):
        jobscript = SKIRTJobScript.from_file(path)
        all_jobscripts[name] = jobscript

    # From names
    if config.names is not None:

        #jobscripts = containers.sequence_from_dict(all_jobscripts, config.names)
        jobscripts = containers.create_subdict(all_jobscripts, config.names)

    # From job ID
    elif config.job_ids is not None:

        #jobscripts = []
        jobscripts = dict()

        for simulation_name in all_jobscripts:
            simulation_id = get_simulation_id(config.host.id, simulation_name)
            simulation = get_simulation_for_host(config.host.id, simulation_id)
            #if simulation.handle.value == config.job_id:
            if simulation.handle.value in config.job_ids:
                jobscript = all_jobscripts[simulation_name]
                #break
                #jobscripts.append(jobscript)
                jobscripts[simulation_name] = jobscript

    # From simulation ID
    elif config.ids is not None:

        #jobscripts = []
        jobscripts = dict()

        for simulation_name in all_jobscripts:
            simulation_id = get_simulation_id(config.host.id, simulation_name)
            #if simulation_id == config.simulation_id:
            if simulation_id in config.ids:
                jobscript = all_jobscripts[simulation_name]
                #break
                #jobscripts.append(jobscript)
                jobscripts[simulation_name] = jobscript

    # From nothing?
    else: jobscripts = all_jobscripts #raise ValueError("Not enough input")

# Load
elif config.jobscript is not None:

    jobscripts = dict()
    filename = fs.strip_extension(fs.name(config.jobscript))
    jobscript = SKIRTJobScript.from_file(config.jobscript)
    jobscripts[filename]= jobscript

# Not specified
else: raise ValueError("Jobscript file is not specified")

# -----------------------------------------------------------------

# Create and setup the remote
remote = Remote(host_id=config.host)

# -----------------------------------------------------------------

# Loop over the jobscripts
for name in jobscripts:

    # Get jobscript
    jobscript = jobscripts[name]

    # Get properties
    output_path = jobscript.pbs_options["o"]
    error_path = jobscript.pbs_options["e"]
    skirt_command = jobscript.commands[0][0]
    simulation_output_path = skirt_command.split("-o ")[1].split(" ")[0]
    simulation_prefix = fs.strip_extension(fs.name(skirt_command.split(".ski")[0].split(" ")[-1]))
    logfilepath = fs.join(simulation_output_path, simulation_prefix + "_log.txt")

    print("--------------------------------------------")
    print("Job script '" + name + "'")
    print(" - output path: " + output_path)
    print(" - error path: " + error_path)
    print(" - simulation output path: " + simulation_output_path)
    print(" - simulation prefix: " + simulation_prefix)
    print(" - logfile path: " + logfilepath)
    print("--------------------------------------------")

    # Get filepath
    if config.show == "error": filepath = error_path
    elif config.show == "output": filepath = output_path
    elif config.show == "log": filepath = logfilepath
    else: raise ValueError("Invalid value for 'output_or_error'")

    # Check whether the log file exists
    if not remote.is_file(filepath): raise RuntimeError("The output file does not (yet) exist remotely")

    # Show the filepath
    log.debug("Reading the remote file '" + filepath + "' ...")

    # Read the log file
    lines = remote.read_lines(filepath)

    # Print the lines of the log file
    if config.show == "error":
        for line in lines: print(fmt.red + line + fmt.reset)
    else:
        for line in lines: print(line)

# -----------------------------------------------------------------
