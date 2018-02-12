#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.set_execution_handles Set the execution handles of simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_string, prompt_integer
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.remote.host import find_host_ids
from pts.core.tools import introspection
from pts.core.simulation.remote import get_simulations_for_host, get_simulation_for_host
from pts.core.basics.handle import ExecutionHandle, handle_types
from pts.core.launch.options import LoggingOptions, SchedulingOptions
from pts.core.tools import sequences

# -----------------------------------------------------------------

host_ids = find_host_ids(schedulers=True)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("host", "host", "name of the remote host", choices=host_ids)
definition.add_positional_optional("ids", "integer_list", "simulation IDs")
definition.add_optional("names", "string_list", "simulation names")
definition.add_flag("from_directories", "use directory names as simulation names")

# Type of execution handle
definition.add_optional("type", "string", "type of handles", choices=handle_types)

# From qstat output
definition.add_optional("from_qstat", "file_path", "from qstat output (in text file)")

# Logging options
definition.import_section_from_composite_class("logging", "logging options", LoggingOptions)

# Read the command line arguments
config = parse_arguments("set_execution_handles", definition, description="Set the execution handles of simulations")

# -----------------------------------------------------------------

# Set simulation names
if config.from_directories:
    if config.names is not None: raise ValueError("Cannot specify names with 'from_directories' enabled")
    config.names = fs.directories_in_path(returns="name")

# -----------------------------------------------------------------

# Initialize list for the simulations
simulations = []

# Set simulation IDs
if config.names is None and config.ids is None: config.ids = introspection.simulation_ids_for_host(config.host.id)

# From names
if config.names is not None:

    # Get simulations for host
    all_simulations = get_simulations_for_host(config.host.id, as_dict=True)

    # Loop over the names
    for name in config.names:
        simulation_id = all_simulations[name].id
        simulations.append(all_simulations[name])

# From IDS
elif config.ids is not None:

    # Load the simulations and put them in the dictionary
    for simulation_id in config.ids: simulations.append(get_simulation_for_host(config.host.id, simulation_id))

# Nothing
else: raise ValueError("Names or IDs must be specified")

# -----------------------------------------------------------------

simulation_names = [simulation.name for simulation in simulations]

# -----------------------------------------------------------------

if config.from_qstat:

    # EXAMPLE:
    # 4034022.master15.delcatty.gen ...-59-14-479486
    # 4034023.master15.delcatty.gen ...-59-14-470028
    # 4034024.master15.delcatty.gen ...-59-14-464181
    # 4034025.master15.delcatty.gen ...-59-14-456412
    # 4034026.master15.delcatty.gen ...-59-14-447609
    # 4034027.master15.delcatty.gen ...-59-14-441748
    # 4034028.master15.delcatty.gen ...-59-14-435743
    # ...

    job_ids = dict()

    for line in fs.read_lines(config.from_qstat):

        job_id = int(line.split(".")[0])
        end = line.split("...")[1]

        #print(job_id, end)

        # Find the simulation name
        simulation_name = sequences.find_unique_endswith(simulation_names, end)

        # Set the job ID
        job_ids[simulation_name] = job_id

else: job_ids = None

# -----------------------------------------------------------------

# Set the handles, save the simulation objects
for simulation in simulations:

    # Show simulation name
    print("Simulation: " + simulation.name)

    if job_ids is not None:

        if simulation.name not in job_ids: continue

        job_id = job_ids[simulation.name]
        handle = ExecutionHandle.job(job_id, config.host.id)

    else:

        # Create handle
        if config.type is None: config.type = prompt_string("handle_type", "type of execution handle", choices=handle_types)

        # Screen
        if config.type == "screen":

            screen_name = prompt_string("screen_name", "screen session name")
            handle = ExecutionHandle.screen(screen_name, config.host.id)

        # Job
        elif config.type == "job":

            job_id = prompt_integer("job_id", "job ID")
            handle = ExecutionHandle.job(job_id, config.host.id)

        # Not implemented
        else: raise NotImplementedError("Not yet implemented")

    #print(handle)
    #continue

    # Set the handle
    simulation.handle = handle

    # Save the simulation
    simulation.save()

# -----------------------------------------------------------------
