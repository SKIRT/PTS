#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.relaunch_simulations_scheduler Relaunch simulations as jobs on a scheduling system.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.remote import SKIRTRemote
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.remote.host import find_host_ids
from pts.core.tools import introspection
from pts.core.simulation.remote import get_simulations_for_host, get_simulation_for_host
from pts.core.simulation.remote import is_finished_status
from pts.core.basics.handle import ExecutionHandle
from pts.core.launch.options import LoggingOptions, SchedulingOptions

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

# Flags
definition.add_flag("finished", "relaunch already finished simulations", True)
definition.add_flag("dry", "dry run")
definition.add_flag("mail", "send mail when jobs start")

# Walltime
definition.add_optional("walltime", "duration", "walltime")

# Jobscripts path
definition.add_optional("jobscripts_path", "directory_path", "path of the directory for the jobscripts to be saved locally")

# Parallelization
definition.add_optional("parallelization", "parallelization", "parallelization scheme for the simulations")

# Logging options
definition.import_section_from_composite_class("logging", "logging options", LoggingOptions)

# Read the command line arguments
config = parse_arguments("relaunch_simulations_scheduler", definition, description="Relaunch simulations as jobs on a scheduling system")

# -----------------------------------------------------------------

# Set simulation names
if config.from_directories:
    if config.names is not None: raise ValueError("Cannot specify names with 'from_directories' enabled")
    config.names = fs.directories_in_path(returns="name")

# -----------------------------------------------------------------

simulations = []

# -----------------------------------------------------------------

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

# Create SKIRT remote
remote = SKIRTRemote(host_id=config.host)

# -----------------------------------------------------------------

# Get the simulation status
simulation_status = remote.get_status(as_dict=True, returns=("name", "status"), simulations=simulations)

# -----------------------------------------------------------------

# Create logging options
logging = LoggingOptions()
logging.set_options(config.logging)

# -----------------------------------------------------------------

def add_to_queue(simulation):

    """
    This function ...
    :param simulation:
    :return:
    """

    # Set scheduling options
    scheduling_options = SchedulingOptions()
    scheduling_options.nodes = 1
    scheduling_options.ppn = remote.host.cluster.sockets_per_node * remote.host.cluster.cores_per_socket  # full node
    scheduling_options.full_node = True
    scheduling_options.walltime = config.walltime
    scheduling_options.mail = config.mail
    if config.jobscripts_path is not None: scheduling_options.local_jobscript_path = fs.join(config.jobscripts_path, simulation_name + ".sh")

    # Set remote input path
    remote_input_path = simulation.remote_input_path

    # Add existing simulation to the queue
    remote.add_to_queue(simulation.definition, name=simulation_name, logging_options=logging,
                        parallelization=config.parallelization, simulation=simulation, simulation_id=simulation_id,
                        clear_existing=True, dry=config.dry, scheduling_options=scheduling_options,
                        remote_input_path=remote_input_path)

# -----------------------------------------------------------------

# Loop over the simulations
for simulation in simulations:

    # Get simulation name
    simulation_id = simulation.id
    simulation_name = simulation.name

    # Get the status
    status = simulation_status[simulation_name]

    # Check status
    if is_finished_status(status):

        # Relaunch finished simulations
        if config.finished:

            # Remove local output, if present
            if simulation.retrieved:
                fs.clear_directory(simulation.output_path)
                simulation.set_retrieved(False)
                simulation.save()
            if simulation.analysed_any_extraction:
                fs.clear_directory(simulation.analysis.extraction.path)
                simulation.unset_analysed_extraction()
                simulation.save()
            if simulation.analysed_any_plotting:
                fs.clear_directory(simulation.analysis.plotting.path)
                simulation.unset_analysed_plotting()
                simulation.save()
            if simulation.analysed_any_misc:
                fs.clear_directory(simulation.analysis.misc.path)
                simulation.unset_analysed_misc()
                simulation.save()

            # Add to queue
            add_to_queue(simulation)

        # Don't relaunch finished simulations
        else: log.success("Simulation '" + simulation_name + "' is already finished")

    # Not finished: add to the queue again
    else:

        # Running?
        if status == "queued": raise NotImplementedError("Needs to be implemented")
        elif status == "running": raise NotImplementedError("Needs to be implemented")

        # Add the simulation to the queue
        add_to_queue(simulation)

# -----------------------------------------------------------------

# Launch the queue
handles = remote.start_queue(dry=config.dry)

# -----------------------------------------------------------------

# Set the handles, save the simulation objects
for simulation in simulations:

    # If all simulations should have the same handle
    if isinstance(handles, ExecutionHandle): simulation.handle = handles
    else: simulation.handle = handles[simulation.name]  # get the handle for this particular simulation

    # Save the simulation
    if not config.dry: simulation.save()

# -----------------------------------------------------------------
