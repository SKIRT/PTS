#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.reschedule_simulations Move simulations from the queue of one remote to the queue of another remote.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.screen import ScreenScript
from pts.core.simulation.remote import SKIRTRemote
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.remote.host import find_hosts, find_host_ids
from pts.core.simulation.remote import get_simulation_for_host, get_simulations_for_host
from pts.core.tools import sequences
from pts.core.launch.options import LoggingOptions
from pts.core.advanced.runtimeestimator import RuntimeEstimator
from pts.core.launch.timing import TimingTable
from pts.core.basics.handle import ExecutionHandle
from pts.core.launch.batchlauncher import SimulationAssignmentTable
from pts.core.launch.options import SchedulingOptions
from pts.core.simulation.remote import is_finished_status

# -----------------------------------------------------------------

all_host_ids = find_host_ids()
#all_hosts = find_hosts()

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Remote host
definition.add_required("host", "host", "remote host on which to schedule the simulations", choices=all_host_ids)

# Required arguments
definition.add_positional_optional("names", "string_list", "simulation names")
definition.add_optional("from_host", "host", "remote host from which to reschedule simulations", choices=all_host_ids)
definition.add_flag("from_directories", "use directory names as simulation names")
definition.add_optional("screen", "file_path", "script file path from where to get the simulation names at the end of the queue")
definition.add_optional("nsimulations", "positive_integer", "number of simulations to reschedule")

# Parallelization
definition.add_optional("parallelization", "parallelization", "parallelization scheme for the simulations")

# Assignment file
definition.add_optional("assignment", "file_path", "path of assignment file")

# Timing table
definition.add_optional("timing", "file_path", "path of timing table")
definition.add_optional("runtime", "duration", "estimated runtime for the simulations")

# Logging options
definition.import_section_from_composite_class("logging", "logging options", LoggingOptions)

# Flags
definition.add_flag("backup", "make backup of the screen script, assignment scheme, and the original simulation objects")
definition.add_optional("backup_path", "directory_path", "backup directory path")

# Jobscripts path
definition.add_optional("jobscripts_path", "directory_path", "path of the directory for the jobscripts to be saved locally")

# Flags
definition.add_flag("dry", "run in dry mode")
definition.add_flag("mail", "send mail when jobs start")

# Give remote input directory
definition.add_optional("remote_input_path", "string", "remote input directory")

# Read the command line arguments
config = parse_arguments("reschedule_simulations", definition, description="Relaunch simulations in a screen from a local script file")

# -----------------------------------------------------------------

logging_options = None

# -----------------------------------------------------------------

# From screen script
if config.screen is not None:

    # Load screen script
    screen = ScreenScript.from_file(config.screen)

    # Make selection
    # Specific number of simulations
    if config.nsimulations is not None:

        # Get specific number of simulations from the end of the queue
        simulation_names = sequences.get_last_values(screen.simulation_names, config.nsimulations)

        # Give warning that simulations are still in the queue (cannot stop the queue without stopping other simulations)
        log.warning("The simulations will still be in the queue of remote host '" + screen.host_id + "'")

    # Not-finished simulations from the screen
    else:

        # Create remote
        remote = SKIRTRemote(host_id=screen.host_id)

        # Get status
        status = remote.get_status(as_dict=True, returns=["name", "status"], simulation_names=screen.simulation_names)

        # Create a list of only the not-finished simulations
        simulation_names = [simulation_name for simulation_name in status if is_finished_status(status[simulation_name])]
        #print("finished simulations: " + tostr(simulation_names))

        # Check screen
        # if remote.is_active_screen(screen.name):
        #    log.info("Killing screen '" + screen.name + "' ...")
        #    remote.kill_screen(screen.name)

    # Get the simulations
    simulations = get_simulations_for_host(screen.host_id, names=simulation_names, as_dict=True)

    # Get logging options
    logging_options = dict()
    for simulation_name in simulations:
        logging = screen.get_logging_options(simulation_name)
        logging_options[simulation_name] = logging
        #print(logging)

# Simulation names are given
elif config.names is not None:

    # Check that the remote host has been specified
    if config.from_host is None: raise ValueError("From host must be specified")

    # Set the simulation names
    simulation_names = config.names

    # Get the simulations
    simulations = get_simulations_for_host(config.from_host, names=simulation_names, as_dict=True)

# From directory names
elif config.from_directories:

    # Check that the remote host has been specified
    if config.from_host is None: raise ValueError("From host must be specified")

    # Get the names
    simulation_names = fs.directories_in_path(returns="name")

    # Get the simulations
    simulations = get_simulations_for_host(config.from_host, names=simulation_names, as_dict=True)

# Simulation IDs are specified
elif config.ids is not None:

    # Check that the remote host has been specified
    if config.from_host is None: raise ValueError("From host must be specified")

    # Get the simulations
    simulations = dict()
    for simulation_id in config.ids:

        # Load the simulation
        simulation = get_simulation_for_host(config.from_host, simulation_id)

        # Add the simulation under its name
        simulations[simulation.name] = simulation

    # Set the simulation names
    simulation_names = [simulation.name for simulation in simulations]

# Nothing
else: raise ValueError("Too little input")

# -----------------------------------------------------------------

# Create SKIRT remote
remote = SKIRTRemote(host_id=config.host)

# -----------------------------------------------------------------

# Remote works with scheduling system
if remote.scheduler:

    # Runtime is passed
    if config.runtime is not None: runtime = config.runtime

    # Create runtime estimator
    elif config.timing is not None:

        # Create runtime estimator
        timing = TimingTable.from_file(config.timing)
        runtime_estimator = RuntimeEstimator(timing)

        # Get

        # Estimate the runtime
        first_simulation_name = simulations.keys()[0]
        first_simulation = simulations[first_simulation_name]
        first_ski = first_simulation.ski_file

        # Get number of wavelengths and number of dust cells
        nwavelengths = first_ski.get_nwavelengths(input_path=first_simulation.input)
        ncells = first_ski.get_ncells(input_path=first_simulation.input)
        #print("nwavelengths:", nwavelengths)
        #print("ncells:", ncells)

        plot_path = None
        runtime = runtime_estimator.runtime_for(first_ski, config.parallelization, config.host.id, config.host.cluster_name, nwavelengths=nwavelengths, ncells=ncells, plot_path=plot_path)

        # Debugging
        log.debug("The estimated runtime for this host is " + str(runtime) + " seconds (" + str(runtime/3600) + " hours)")

    # No timing info
    else: raise ValueError("When scheduling to remote with scheduling system, runtime or timing table path must be specified")

# No scheduling system: runtime doesn't have to be calculated
else: runtime = None

# -----------------------------------------------------------------

shared_input = defaultdict(list)

# Check which simulations share input
for simulation_name in simulations:

    simulation = simulations[simulation_name]
    remote_input_path = simulation.remote_input_path
    if remote_input_path is None: continue

    # Add name for this input directory
    shared_input[remote_input_path].append(simulation_name)

# -----------------------------------------------------------------

shared_groups = shared_input.values()
#print(shared_groups)

# -----------------------------------------------------------------

shared_input_paths = dict()

# -----------------------------------------------------------------

# Create list of the new simulations
new_simulations = []

# Loop over the simulations
for simulation_name in simulations:

    # Inform the user
    log.info("Adding simulation '" + simulation_name + "' to the remote queue ...")

    # Get simulation
    simulation = simulations[simulation_name]
    #handle = simulation.handle

    # Get logging options
    logging = logging_options[simulation_name] if logging_options is not None else None

    # Is this simulation shared?
    if config.remote_input_path is not None:
        shared = False
        remote_input_path = config.remote_input_path
    else:
        shared = sequences.in_one(simulation_name, shared_groups, allow_more=False)
        if shared:
            # Determine key
            shared_key = tuple(sequences.pick_contains(shared_groups, simulation_name))
            if shared_key in shared_input_paths: remote_input_path = shared_input_paths[shared_key]
            else: remote_input_path = None
        else: shared_key = remote_input_path = None

    # Set scheduling options if necessary
    if remote.scheduler:
        scheduling_options = SchedulingOptions()
        scheduling_options.nodes = 1
        scheduling_options.ppn = remote.host.cluster.sockets_per_node * remote.host.cluster.cores_per_socket # full node
        scheduling_options.full_node = True
        scheduling_options.walltime = runtime
        if config.jobscripts_path is not None: scheduling_options.local_jobscript_path = fs.join(config.jobscripts_path, simulation_name + ".sh")
        scheduling_options.mail = config.mail
    else: scheduling_options = None

    # Show
    #print(logging)
    #print(config.parallelization)
    #print(scheduling_options)
    #continue

    # Add to the queue
    new_simulation = remote.add_to_queue(simulation.definition, name=simulation_name, logging_options=logging,
                                         parallelization=config.parallelization, remote_input_path=remote_input_path,
                                         clear_existing=True, dry=config.dry, scheduling_options=scheduling_options)

    #print(new_simulation.remote_input_path)

    ## SET PROPERTIES FROM ORIGINAL SIMULATION

    # ANALYSIS OPTIONS
    new_simulation.analysis = simulation.analysis
    new_simulation.analyser_paths = simulation.analyser_paths

    # Options for retrieval
    new_simulation.retrieve_types = simulation.retrieve_types

    # Options for removing remote or local input and output
    new_simulation.remove_remote_input = simulation.remove_remote_input
    new_simulation.remove_remote_output = simulation.remove_remote_output  # After retrieval
    new_simulation.remove_remote_simulation_directory = simulation.remove_remote_simulation_directory  # After retrieval
    new_simulation.remove_local_output = simulation.remove_local_output  # After analysis

    ##

    # If the input directory is shared between the different simulations
    if shared and shared_key not in shared_input_paths: shared_input_paths[shared_key] = new_simulation.remote_input_path

    # Add the new simulation
    new_simulations.append(new_simulation)

# -----------------------------------------------------------------

# Launch the queue
handles = remote.start_queue(dry=config.dry)

# -----------------------------------------------------------------

# Set the handles
for simulation in new_simulations:

    # If all simulations should have the same handle
    if isinstance(handles, ExecutionHandle): simulation.handle = handles
    else: simulation.handle = handles[simulation.name]  # get the handle for this particular simulation

    # Save the simulation
    if not config.dry: simulation.save()

# -----------------------------------------------------------------

# Make backups
if config.backup:

    # Create backup directory
    if config.backup_path is None: backup_path = fs.create_directory_in_cwd("backup_reschedule")
    else: backup_path = config.backup_path

    # Backup screen script
    if config.screen is not None: fs.copy_file(config.screen, backup_path)

    # Backup assignment table
    if config.assignment is not None: fs.copy_file(config.assignment, backup_path)

    # Backup the simulation files
    for simulation_name in simulations: fs.copy_file(simulations[simulation_name].path, backup_path)

# -----------------------------------------------------------------

# Remove the original simulation files
for simulation_name in simulations:

    simulation = simulations[simulation_name]
    if not config.dry: fs.remove_file(simulation.path)
    else: log.warning("[DRY] Not removing the simulation from '" + simulation.path + "' ...")

# -----------------------------------------------------------------

# Clear the list of simulations
simulations = []

# -----------------------------------------------------------------

# Adapt the screen script
if config.screen is not None:

    # Load screen script
    screen = ScreenScript.from_file(config.screen)

    # Remove the simulations that have been scheduled to the other remote
    for simulation_name in simulation_names:

        # Remove the simulation from the screen script
        #if not config.dry: screen.remove_simulation(simulation_name)
        #else: log.warning("[DRY] Not removing the simulation entry '" + simulation_name + "' from the screen script ...")
        screen.remove_simulation(simulation_name)

    # Save the screen script
    if not config.dry: screen.save()

# Adapt the assignment scheme
if config.assignment is not None:

    # Load the assignment scheme
    assignment = SimulationAssignmentTable.from_file(config.assignment)

    # Change the assignment
    for simulation in new_simulations:

        # Set properties
        assignment.set_host_for_simulation(simulation.name, config.host.id, config.host.cluster_name)
        assignment.set_id_for_simulation(simulation.name, simulation.id)
        assignment.set_success_for_simulation(simulation.name)

    # Save the assignment scheme
    if not config.dry: assignment.save()

# -----------------------------------------------------------------
