#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_status View the status of the simulations of a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr
from pts.core.simulation.remote import SKIRTRemote
from pts.core.basics.log import log
from pts.core.tools import numbers
from pts.core.launch.manager import SimulationManager
from pts.core.launch.batchlauncher import SimulationStatusTable
from pts.core.tools import filesystem as fs
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulations_for_host

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Plotting
definition.add_flag("plot_runtimes", "plot runtimes")
definition.add_flag("plot_memory", "plot memory usage")
definition.add_flag("plot_chisquared", "plot chi squared")

# Show
definition.add_flag("show_runtimes", "show runtimes")
definition.add_flag("show_memory", "show memory")

# Show parameters
definition.add_flag("parameters", "show the parameter values")
definition.add_flag("extra", "show extra info")

# Flags
definition.add_flag("offline", "offline mode")
definition.add_flag("lazy", "lazy mode")
definition.add_flag("interactive", "interactive simulation management mode")
definition.add_flag("find_simulations", "find missing simulations by searching on simulation name")
definition.add_optional("find_remotes", "string_list", "find missing simulations in these remote hosts", default=all_host_ids, choices=all_host_ids)
definition.add_flag("dry", "run in dry mode")
definition.add_flag("write_status", "write the status", False)
definition.add_flag("write_commands", "write the commands", False)

# Get configuration
config = parse_arguments("generation_status", definition, "View the status of the simulations of a certain generation")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Get the chi squared table
chi_squared = generation.chi_squared_table

# -----------------------------------------------------------------

# Check
if not generation.has_assignment_table: raise RuntimeError("No assignment for this generation")

# -----------------------------------------------------------------

# Get number of simulations
nsimulations = generation.nsimulations

# -----------------------------------------------------------------

if len(chi_squared) > 0:

    # Get the maximum chi squared value
    min_chi_squared = chi_squared.best_chi_squared
    max_chi_squared = chi_squared.worst_chi_squared
    max_chi_squared_magnitude = numbers.order_of_magnitude(max_chi_squared)

    # Set number of digits for chi squared values
    chisq_ndecimal = 1
    chisq_ndigits = max_chi_squared_magnitude + 1 + chisq_ndecimal

    print(fmt.bold + "Best chi squared: " + fmt.reset + tostr(min_chi_squared))
    print(fmt.bold + "Worst chi squared: " + fmt.reset + tostr(max_chi_squared))
    print("")

# -----------------------------------------------------------------

remotes = dict()
screen_states = dict()
jobs_status = dict()

# -----------------------------------------------------------------

if not config.offline:

    # Inform the user
    log.info("Loading the remotes ...")

    # Create the remote instances
    for host_id in generation.host_ids: remotes[host_id] = SKIRTRemote(host_id=host_id)

    # Get screen states
    for host_id in remotes:

        if remotes[host_id].scheduler: jobs_status[host_id] = remotes[host_id].get_jobs_status()
        else: screen_states[host_id] = remotes[host_id].screen_states()

# -----------------------------------------------------------------

# Get the parameters table
parameters = generation.parameters_table

# Get parameter units
parameter_units = [parameters.unit_for(label) for label in fitting_run.free_parameter_labels]

# -----------------------------------------------------------------

# Number of digits for parameter values
parameters_ndigits = 3

# -----------------------------------------------------------------

# Set paths
manage_path = fs.create_directory_in(fitting_run.generations_path, "manage__" + generation.name)
current_indices = fs.directories_in_path(manage_path, returns="name", convert=int)
manage_current_path = fs.create_directory_in(manage_path, str(numbers.lowest_missing_integer(current_indices)))
backup_path = fs.join(manage_current_path, "backup")

# -----------------------------------------------------------------

status_list = []

# -----------------------------------------------------------------

# Initialize dictionary
simulations_hosts = dict()

# Get simulations for a remote host
def get_simulations(host_id):

    # Get the simulations
    simulations_host = get_simulations_for_host(host_id, as_dict=True)

    # Set the simulations
    simulations_hosts[host_id] = simulations_host

    # Return the simulations
    return simulations_host

# -----------------------------------------------------------------

def find_simulation(simulation_name, host_ids):

    """
    This function ...
    :param simulation_name:
    :param host_ids:
    :return:
    """

    simulation = None

    # Loop over the remotes
    for host_id in host_ids:

        # Check whether the simulation name is in the
        if simulation_name in get_simulations(host_id):
            #the_host_id = host_id
            simulation = simulations_hosts[host_id][simulation_name]
            assert simulation.host_id == host_id # check that the simulation object has the correct host ID as attribute
            #the_simulation_id = simulation.id
            break

    # Return the host ID and simulation ID
    #return the_host_id, the_simulation_id

    # Return the simulation
    return simulation

# -----------------------------------------------------------------

# Initialize flag for when assignment scheme is changed
changed_assignment = False

# -----------------------------------------------------------------

# Loop over the simulations
for simulation_name in generation.simulation_names:

    load_simulation = True
    if config.lazy: load_simulation = False
    if not generation.has_simulation(simulation_name):
        log.warning("The simulation file for '" + simulation_name + "' is not present anymore")
        host_id = generation.get_host_id(simulation_name)
        simulation_id = generation.get_simulation_id(simulation_name)
        log.warning("Simulation ID was '" + str(simulation_id) + "' and host ID is '" + host_id + "'")
        if config.find_simulations:
            log.warning("Looking for the simulation '" + simulation_name + "' on other remote hosts ...")
            #host_id, simulation_id = find_simulation(simulation_name, config.find_remotes)
            simulation = find_simulation(simulation_name, config.find_remotes)
            if simulation is None: raise RuntimeError("Simulation '" + simulation_name + "' was not found")
            actual_host_id = simulation.host_id
            actual_id = simulation.id
            cluster_name = simulation.cluster_name
            log.warning("Simulation identified as '" + str(actual_id) + "' for host '" + actual_host_id + "'")
            log.warning("Fixing assignment table ...")
            generation.set_id_and_host_for_simulation(simulation_name, actual_id, actual_host_id, cluster_name=cluster_name)
            changed_assignment = True
            load_simulation = True
        else: load_simulation = False

    # Load the simulation
    if load_simulation:

        # Load the simulation
        simulation = generation.get_simulation(simulation_name)

        # Get properties
        host_id = simulation.host_id
        simulation_id = simulation.id

    # Don't load the simulation
    else:

        # No simulation object
        simulation = None

        # Get host ID and simulation ID from generation assignment table
        host_id = generation.get_host_id(simulation_name)
        simulation_id = generation.get_simulation_id(simulation_name)

    # No simulation object
    if simulation is None:

        has_misc = generation.has_misc_output(simulation_name)
        has_plotting = generation.has_plotting_output(simulation_name)
        has_extraction = generation.has_extraction_output(simulation_name)

        # Has chi squared
        if generation.is_analysed(simulation_name): simulation_status = "analysed"

        # Has any analysis output
        elif has_extraction or has_plotting or has_misc:

            analysed = []
            if has_extraction: analysed.append("extraction")
            if has_plotting: analysed.append("plotting")
            if has_misc: analysed.append("misc")
            simulation_status = "analysed: " + ", ".join(analysed)

        # Has simulation output
        elif generation.is_retrieved(simulation_name): simulation_status = "retrieved"

        # No simulation output
        else: simulation_status = "unknown"

    # Simulation object
    else:

        # Already analysed
        if simulation.analysed: simulation_status = "analysed"

        # Partly analysed
        elif simulation.analysed_any:

            analysed = []
            if simulation.analysed_all_extraction: analysed.append("extraction")
            if simulation.analysed_all_plotting: analysed.append("plotting")
            if simulation.analysed_all_misc: analysed.append("misc")
            if simulation.analysed_batch: analysed.append("batch")
            if simulation.analysed_scaling: analysed.append("scaling")
            if simulation.analysed_all_extra: analysed.append("extra")

            if len(analysed) > 0: simulation_status = "analysed: " + ", ".join(analysed)
            else: simulation_status = "analysed: started"

        # Retrieved
        elif simulation.retrieved: simulation_status = "retrieved"

        # Not retrieved
        else:

            # Not yet retrieved, what is the status?
            screen_states_host = screen_states[host_id] if host_id in screen_states else None
            jobs_status_host = jobs_status[host_id] if host_id in jobs_status else None
            if host_id in remotes: simulation_status = remotes[host_id].get_simulation_status(simulation, screen_states=screen_states_host, jobs_status=jobs_status_host)
            else: simulation_status = "unknown"

    # Add the status
    status_list.append(simulation_status)

# -----------------------------------------------------------------

# Save the assignment table if it has been adapted
if changed_assignment:
    log.debug("Saving the changed assignment table for the generation ...")
    generation.assignment_table.save()

# -----------------------------------------------------------------

# Create the simulation status table
status = SimulationStatusTable.from_columns(generation.simulation_names, status_list)

# -----------------------------------------------------------------

# Create simulation manager
manager = SimulationManager()

# Set options
manager.config.show_status = not config.interactive
manager.config.show_runtimes = config.show_runtimes
manager.config.show_memory = config.show_memory
manager.config.plot_runtimes = config.plot_runtimes
manager.config.plot_memory = config.plot_memory
manager.config.interactive = config.interactive
if config.interactive: manager.config.commands = ["status"]
manager.config.dry = config.dry
manager.config.shared_input = True
manager.config.write_status = config.write_status
manager.config.write_moved = True
manager.config.write_relaunched = True
manager.config.write_commands = config.write_commands

# Set paths
manager.config.path = manage_current_path
manager.config.backup_dir_path = manage_current_path
manager.config.backup_dirname = "backup"

# Run the manager
manager.run(assignment=generation.assignment_table, timing=fitting_run.timing_table, memory=fitting_run.memory_table,
            status=status, info_tables=[parameters, chi_squared], remotes=remotes)

# -----------------------------------------------------------------

# Plot chi squared?
if config.plot_chisquared: chi_squared.distribution.plot(title="Chi squared values", xlogscale=True)

# -----------------------------------------------------------------

# Remove output directory if nothing was written
if fs.is_empty(manage_current_path, recursive=True): fs.remove_directory(manage_current_path)

# -----------------------------------------------------------------
