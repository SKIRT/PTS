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

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

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

# # Determine number of retrieved and analysed simulations
# if config.lazy:
#     nretrieved = generation.nretrieved_simulations_basic
#     nanalysed = generation.nanalysed_simulations_basic
# else:
#     nretrieved = generation.nretrieved_simulations
#     nanalysed = generation.nanalysed_simulations
#
# # -----------------------------------------------------------------
#
# # Determine fraction analysed
# fraction_retrieved = float(nretrieved) / nsimulations
# fraction_analysed = float(nanalysed) / nsimulations
#
# # Show
# print("")
# print(fmt.bold + "Total number of simulations: " + fmt.reset + str(nsimulations))
# print(fmt.bold + "Number of retrieved simulations: " + fmt.reset + str(nretrieved) + " (" + tostr(fraction_retrieved*100, round=True, ndigits=2) + "%)")
# print(fmt.bold + "Number of analysed simulations: " + fmt.reset + str(nanalysed) + " (" + tostr(fraction_analysed*100, round=True, ndigits=2) + "%)")
# print("")

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
states = dict()
if not config.offline:

    # Inform the user
    log.info("Loading the remotes ...")

    # Create the remote instances
    for host_id in generation.host_ids: remotes[host_id] = SKIRTRemote(host_id=host_id)

    # Get screen states
    for host_id in remotes: states[host_id] = remotes[host_id].screen_states()

# -----------------------------------------------------------------

# Get the parameters table
parameters = generation.parameters_table

# Get parameter units
parameter_units = [parameters.unit_for(label) for label in fitting_run.free_parameter_labels]

# -----------------------------------------------------------------

# Number of digits for parameter values
parameters_ndigits = 3

# -----------------------------------------------------------------

# # Print header
# if config.parameters:
#
#     if config.extra: nspaces = 52
#     else: nspaces = 40
#
#     # Print
#     print(" " * nspaces + "\t" + "\t".join(fitting_run.free_parameter_labels))
#     print(" " * nspaces + "\t" + "\t".join(["[" + tostr(unit) + "]   " for unit in parameter_units]))

# -----------------------------------------------------------------

# # Loop over the simulations
# nfinished = 0
# for simulation_name in generation.simulation_names:
#
#     # Get the simulation
#     if config.lazy or not generation.has_simulation(simulation_name):
#
#         # No simulation object
#         simulation = None
#
#         # Get host ID and simulation ID from generation assignment table
#         host_id = generation.get_host_id(simulation_name)
#         simulation_id = generation.get_simulation_id(simulation_name)
#
#     # Simulation file is not present anymore
#     else:
#
#         # Load the simulation
#         simulation = generation.get_simulation(simulation_name)
#
#         # Get properties
#         host_id = simulation.host_id
#         simulation_id = simulation.id
#
#     # Get the parameter values
#     if config.parameters: parameter_values = parameters.parameter_values_for_simulation(simulation_name)
#     else: parameter_values = None
#
#     # Get parameter string
#     if parameter_values is not None:
#         parameter_strings = [tostr(parameter_values[label].value, scientific=True, decimal_places=3) + "  " for label in fitting_run.free_parameter_labels]
#         parameters_string = "\t".join(parameter_strings)
#     else: parameters_string = ""
#
#     # Get extra info
#     if config.extra:
#         ndigits = 3
#         id_string = strings.integer(simulation_id, ndigits)
#         host_string = strings.to_length(host_id, 5)
#         extra_string = " (" + host_string + " " + id_string + ")"
#     else: extra_string = ""
#
#     # No simulation object
#     if simulation is None:
#
#         # Has chi squared
#         if generation.is_analysed(simulation_name):
#
#             nfinished += 1
#             chisq = chi_squared.chi_squared_for(simulation_name)
#             print(" - " + fmt.green + simulation_name + extra_string + ": " + strings.number(chisq, chisq_ndecimal, chisq_ndigits, fill=" ") + "\t" + parameters_string + fmt.reset)
#
#         # Has miscellaneous output, but no chi squared
#         elif generation.has_misc_output(simulation_name):
#             nfinished += 1
#             print(" - " + fmt.yellow + simulation_name + extra_string + ": no chisq" + "\t" + parameters_string + fmt.reset)
#
#         # Has plotting output
#         elif generation.has_plotting_output(simulation_name):
#             nfinished += 1
#             print(" - " + fmt.yellow + simulation_name + extra_string + ": plotted" + "\t" + parameters_string + fmt.reset)
#
#         # Has extraction output
#         elif generation.has_extraction_output(simulation_name):
#             nfinished += 1
#             print(" - " + fmt.yellow + simulation_name + extra_string + ": extracted" + "\t" + parameters_string + fmt.reset)
#
#         # Has simulation output
#         elif generation.is_retrieved(simulation_name):
#             nfinished += 1
#             print(" - " + fmt.yellow + simulation_name + extra_string + ": retrieved" + "\t" + parameters_string + fmt.reset)
#
#         # No simulation output
#         else: print(" - " + fmt.red + simulation_name + extra_string + ": unknown" + "\t"  + parameters_string + fmt.reset)
#
#     # Simulation object
#     else:
#
#         # Already analysed
#         if simulation.analysed:
#
#             # Finished
#             nfinished += 1
#
#             # Get chi squared
#             chisq = chi_squared.chi_squared_for(simulation_name)
#             print(" - " + fmt.green + simulation_name + extra_string +  ": " + strings.number(chisq, chisq_ndecimal, chisq_ndigits, fill=" ") + "\t" + parameters_string + fmt.reset)
#
#         elif simulation.retrieved:
#
#             nfinished += 1
#             print(" - " + fmt.yellow + simulation_name + ": not analysed" + fmt.reset)
#
#         else:
#
#             # Not yet retrieved, what is the status?
#             if host_id in remotes: simulation_status = remotes[host_id].get_simulation_status(simulation, screen_states=states[host_id])
#             else: simulation_status = " unknown"
#
#             # Show
#             if simulation_status == "finished":
#                 nfinished += 1
#                 print(" - " + fmt.yellow + simulation_name + extra_string + ": " + simulation_status + "\t" + parameters_string + fmt.reset)
#             else: print(" - " + fmt.red + simulation_name + extra_string + ": " + simulation_status + "\t" + parameters_string + fmt.reset)
#
# # -----------------------------------------------------------------
#
# # Show number of finished simulations
# print("")
# fraction_finished = float(nfinished) / nsimulations
# print(fmt.bold + "Number of finished simulations: " + fmt.reset + str(nfinished) + " (" + tostr(fraction_finished*100, round=True, ndigits=2) + "%)")

# -----------------------------------------------------------------

#simulation_names = []
status_list = []

# -----------------------------------------------------------------

# Loop over the simulations
for simulation_name in generation.simulation_names:

    # Get the simulation
    if config.lazy or not generation.has_simulation(simulation_name):

        # No simulation object
        simulation = None

        # Get host ID and simulation ID from generation assignment table
        host_id = generation.get_host_id(simulation_name)
        simulation_id = generation.get_simulation_id(simulation_name)

    # Simulation file is not present anymore
    else:

        # Load the simulation
        simulation = generation.get_simulation(simulation_name)

        # Get properties
        host_id = simulation.host_id
        simulation_id = simulation.id

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
            if host_id in remotes: simulation_status = remotes[host_id].get_simulation_status(simulation, screen_states=states[host_id])
            else: simulation_status = "unknown"

    # Add the status
    status_list.append(simulation_status)

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

# Set backup path
backup_path = fs.create_directory_in(fitting_run.generations_path, "backup_" + generation.name)
manager.config.backup_dir_path = backup_path

# Run the manager
manager.run(assignment=generation.assignment_table, timing=fitting_run.timing_table, memory=fitting_run.memory_table,
            status=status, info_tables=[parameters, chi_squared])

# -----------------------------------------------------------------

# Plot chi squared?
if config.plot_chisquared: chi_squared.distribution.plot(title="Chi squared values", xlogscale=True)

# -----------------------------------------------------------------
