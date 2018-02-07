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

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.screen import ScreenScript
from pts.core.simulation.remote import SKIRTRemote
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.remote.host import find_hosts, find_host_ids
from pts.core.simulation.remote import get_simulation_for_host, get_simulations_for_host
from pts.core.tools import sequences
from pts.core.tools.stringify import tostr
from pts.core.launch.options import LoggingOptions
from pts.core.advanced.runtimeestimator import RuntimeEstimator
from pts.core.launch.timing import TimingTable

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
definition.add_optional("runtime", "time_quantity", "estimated runtime for the simulations")

# Logging options
definition.import_section_from_composite_class("logging", "logging options", LoggingOptions)

# Read the command line arguments
config = parse_arguments("reschedule_simulations", definition, description="Relaunch simulations in a screen from a local script file")

# -----------------------------------------------------------------

#fs.set_nallowed_open_files(2000)

# -----------------------------------------------------------------

def is_finished_status(stat):

    """
    This function ...
    :param stat:
    :return:
    """

    return stat in ["finished", "retrieved", "analysed"]

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
        log.warning("The simulations are still in the queue")

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

    # Get the names
    #simulation_names =

    # Get the simulations
    simulations = dict()
    for simulation_id in config.ids:

        # Load the simulation
        simulation = get_simulation_for_host(config.from_host, simulation_id)

        # Add the simulation under its name
        simulations[simulation.name] = simulation

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
        log.debug("The estimated runtime for this host is " + str(runtime) + " seconds")

    # No timing info
    else: raise ValueError("When scheduling to remote with scheduling system, runtime or timing table path must be specified")

# -----------------------------------------------------------------

# Loop over the simulations
for simulation_name in simulations:

    print(simulation_name)

    simulation = simulations[simulation_name]
    handle = simulation.handle
    #print(handle)

    # Get definition
    #definition = simulation.definition
    #print(definition)
    logging = logging_options[simulation_name] if logging_options is not None else None
    #arguments = simulation.get_arguments(logging_options=logging, parallelization=config.parallelization)
    #print(arguments)

    # Add to the queue
    new_simulation = remote.add_to_queue(simulation.definition, name=simulation_name, logging_options=logging,
                                         parallelization=config.parallelization, simulation_id=simulation.id, save_simulation=False)

# -----------------------------------------------------------------

# Launch the queue
remote.start_queue()

# -----------------------------------------------------------------
