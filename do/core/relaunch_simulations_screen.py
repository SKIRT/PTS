#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.relaunch_simulations_screen Relaunch simulations in a screen from a local script file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.screen import ScreenScript
from pts.core.simulation.remote import SKIRTRemote
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.simulation.simulation import RemoteSimulation

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "script file path")

# Optional arguments
definition.add_optional("parallelization", "parallelization", "new parallelization scheme")
definition.add_optional("ncores", "positive_integer", "new number of cores")
definition.add_optional("nprocesses", "positive_integer", "new number of processes")
definition.add_optional("nthreads", "positive_integer", "new number of threads")
definition.add_optional("threads_per_core", "positive_integer", "new number of threads")
definition.add_flag("data_parallel", "new data parallelization flag", None)

# Flags
definition.add_flag("finished", "relaunch already finished simulations", True)
definition.add_flag("launch", "launch the screen again (instead of just adapting the screen script", True)
definition.add_flag("attached", "run in attached mode", False)
#definition.add_flag("direct", "use input screen script file directly", False)

# Read the command line arguments
config = parse_arguments("relaunch_simulations_screen", definition, description="Relaunch simulations in a screen from a local script file")

# -----------------------------------------------------------------

fs.set_nallowed_open_files(2000)

# -----------------------------------------------------------------

# Check arguments
if config.parallelization is not None:
    if config.ncores is not None or config.nprocesses is not None or config.nthreads is not None or config.threads_per_core is not None or config.data_parallel is not None:
        raise ValueError("Cannot specify parallelization and other parallelization specifications")

# -----------------------------------------------------------------

# Load screen script
screen = ScreenScript.from_file(config.filename)

# -----------------------------------------------------------------

# Direct?
#if config.direct: new_path = config.filename

# -----------------------------------------------------------------

# Adapt
#else:

# Loop over the arguments
for simulation_name in screen.simulation_names:

    # Set parallelization
    if config.parallelization is not None: screen.set_parallelization(simulation_name, config.parallelization)

    # Get arguments
    arguments = screen.get_arguments(simulation_name)
    parallelization = screen.get_parallelization(simulation_name)

    # Set parallelization settings
    if config.ncores is not None: parallelization.ncores = config.ncores
    if config.nprocesses is not None: parallelization.nprocesses = config.nprocesses
    if config.nthreads is not None: parallelization.nthreads = config.nthreads
    if config.threads_per_core is not None: parallelization.threads_per_core = config.threads_per_core
    if config.data_parallel is not None: parallelization.data_parallel = config.data_parallel
    arguments.parallelization = parallelization

# Save screen
new_path = fs.added_suffix(config.filename, "_relaunch")
screen.saveto(new_path)

# -----------------------------------------------------------------

if not config.launch: exit()

# -----------------------------------------------------------------

# Create SKIRT remote
remote = SKIRTRemote(host_id=screen.host_id)

# -----------------------------------------------------------------

# Check screen
if remote.is_active_screen(screen.name):
    log.info("Killing screen '" + screen.name + "' ...")
    remote.kill_screen(screen.name)

# -----------------------------------------------------------------

# Get simulations on the remote and their status
status = remote.get_status()
simulations = dict()
for filepath, stat in status:
    simulation = RemoteSimulation.from_file(filepath)
    simulation_name = simulation.name
    simulations[simulation_name] = (simulation, stat)

# -----------------------------------------------------------------

# Loop over the simulations of the screen script
for simulation_name in screen.simulation_names:

    # Check
    if simulation_name not in simulations: raise ValueError("Simulation '" + simulation_name + "' could not be found")

    # Get the simulation
    simulation = simulations[simulation_name][0]

    # Check status
    status = simulations[simulation_name][1]
    if status == "finished" or status == "retrieved" or status == "analysed":

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

        # Don't relaunch finished simulations
        else:

            # Remove from screen: not needed anymore
            screen.remove_simulation(simulation_name)
            screen.save()

    # Not finished
    else: pass

# -----------------------------------------------------------------

# Check whether any are left
if not screen.has_simulations:
    log.success("Not relaunching any simulation")
    exit()

# -----------------------------------------------------------------

# Start screen on the remote host
remote.start_screen(screen.name, new_path, remote.skirt_run_dir, screen_output_path=screen.output_path, keep_remote_script=True, attached=config.attached)

# -----------------------------------------------------------------
