#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.clear_generation_simulations Remove the simulation files for a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_proceed
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.core.tools.stringify import tostr
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.tools import sequences
from pts.core.simulation.remote import get_simulation_paths_for_host

# -----------------------------------------------------------------

# Get the modeling commands
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_positional_optional("generations", "string_list", "generation(s) for which to remove the simulations")
definition.add_optional("remote", "string", "remote host for which to remove the simulations (none means all")

# Remove simulations of generations other than certain generations
definition.add_optional("other_generations", "string_list", "remove all other generations than these generations")
definition.add_flag("hard", "remove ALL simulations objects that are not from the 'other_generations' (use with CARE!)")

# Make backup?
definition.add_flag("backup", "make backup of simulation files")
definition.add_optional("backup_path", "directory_path", "backup directory path")

# Get configuration
config = parse_arguments("clear_generation_simulations", definition, "Remove the simulation files for a certain generation")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# HARD delete
if config.hard:

    # Initialize a dictionary with the simulation filepaths to keep per host
    keep_paths = defaultdict(list)

    # Loop over the generations of which we CANNOT remove the simulations
    for generation_name in config.other_generations:

        # Check
        if not fitting_run.is_generation(generation_name): raise ValueError("Generation '" + generation_name + "' doesn't exist")

        # Get generation path
        generation_path = fitting_run.get_generation_path(generation_name)

        # Get the generation
        generation = fitting_run.get_generation(generation_name)

        # Loop over the remote hosts
        for host_id in generation.host_ids:

            # Get paths of the simulation files
            filepaths = generation.get_simulation_paths_for_host(host_id, id_or_name="name", not_exist="ignore", as_dict=False)
            nfilepaths = len(filepaths)

            # Add the filepaths
            keep_paths[host_id].extend(filepaths)

    # Loop over the remote hosts
    if config.remote is not None: host_ids = [config.remote]
    else: host_ids = keep_paths.keys()

    # Loop over the hosts
    for host_id in host_ids:

        # Get simulation paths
        paths = get_simulation_paths_for_host(host_id)

        #keep = keep_paths[host_id]
        #print(keep)

        # Get list of filepaths to remove
        remove_paths = sequences.get_other(paths, keep_paths[host_id])
        remove_names = [fs.strip_extension(fs.name(filepath)) for filepath in remove_paths]

        # Proceed?
        if not prompt_proceed("proceed removing the simulation files: " + tostr(remove_names)): continue

        # Backup
        if config.backup:
            log.debug("Creating backup of the simulation files ...")
            if config.remote is not None: backup_path = config.backup_path
            else: backup_path = fs.create_directory_in(config.backup_path, host_id)
            fs.copy_files(remove_paths, backup_path)

        # Remove
        fs.remove_files(remove_paths)

# Look only at other generations of the same fitting run
else:

    # Determine generation names
    if config.generations is not None:

        if config.other_generations is not None: raise ValueError("Cannot specify also 'other_generations'")
        generation_names = config.generations

    # Other generations
    elif config.other_generations is not None:

        # Get the other names than those specified
        generation_names = sequences.get_other(fitting_run.generation_names, config.other_generations)

    # Which generations?
    else: raise ValueError("Generation names cannot be determined")

    # Show
    log.debug("Removing the simulation files for generations: '" + tostr(generation_names) + "' ...")

    # Loop over the generations
    for generation_name in generation_names:

        # Check
        if not fitting_run.is_generation(generation_name): raise ValueError("Generation '" + generation_name + "' doesn't exist")

        # Get generation path
        generation_path = fitting_run.get_generation_path(generation_name)

        # Get the generation
        generation = fitting_run.get_generation(generation_name)

        # Loop over the remote hosts
        for host_id in generation.host_ids:

            # Clear for this host?
            if config.remote is not None and host_id != config.remote: continue

            # Get paths of the simulation files
            filepaths = generation.get_simulation_paths_for_host(host_id, id_or_name="name", not_exist="ignore")
            nfilepaths = len(filepaths)
            if nfilepaths == 0:
                log.warning("No simulation objects anymore for generation '" + generation_name + "'")
                continue # skip generation

            # Proceed?
            remove_names = [fs.strip_extension(fs.name(filepath)) for filepath in filepaths.keys()]
            if not prompt_proceed("proceed removing the simulation files of generation '" + generation_name + "': " + tostr(remove_names)): continue

            # Get the filepaths
            paths = filepaths.values()

            # Backup
            if config.backup:
                log.debug("Creating backup of the simulation files of generation '" + generation_name + "' ...")
                fs.copy_files(paths, config.backup_path)

            # Remove
            fs.remove_files(paths)

# -----------------------------------------------------------------
