#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.backup_generation_cache Make a backup of the cached generation simulation data.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.remote.host import find_host_ids
from pts.core.remote.remote import Remote
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Load the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generation name
definition.add_required("generation", "string", "generation name")

# Remote host ID
definition.add_required("host_id", "string", "remote host ID to cache to", choices=all_host_ids)

# Create configuration
config = parse_arguments("backup_generation_cache", definition, "Make a backup of the cached generation simulation data")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Get the cache path
cache_path = generation.single_cache_path

# -----------------------------------------------------------------

# Create remote
remote = Remote(host_id=config.host_id)

# -----------------------------------------------------------------

# Determine RT Modeling backup directory path
remote_backup_path = fs.join(remote.home_directory, "RT Modeling Backup")
remote_backup_galaxy_path = fs.join(remote_backup_path, environment.galaxy_name)

# -----------------------------------------------------------------

# Create path
if not remote.is_directory(remote_backup_galaxy_path): remote.create_directory(remote_backup_galaxy_path, recursive=True)

# -----------------------------------------------------------------

# Synchronize
remote.synchronize(cache_path, remote_backup_galaxy_path, show_output=True)

# -----------------------------------------------------------------
