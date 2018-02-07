#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.cache_generation_images Cache datacubes and mock observed images created for the simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.core.remote.remote import Remote
from pts.modeling.core.steps import cached_directory_path_for_single_command
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd

# -----------------------------------------------------------------

# Load the modeling environment
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# Fitting run
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation
definition.add_required("generation", "string", "name of the generation")

# The remote host
definition.add_optional("remote", "string", "remote host on which to cache the images")

# Flags
#definition.add_flag("reset", "reset everything: retrieve images and remove remote directories", False)
definition.add_flag("datacubes", "also cache the simulated datacubes")

# Create configuration
config = parse_arguments("cache_generation_images", definition, "Cache datacubes and mock observed images created for the simulations of a generation")

# -----------------------------------------------------------------

# Load the generation
fitting_run = runs.load(config.run)
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Setup the remote
if config.remote is not None: host_id = config.remote
else: host_id = environment.cache_host_id
remote = Remote(host_id=host_id)

# -----------------------------------------------------------------

# Create remote fit path
remote_fit_path = fs.join(remote.home_directory, environment.galaxy_name + "_fit")
if not remote.is_directory(remote_fit_path): remote.create_directory(remote_fit_path)

# Create remote path for the generation
remote_run_path = fs.join(remote_fit_path, fitting_run.name)
if not remote.is_directory(remote_run_path): remote.create_directory(remote_run_path)
remote_generation_path = fs.join(remote_run_path, generation.name)
if not remote.is_directory(remote_generation_path): remote.create_directory(remote_generation_path)

# -----------------------------------------------------------------

# Reset?
if config.reset: raise NotImplementedError("Not yet implemented")

# -----------------------------------------------------------------

# Loop over the analysed simulations
simulations = generation.analysed_simulations
for simulation in simulations:

    # Inform the user
    log.info("Caching simulation '" + simulation.name + "' ...")

    # Get the simulation name
    simulation_name = simulation.name

    # Get the datacube path
    datacube_path = generation.get_simulation_datacube_path(simulation_name)

    # Get the images path
    images_path = generation.get_simulation_misc_image_fluxes_images_earth_path(simulation_name)

    # Make remote directory for the simulation
    remote_simulation_path = fs.join(remote_generation_path, simulation_name)
    if not remote.is_directory(remote_simulation_path): remote.create_directory(remote_simulation_path)

    # Datacubes?
    if config.datacubes:

        # Debugging
        log.debug("Caching simulation datacube ...")

        # Upload the datacube
        remote.upload_file_to(datacube_path, remote_simulation_path, remove=True, show_output=True, replace=True)

    # Debugging
    log.debug("Caching the mock observed images ...")

    # Upload the images
    remote_images_path = fs.join(remote_simulation_path, "images")
    if not remote.is_directory(remote_images_path): remote.create_directory(remote_images_path)

    # Upload
    remote.upload_files_in_path_to(images_path, remote_images_path, remove=True, show_output=True, extension="fits", replace=True)

# -----------------------------------------------------------------
