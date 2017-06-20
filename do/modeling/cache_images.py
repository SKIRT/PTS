#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.cache_images Cache original images for which initialization has been performed.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.core.tools import filesystem as fs
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.core.remote.remote import Remote
from pts.modeling.core.steps import cached_directory_name_for_single_command
from pts.modeling.component.galaxy import get_data_image_and_error_paths, get_initial_dataset

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Setup the remote
remote = Remote(host_id=environment.cache_host_id)

# -----------------------------------------------------------------

directory_name = cached_directory_name_for_single_command(environment, "initialize_preparation")
remote_data_path = fs.join(remote.home_directory, directory_name)
if not remote.is_directory(remote_data_path): remote.create_directory(remote_data_path)

# -----------------------------------------------------------------

# Get the image paths
paths, error_paths = get_data_image_and_error_paths(modeling_path)

# -----------------------------------------------------------------

# Get the initial dataset
dataset = get_initial_dataset(modeling_path)

# -----------------------------------------------------------------

# Loop over the images
for name in paths:

    # Check whether the name is present in the dataset
    if name not in dataset.names:
        log.warning("The '" + name + "' image is not yet in the initial dataset, skipping ...")
        continue
    else:
        prep_path = fs.join(environment.prep_path, name)
        if not fs.is_directory(prep_path):
            log.warning("The '" + name + "' preparation directory is not yet present, skipping ...")
            continue
        elif fs.is_empty(prep_path):
            log.warning("The '" + name + "' preparation directory is empty, skipping ...")
            continue
        initialized_path = fs.join(prep_path, "initialized.fits")
        if not fs.is_file(initialized_path):
            log.warning("The initialized '" + name + "' image is not found, skipping ...")
            continue

    # Cache

    # Create directory for the image
    # Get the path
    path = paths[name]

    # Get filename
    filename = fs.strip_extension(fs.name(path))

    # Get directory name
    name = fs.name(fs.directory_of(path))

    # Create remote directory
    remote_image_directory = fs.join(remote_data_path, name)
    if not remote.is_directory(remote_image_directory): remote.create_directory(remote_image_directory)

    # Debugging
    log.info("Caching " + name + " image ...")

    # Debugging
    log.info("Uploading the file to '" + remote_image_directory + "' ...")

    # Upload
    remote_path = remote.upload_file_to(path, remote_image_directory, remove=True)

    # Succes
    log.success("File uploaded to '" + remote_path + "'")

    # Cache error path
    if name in error_paths:

        error_path = error_paths[name]

        # Info
        log.info("Caching the " + name + " error map ...")

        # Upload
        remote_path = remote.upload_file_to(error_path, remote_image_directory, remove=True)

        # Success
        log.success("File uploaded to '" + remote_path + "'")

# -----------------------------------------------------------------