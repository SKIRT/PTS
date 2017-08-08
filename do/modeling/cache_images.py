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
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.core.remote.remote import Remote
from pts.modeling.core.steps import cached_directory_path_for_single_command
from pts.modeling.component.galaxy import get_data_image_and_error_paths, get_initial_dataset
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.data.component import instrument_to_origin
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_flag("reset", "reset everything: retrieve images and remove remote data folder", False)

# Create configuration
config = parse_arguments("cache_images", definition)

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Setup the remote
remote = Remote(host_id=environment.cache_host_id)

# -----------------------------------------------------------------

remote_data_path = cached_directory_path_for_single_command(environment, "initialize_preparation", remote)
if not remote.is_directory(remote_data_path): remote.create_directory(remote_data_path)

# -----------------------------------------------------------------

# Reset?
if config.reset:

    # Loop over all FITS files in the remote data directory
    for name, path in remote.files_in_path(remote_data_path, recursive=True, extension="fits", returns=["name", "path"]):

        # Determine origin
        origin = instrument_to_origin(name.split("_")[1])

        # Determine local directory for this image
        origin_path = fs.join(environment.data_images_path, origin)
        if not fs.is_directory(origin_path): fs.create_directory(origin_path)

        # Determine local path
        local_path = fs.join(origin_path, name)

        #print("local_path")

        # Check whether the image is not present
        if fs.is_file(local_path):
            log.warning("The '" + name + "' remotely cached image is still present locally. Keeping this file and throwing the remote file away.")
            continue
        else:
            # Infomr
            log.info("Downloading the '" + name + "' file to [" + origin_path + "] ...")

            # Download
            remote.download_file_to(path, origin_path, remove=True)

            # Success
            log.success("Succesfully retrieved the '" + name + "' image ...")

    # Inform the user
    log.info("Clearing the remote data structure ...")

    # Clear the entire data structure now
    remote.clear_directory(remote_data_path)

# -----------------------------------------------------------------

# Get the image paths
paths, error_paths = get_data_image_and_error_paths(modeling_path)

# -----------------------------------------------------------------

# Get the initial dataset
dataset = get_initial_dataset(modeling_path, check=False)

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
            # Check whether cached image is present
            remote_prep_path = cached_directory_path_for_single_command(environment, "prepare_data", remote)
            #print(remote_prep_path)
            remote_image_directory_path = fs.join(remote_prep_path, name)
            #print(remote_image_directory_path)
            remote_initialized_path = fs.join(remote_image_directory_path, "initialized.fits")
            if not remote.is_directory(remote_prep_path) or not remote.is_directory(remote_image_directory_path) or not remote.is_file(remote_initialized_path):
                log.warning("The initialized '" + name + "' image is not found, skipping ...")
                continue

    # Cache

    # Create directory for the image
    # Get the path
    path = paths[name]

    # Get filename
    filename = fs.strip_extension(fs.name(path))

    #print(filename)

    # Determine origin from filename
    #origin = instrument_to_origin(filename.split("_")[1])

    # Other way
    origin = fs.name(fs.directory_of(path))

    # CHeck
    #if origin != origin_alt: raise ValueError(origin + " != " + origin_alt)

    # Create remote directory
    remote_image_directory = fs.join(remote_data_path, origin)
    if not remote.is_directory(remote_image_directory): remote.create_directory(remote_image_directory)

    # Check whether the file exists
    if not fs.is_file(path):

        # Check whether already cached to the remote
        remote_image_path = fs.join(remote_image_directory, fs.name(path))

        if not remote.is_file(remote_image_path):
            raise IOError("The '" + name + "' image could not be found")
        else: log.success("The '" + name + "' image is already cached")

    # Debugging
    log.info("Caching " + name + " image ...")

    # Debugging
    log.info("Uploading the file to '" + remote_image_directory + "' ...")

    # Upload
    #remote_path = remote.upload_file_to(path, remote_image_directory, remove=True)
    remote.upload_file_to(path, remote_image_directory, remove=True)

    # Succes
    #log.success("File uploaded to '" + remote_path + "'")

    # Cache error path
    if name in error_paths:

        # Get the path to the error file
        error_path = error_paths[name]

        # Info
        log.info("Caching the " + name + " error map ...")

        # Upload
        #remote_path = remote.upload_file_to(error_path, remote_image_directory, remove=True)
        remote.upload_file_to(error_path, remote_image_directory, remove=True)

        # Success
        #log.success("File uploaded to '" + remote_path + "'")

    # Succes
    log.success("The '" + name + "' image is succesfully cached to the remote host '" + remote.host_id + "'")

# -----------------------------------------------------------------