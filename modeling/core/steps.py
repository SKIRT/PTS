#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.steps Modeling steps.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

all_modeling = "all"
galaxy_modeling = "galaxy"
sed_modeling = "sed"
images_modeling = "images"
modeling_types = [galaxy_modeling, sed_modeling, images_modeling]

# -----------------------------------------------------------------

single_commands = OrderedDict()
single_commands["fetch_properties"] = galaxy_modeling
single_commands["fetch_seds"] = galaxy_modeling
single_commands["fetch_images"] = galaxy_modeling
single_commands["inspect_data"] = galaxy_modeling
single_commands["initialize_preparation"] = galaxy_modeling
single_commands["inspect_initialization"] = galaxy_modeling
single_commands["prepare_data"] = galaxy_modeling
single_commands["inspect_preparation"] = galaxy_modeling
single_commands["decompose"] = galaxy_modeling
single_commands["truncate"] = galaxy_modeling
single_commands["set_significance_levels"] = galaxy_modeling
single_commands["photometry"] = galaxy_modeling
single_commands["make_colours_maps"] = galaxy_modeling
single_commands["make_ssfr_maps"] = galaxy_modeling
single_commands["make_tir_maps"] = galaxy_modeling
single_commands["make_attenuation_maps"] = galaxy_modeling
single_commands["make_dust_map"] = galaxy_modeling
single_commands["make_old_stellar_maps"] = galaxy_modeling
single_commands["make_young_stellar_maps"] = galaxy_modeling
single_commands["make_ionizing_stellar_maps"] = galaxy_modeling
single_commands["make_component_maps"] = galaxy_modeling
single_commands["plot_sed"] = sed_modeling # only for SEDModeler
single_commands["build_model_galaxy"] = galaxy_modeling
single_commands["build_model_sed"] = sed_modeling
single_commands["build_model_images"] = images_modeling
single_commands["generate_representations"] = galaxy_modeling
single_commands["build_representation_galaxy"] = galaxy_modeling
single_commands["build_representation_sed"] = sed_modeling
single_commands["build_representation_images"] = images_modeling
single_commands["configure_fit"] = all_modeling
single_commands["initialize_fit_sed"] = sed_modeling # for sed modeling
single_commands["initialize_fit_galaxy"] = galaxy_modeling # for galaxy modeling

# -----------------------------------------------------------------

def commands_for_modeling_type(modeling_type):

    """
    This function ...
    :param modeling_type:
    :return:
    """

    # Initialize a list for the
    commands = []
    for command in single_commands:

        type_for_command = single_commands[command]
        if type_for_command == modeling_type or type_for_command == all_modeling:
            commands.append(command)

    # Return the commands
    return commands

# -----------------------------------------------------------------

def galaxy_modeling_commands():

    """
    This function ...
    :return:
    """

    return commands_for_modeling_type(galaxy_modeling)

# -----------------------------------------------------------------

def sed_modeling_commands():

    """
    This function ...
    :return:
    """

    return commands_for_modeling_type(sed_modeling)

# -----------------------------------------------------------------

def images_modeling_commands():

    """
    This function ...
    :return:
    """

    return commands_for_modeling_type(images_modeling)

# -----------------------------------------------------------------

def output_paths_for_single_command(environment, command_name):

    """
    This function ...
    :param environment: the modeling environment object
    :param command_name: the command
    :return:
    """

    # Fetch properties
    if command_name == "fetch_properties": return [environment.galaxy_properties_path]

    # Fetch SEDs
    elif command_name == "fetch_seds": return [environment.data_seds_path]

    # Fetch images
    elif command_name == "fetch_images": return [environment.data_images_path]

    # Inspect data
    elif command_name == "inspect_data": return []

    # Initialize preparation
    elif command_name == "initialize_preparation": return []

    # Inspect initialization
    elif command_name == "inspect_initialization": return []

    # Prepare data
    elif command_name == "prepare_data": return []

    # Inspect preparation
    elif command_name == "inspect_preparation": return []

    # Decompose
    elif command_name == "decompose": return []

    # Trunate
    elif command_name == "truncate": return [environment.truncation_ellipse_path]

    # Significance levels
    elif command_name == "set_significance_levels": return [environment.significance_levels_path]

    # Photometry
    elif command_name == "photometry": return [environment.phot_path]

    # Make colour maps
    elif command_name == "make_colours_maps": return [environment.maps_colours_path]

    # sSFR MAPS
    elif command_name == "make_ssfr_maps": return [environment.maps_ssfr_path]

    # TIR maps
    elif command_name == "make_tir_maps": return [environment.maps_tir_path]

    # Attenuation maps
    elif command_name == "make_attenuation_maps": return [environment.maps_attenuation_path]

    # Dust map
    elif command_name == "make_dust_map": return [environment.maps_dust_path]

    # Old stars
    elif command_name == "make_old_stellar_maps": return [environment.maps_old_path]

    # Young stars
    elif command_name == "make_young_stellar_maps": return [environment.maps_young_path]

    # Ionizing stars
    elif command_name == "make_ionizing_stellar_maps": return [environment.maps_ionizing_path]

    # Component maps
    elif command_name == "make_component_maps": return [environment.maps_components_path]

    # Plot SED
    elif command_name == "plot_sed": return []

    # Build model
    elif command_name == "build_model_galaxy": return []

    # Build model for SED modeling
    elif command_name == "build_model_sed": return []

    # Build model for images modeling
    elif command_name == "build_model_images": return []

    #
    elif command_name == "generate_representations": return []

    # Build representation for galaxy model
    elif command_name == "build_representation_galaxy": return []

    # Build representation for SED modeling
    elif command_name == "build_representation_sed": return []

    # Build representation for images modeling
    elif command_name == "build_representation_images": return []

    # Configure fit
    elif command_name == "configure_fit": return []

    # Initialize fit for SED modeling
    elif command_name == "initialize_fit_sed": return []

    # Initialize fit for galaxy modeling
    elif command_name == "initialize_fit_galaxy": return []

    # Other
    else: raise ValueError("Invalid command: '" + command_name + "'")

# -----------------------------------------------------------------

def cached_directory_name_for_single_command(environment, command_name):

    """
    THis function ...
    :param environment:
    :param command_name:
    :return:
    """

    galaxy_name = environment.galaxy_name

    # Fetch properties
    if command_name == "fetch_properties": return None

    # Fetch SEDs
    elif command_name == "fetch_seds": return None

    # Fetch images
    elif command_name == "fetch_images": return None

    # Inspect data
    elif command_name == "inspect_data": return None

    # Initialize preparation
    elif command_name == "initialize_preparation": return galaxy_name + "_data"

    # Inspect initialization
    elif command_name == "inspect_initialization": return None

    # Prepare data
    elif command_name == "prepare_data": return galaxy_name + "_preparation"

    # Inspect preparation
    elif command_name == "inspect_preparation": return None

    # Decompose
    elif command_name == "decompose": return None

    # Trunate
    elif command_name == "truncate": return galaxy_name + "_truncation"

    # Significance levels
    elif command_name == "set_significance_levels": return None

    # Photometry
    elif command_name == "photometry": return None

    # Make colour maps
    elif command_name == "make_colours_maps": return None

    # sSFR MAPS
    elif command_name == "make_ssfr_maps": return None

    # TIR maps
    elif command_name == "make_tir_maps": return None

    # Attenuation maps
    elif command_name == "make_attenuation_maps": return None

    # Dust map
    elif command_name == "make_dust_map": return None

    # Old stars
    elif command_name == "make_old_stellar_maps": return None

    # Young stars
    elif command_name == "make_young_stellar_maps": return None

    # Ionizing stars
    elif command_name == "make_ionizing_stellar_maps": return None

    # RT Model component maps
    elif command_name == "make_component_maps": return None

    # Plot SED
    elif command_name == "plot_sed": return None

    # Build model
    elif command_name == "build_model_galaxy": return None

    # Build model for SED modeling
    elif command_name == "build_model_sed": return None

    # Build model for images modeling
    elif command_name == "build_model_images": return None

    #
    elif command_name == "generate_representations": return None

    # Build representation
    elif command_name == "build_representation_galaxy": return None

    # Build representation for SED modeling
    elif command_name == "build_representation_sed": return None

    # Build representation for images modeling
    elif command_name == "build_representation_images": return None

    # Configure fit
    elif command_name == "configure_fit": return None

    # Initialize fit for SED modeling
    elif command_name == "initialize_fit_sed": return None

    # Initialize fit for galaxy modeling
    elif command_name == "initialize_fit_galaxy": return None

    # Other
    else: raise ValueError("Invalid command: '" + command_name + "'")

# -----------------------------------------------------------------

def cached_directory_path_for_single_command(environment, command_name, remote):

    """
    This function ...
    :param environment:
    :param command_name:
    :param remote:
    :return:
    """

    name = cached_directory_name_for_single_command(environment, command_name)
    if name is None: return None
    else: return fs.join(remote.home_directory, name)

# -----------------------------------------------------------------

def single_commands_for_galaxy_modeling(command_name):

    """
    This function ...
    :param command_name:
    :return:
    """

    return single_commands[command_name] == galaxy_modeling or single_commands[command_name] == all_modeling

# -----------------------------------------------------------------

def single_commands_for_sed_modeling(command_name):

    """
    This function ...
    :param command_name:
    :return:
    """

    return single_commands[command_name] == sed_modeling or single_commands[command_name] == all_modeling

# -----------------------------------------------------------------

def single_commands_for_images_modeling(command_name):

    """
    This function ...
    :return:
    """

    return single_commands[command_name] == images_modeling or single_commands[command_name] == all_modeling

# -----------------------------------------------------------------

single_commands_galaxy = [command for command in single_commands if single_commands_for_galaxy_modeling(command)]

# -----------------------------------------------------------------

single_commands_sed = [command for command in single_commands if single_commands_for_sed_modeling(command)]

# -----------------------------------------------------------------

single_commands_images = [command for command in single_commands if single_commands_for_images_modeling(command)]

# -----------------------------------------------------------------

repeated_commands = ["explore", "fit_sed", "finish_exploration"]

# -----------------------------------------------------------------

def commands_before(command_name):

    """
    This function ...
    :param command_name:
    :return:
    """

    commands = []
    for command in single_commands:
        if command == command_name: break
        commands.append(command)
    return commands

# -----------------------------------------------------------------

def commands_after(command_name):

    """
    This function ...
    :param command_name:
    :return:
    """

    commands = []
    for command in reversed(single_commands):
        if command == command_name: break
        commands.append(command)
    return list(reversed(commands))

# -----------------------------------------------------------------

def commands_before_and_including(command_name):

    """
    This function ...
    :param command_name:
    :return:
    """

    return commands_before(command_name) + [command_name]

# -----------------------------------------------------------------

def commands_after_and_including(command_name):

    """
    This function ...
    :param command_name:
    :return:
    """

    return [command_name] + commands_after(command_name)
    
# -----------------------------------------------------------------
