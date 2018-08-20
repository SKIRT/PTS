#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.setup Contains a function that performs a setup.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from ..core.tools import introspection
from ..core.tools.strings import add_quotes_if_spaces
from ..core.basics.log import log
from ..core.tools import filesystem as fs

# -----------------------------------------------------------------

# Determine commands to be ignored for the history
ignore_titles = ["MODELING", "OTHER", "PLOTTING", "HIDDEN", "SHOW", "EXTRA", "HTML"]
ignore_commands = []
table = introspection.get_argument_table("modeling")
# Loop over all commands for which the title is one of the ignore titles
for index in range(len(table["Command"])):
    command_name = table["Command"][index]
    title = table["Title"][index]
    #print(list(title))
    if title in ignore_titles: ignore_commands.append(command_name)

# -----------------------------------------------------------------

# Determine commands for which the cwd doesn't need to be a modeling directory
not_modeling_cwd_commands = ["setup", "calculate_weights", "project_data"]

# -----------------------------------------------------------------

def setup(command_name, cwd, configuration_method_argument=None):

    """
    This function ...
    :param command_name:
    :param cwd:
    :param configuration_method_argument:
    :return: 
    """

    # Check directory
    # NOW CALLED IN DO/RUN.PY: TO AVOID CONFIG BEING WRITTEN OUT BEFORE
    #check_modeling_cwd(command_name, cwd)

    if command_name in not_modeling_cwd_commands: return

    # Add command to history
    mark_start(command_name, cwd)

    # Add command to commands file
    add_command(command_name, cwd, configuration_method_argument=configuration_method_argument)

# -----------------------------------------------------------------

def check_modeling_cwd(command_name, cwd):

    """
    This function ...
    :param command_name:
    :param cwd:
    :return:
    """

    if command_name in not_modeling_cwd_commands: return

    # Check whether this is a modeling directory
    from .core.environment import is_modeling_path
    if not is_modeling_path(cwd): raise ValueError("Not a modeling directory")

# -----------------------------------------------------------------

def mark_start(command_name, cwd):

    """
    This function ...
    :param command_name:
    :param cwd:
    :return:
    """

    # Ignore not-pipeline
    if command_name in ignore_commands: return

    # Load the history
    from .component.component import load_modeling_history
    history = load_modeling_history(cwd)

    # Add entry
    history.add_entry_and_save(command_name)

# -----------------------------------------------------------------

def add_command(command_name, cwd, configuration_method_argument=None):

    """
    This fucntion ...
    :param command_name:
    :param cwd:
    :param configuration_method_argument:
    :return:
    """

    # Get the command-line arguments
    arguments = [add_quotes_if_spaces(string) for string in sys.argv[1:]]

    # Get the argument string
    argument_string = "pts "
    if configuration_method_argument is not None and configuration_method_argument != "": argument_string += configuration_method_argument + " "
    argument_string += command_name + " " + " ".join(arguments)

    # Load the commmands
    from .component.component import load_modeling_commands
    commands = load_modeling_commands(cwd)

    # Add entry
    commands.append(argument_string)

    # Save
    commands.save()

# -----------------------------------------------------------------

def finish(command_name, cwd, config_path=None):

    """
    This function ...
    :param command_name: 
    :param cwd:
    :param config_path:
    :return: 
    """

    if command_name in not_modeling_cwd_commands: return

    # Mark end in history
    mark_end(command_name, cwd)

    # Tag the configuration file as succesful
    if config_path is not None: tag_config(config_path)

# -----------------------------------------------------------------

def mark_end(command_name, cwd):

    """
    This function ...
    :param command_name:
    :param cwd:
    :return:
    """

    # Ignore not-pipeline
    if command_name in ignore_commands: return

    # Load the history
    from .component.component import load_modeling_history
    history = load_modeling_history(cwd)

    # Mark end
    history.mark_end_and_save()

# -----------------------------------------------------------------

def tag_config(config_path, tag="*"):

    """
    This function ...
    :param config_path:
    :param tag:
    :return:
    """

    # Debugging
    log.debug("Tagging configuration file as succesful ...")

    # Rename by adding the tag
    fs.add_suffix(config_path, tag)

# -----------------------------------------------------------------
