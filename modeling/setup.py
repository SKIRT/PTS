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

def setup(command_name, cwd, configuration_method_argument=None):

    """
    This function ...
    :param command_name:
    :param cwd:
    :param configuration_method_argument:
    :return: 
    """

    # Add command to history
    mark_start(command_name, cwd)

    # Add command to commands file
    add_command(command_name, cwd, configuration_method_argument=configuration_method_argument)

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
    if configuration_method_argument is not None: argument_string += configuration_method_argument
    argument_string += command_name + " " + " ".join(arguments)

    # Load the commmands
    from .component.component import load_modeling_commands
    commands = load_modeling_commands(cwd)

    # Add entry
    commands.append(argument_string)

    # Save
    commands.save()

# -----------------------------------------------------------------

def finish(command_name, cwd):

    """
    This function ...
    :param command_name: 
    :param cwd:
    :return: 
    """

    # Mark end in history
    mark_end(command_name, cwd)

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
