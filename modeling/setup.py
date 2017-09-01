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

def setup(command_name, cwd):

    """
    This function ...
    :param command_name:
    :param cwd:
    :return: 
    """

    # Add command to history
    mark_start(command_name, cwd)

    # Add command to commands file
    add_command(command_name, cwd)

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

def add_command(command_name, cwd):

    """
    This fucntion ...
    :param command_name:
    :param cwd:
    :return:
    """

    # Get the argument string
    argument_string = "pts " + command_name + " " + " ".join(sys.argv[1:])

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
