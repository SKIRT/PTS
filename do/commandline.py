#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.commandline

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from operator import itemgetter

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

def show_all_available(scripts, tables=None):

    """
    This function ...
    :param scripts:
    :param tables:
    :return:
    """

    print("Available PTS do commands:")

    # The list that will contain the info about all the scripts / table commands
    info = []

    ## Combine scripts and tables

    # Scripts
    for script in scripts:

        subproject = script[0]
        command = script[1]

        # Determine the path to the script
        path = fs.join(introspection.pts_package_dir, "do", subproject, command)

        # Get the description
        description = get_description(path, subproject, command[:-3])

        # Add entry to the info
        info.append((subproject, command[:-3], description, "arguments"))

    # Tables
    for subproject in tables:

        table = tables[subproject]
        for i in range(len(table["Command"])):

            command = table["Command"][i]
            hidden = False
            if command.startswith("*"):
                hidden = True
                command = command[1:]
            if hidden: continue # skip hidden
            description = table["Description"][i]
            configuration_method = table["Configuration method"][i]
            info.append((subproject, command, description, configuration_method))

    # Sort on the 'do' subfolder name
    info = sorted(info, key=itemgetter(0))

    lengths = dict()
    lengths_config_methods = dict()
    for script in info:
        subproject = script[0]
        command = script[1]
        if subproject not in lengths or len(command) > lengths[subproject]: lengths[subproject] = len(command)

        config_method = script[3] if len(script) == 4 else ""
        if subproject not in lengths_config_methods or len(config_method) > lengths_config_methods[subproject]:
            lengths_config_methods[subproject] = len(config_method)

    current_dir = None
    for script in info:

        from_table = len(script) == 4

        configuration_method = script[3] if from_table else ""
        description = script[2] if from_table else ""

        nspaces = lengths[script[0]] - len(script[1]) + 3
        pre_config_infix = " " * nspaces + "[" if from_table else ""
        after_config_infix = "]   " + " " *(lengths_config_methods[script[0]] - len(configuration_method)) if from_table else ""

        coloured_subproject = fmt.green + script[0] + fmt.reset
        coloured_name = fmt.bold + fmt.underlined + fmt.red + script[1] + fmt.reset
        coloured_config_method = fmt.yellow + pre_config_infix + configuration_method + after_config_infix + fmt.reset

        if current_dir == script[0]:
            print(" " * len(current_dir) + "/" + coloured_name + coloured_config_method + description)
        else:
            print(coloured_subproject + "/" + coloured_name + coloured_config_method + description)
            current_dir = script[0]

# -----------------------------------------------------------------

def show_possible_matches(matches, table_matches=None, tables=None):

    """
    This function ...
    :param matches:
    :param table_matches:
    :param tables:
    :return:
    """

    print("The command you provided is ambiguous. Possible matches:")

    # The list that will contain the info about all the scripts / table commands
    info = []

    ## Combine script and table matches

    # Scripts
    for script in matches:

        subproject = script[0]
        command = script[1]

        # Determine the path to the script
        path = fs.join(introspection.pts_package_dir, "do", subproject, command)

        # Get the description
        description = get_description(path, subproject, command[:-3])

        # Add entry to the info
        info.append((subproject, command[:-3], description, "arguments"))

    # Tables
    for subproject, index in table_matches:

        command = tables[subproject]["Command"][index]
        hidden = False
        if command.startswith("*"):
            hidden = True
            command = command[1:]
        if hidden: continue # skip if hidden
        description = tables[subproject]["Description"][index]
        configuration_method = tables[subproject]["Configuration method"][index]
        info.append((subproject, command, description, configuration_method))

    # Sort on the 'do' subfolder name
    info = sorted(info, key=itemgetter(0))

    lengths = dict()
    lengths_config_methods = dict()
    for script in info:
        subproject = script[0]
        command = script[1]
        if subproject not in lengths or len(command) > lengths[subproject]: lengths[subproject] = len(command)

        config_method = script[3] if len(script) == 4 else ""
        if subproject not in lengths_config_methods or len(config_method) > lengths_config_methods[subproject]:
            lengths_config_methods[subproject] = len(config_method)

    current_dir = None
    for script in info:

        from_table = len(script) == 4

        configuration_method = script[3] if from_table else ""
        description = script[2] if from_table else ""

        nspaces = lengths[script[0]] - len(script[1]) + 3
        pre_config_infix = " " * nspaces + "[" if from_table else ""
        after_config_infix = "]   " + " " * (lengths_config_methods[script[0]] - len(configuration_method)) if from_table else ""

        coloured_subproject = fmt.green + script[0] + fmt.reset
        coloured_name = fmt.bold + fmt.underlined + fmt.red + script[1] + fmt.reset
        coloured_config_method = fmt.yellow + pre_config_infix + configuration_method + after_config_infix + fmt.reset

        if current_dir == script[0]:
            print(" " * len(current_dir) + "/" + coloured_name + coloured_config_method + description)
        else:
            print(coloured_subproject + "/" + coloured_name + coloured_config_method + description)
            current_dir = script[0]

# -----------------------------------------------------------------

def get_description(script_path, subproject, command):

    """
    This function ...
    :param script_path:
    :param subproject:
    :param command:
    :return:
    """

    description = ""

    with open(script_path, 'r') as script:

        for line in script:

            line = line.rstrip("\n")

            if "## \package" in line:

                splitted = line.split(subproject + "." + command + " ")

                if len(splitted) > 1: description = splitted[1]
                break

    # Return description
    return description

# -----------------------------------------------------------------
