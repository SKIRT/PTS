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
import time as _time

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools import formatting as fmt
from pts.core.tools import time

# -----------------------------------------------------------------

pts_settings_names = ["version", "interactive", "arguments", "configfile", "rerun", "remote", "keep", "input", "output", "input_files", "output_files"]

expects_argument = dict()
expects_argument["version"] = False
expects_argument["interactive"] = False
expects_argument["arguments"] = False
expects_argument["configfile"] = True
expects_argument["rerun"] = False
expects_argument["remote"] = True
expects_argument["keep"] = False
expects_argument["input"] = True
expects_argument["output"] = True
expects_argument["input_files"] = True
expects_argument["output_files"] = True

# -----------------------------------------------------------------

class Command(object):

    """
    This class ...
    """

    def __init__(self, command, description, settings, input_dict, cwd=None, finish=None, pts_settings=None):

        """
        This function ...
        :param command:
        :param settings: dictionary
        :param input_dict: dictionary
        :param cwd:
        :param finish:
        :param pts_settings:
        """

        # Set working directory
        if cwd is None:
            self.cwd_specified = False
            cwd = fs.cwd()
        else: self.cwd_specified = True

        # Set attributes
        self.command = command
        self.description = description
        self.settings = settings
        self.input_dict = input_dict
        self.cwd = cwd
        self.finish = finish
        self.pts_settings = pts_settings

    # -----------------------------------------------------------------

    def __str__(self):

        """
        Thisf unction ...
        :return:
        """

        command = "pts "

        if self.pts_settings is not None:
            for name in self.pts_settings:
                command += "--" + name + " "
                if self.pts_settings[name] is not True: command += self.pts_settings[name] + " "

        from ..core.tools.stringify import represent_dict

        command += self.command + " "
        command += represent_dict(self.settings)

        # Return
        return command

# -----------------------------------------------------------------

def start_and_clear(command_name, target, **kwargs):

    """
    This function ...
    :param command_name: 
    :param target: 
    :param kwargs: 
    :return: 
    """

    # Start
    start_target(command_name, target, **kwargs)

    # Remove temp
    print("Clearing temporary data ...")
    introspection.remove_temp_dirs()

# -----------------------------------------------------------------

def start_target(command_name, target, **kwargs):

    """
    This function ...
    :return:
    """

    # Record starting time
    start = _time.time()

    # Run
    target(**kwargs)

    # Record end time
    end = _time.time()
    seconds = end - start

    # Succesfully finished
    print("Finished " + command_name + " in " + time.display_time(seconds))

# -----------------------------------------------------------------

def show_all_available(scripts, tables=None):

    """
    This function ...
    :param scripts:
    :param tables:
    :return:
    """

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
        title = None
        configuration_module_path = None
        info.append((subproject, command[:-3], description, "arguments", title, configuration_module_path))

    # Tables
    for subproject in tables:

        table = tables[subproject]

        # Loop over the entries
        for i in range(len(table["Command"])):

            command = table["Command"][i]
            hidden = False
            if command.startswith("*"):
                hidden = True
                command = command[1:]
            if hidden: continue # skip hidden
            description = table["Description"][i]
            configuration_method = table["Configuration method"][i]

            configuration_name = table["Configuration"][i]
            if configuration_name == "--": configuration_name = command
            configuration_module_path = "pts." + subproject + ".config." + configuration_name
            real_configuration_module_path = introspection.pts_root_dir + "/" + configuration_module_path.replace(".", "/") + ".py"
            if not fs.is_file(real_configuration_module_path): real_configuration_module_path = None

            title = table["Title"][i]
            info.append((subproject, command, description, configuration_method, title, real_configuration_module_path))

    # Sort on the 'do' subfolder name
    info = sorted(info, key=itemgetter(0))

    # Print in columns
    with fmt.print_in_columns(6) as print_row:

        previous_subproject = None
        previous_title = None

        for script in info:

            configuration_method = script[3]
            description = script[2]

            # Get the title
            title = script[4]

            # Get the configuration module path
            has_configuration = script[5] is not None
            if has_configuration: has_configuration_string = fmt.green + "C" + fmt.reset
            else: has_configuration_string = fmt.red + "NC" + fmt.reset

            # Transform title
            if title is None: title = ""

            coloured_subproject = fmt.green + script[0] + fmt.reset
            coloured_name = fmt.bold + fmt.underlined + fmt.red + script[1] + fmt.reset
            coloured_config_method = fmt.yellow + "[" + configuration_method + "]" + fmt.reset

            coloured_title = fmt.blue + fmt.bold + title + fmt.reset

            if previous_subproject is None or previous_subproject != script[0]:

                previous_subproject = script[0]
                print_row(coloured_subproject, "/")

            if previous_title is None or previous_title != title:

                previous_title = title

                if title != "":
                    print_row()
                    print_row("", "", coloured_title)

            # Print command
            print_row("", "/", coloured_name, coloured_config_method, description, has_configuration_string)

# -----------------------------------------------------------------

def show_possible_matches(matches, table_matches=None, tables=None):

    """
    This function ...
    :param matches:
    :param table_matches:
    :param tables:
    :return:
    """

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
        title = None
        configuration_module_path = None
        info.append((subproject, command[:-3], description, "arguments", title, configuration_module_path))

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

        configuration_name = tables[subproject]["Configuration"][index]
        if configuration_name == "--": configuration_name = command
        configuration_module_path = "pts." + subproject + ".config." + configuration_name
        real_configuration_module_path = introspection.pts_root_dir + "/" + configuration_module_path.replace(".", "/") + ".py"
        if not fs.is_file(real_configuration_module_path): real_configuration_module_path = None

        title = tables[subproject]["Title"][index]
        info.append((subproject, command, description, configuration_method, title, real_configuration_module_path))

    # Sort on the 'do' subfolder name
    info = sorted(info, key=itemgetter(0))

    # Print in columns
    with fmt.print_in_columns(6) as print_row:

        previous_subproject = None
        previous_title = None

        for script in info:

            configuration_method = script[3]
            description = script[2]

            # Get the title
            title = script[4]

            # Transform title
            if title is None: title = ""

            # Get the configuration module path
            has_configuration = script[5] is not None
            if has_configuration: has_configuration_string = fmt.green + "C" + fmt.reset
            else: has_configuration_string = fmt.red + "NC" + fmt.reset

            coloured_subproject = fmt.green + script[0] + fmt.reset
            coloured_name = fmt.bold + fmt.underlined + fmt.red + script[1] + fmt.reset
            coloured_config_method = fmt.yellow + "[" + configuration_method + "]" + fmt.reset

            coloured_title = fmt.blue + fmt.bold + title + fmt.reset

            if previous_subproject is None or previous_subproject != script[0]:
                previous_subproject = script[0]
                print_row(coloured_subproject, "/")

            if previous_title is None or previous_title != title:

                previous_title = title

                if title != "":
                    print_row()
                    print_row("", "", coloured_title)

            # Print command
            print_row("", "/", coloured_name, coloured_config_method, description, has_configuration_string)

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

    description = description.replace("\c", "")

    # Return description
    return description

# -----------------------------------------------------------------

def print_welcome():

    """
    This function ...
    :return: 
    """

    print("")
    print("  ### Welcome to PTS ### ")
    print("")

# -----------------------------------------------------------------
