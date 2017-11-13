#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.commands Contains the ModelingCommands class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools import strings
from ...do.commandline import Command, pts_settings_names, expects_argument
from ...core.tools import introspection
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

def get_modeling_commands():

    """
    Thisf unction ...
    :return:
    """

    return introspection.get_all_commands_in_subproject("modeling")

# -----------------------------------------------------------------

class ModelingCommands(list):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelingCommands, self).__init__(*args)

        # Set the path
        self.path = kwargs.pop("path", None)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        # Load the lines
        return cls(fs.get_lines(filepath), path=filepath)

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Write the lines
        fs.write_lines(path, self)

        # Update the path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    def as_lists(self):

        """
        Thisf unction ...
        :return:
        """

        lists = []
        for command in self: lists.append(strings.split_except_within_single_quotes(command))
        return lists

    # -----------------------------------------------------------------

    @lazyproperty
    def all_commands(self):

        """
        This function ...
        :return:
        """

        return get_modeling_commands()

    # -----------------------------------------------------------------

    def get_description(self, command_name):

        """
        Thisf unction ...
        :param command_name:
        :return:
        """

        return self.all_commands[command_name]

    # -----------------------------------------------------------------

    def get_definition(self, command_name, modeling_path):

        """
        This function ...
        :param command_name:
        :param modeling_path:
        :return:
        """

        from ...core.launch.pts import get_definition_for_command
        return get_definition_for_command(command_name, cwd=modeling_path)

    # -----------------------------------------------------------------

    def as_commands(self):

        """
        This function ...
        :return:
        """

        commands = []

        # Loop over the command lines
        for lst in self.as_lists():

            # PTS settings
            pts_settings = OrderedDict()

            # Find the actual command
            index = 1
            expect_argument = False
            for index, string in enumerate(lst[1:]):

                # PTS setting
                if string.startswith("-"):
                    setting_name = string.split("--")[1] if string.startswith("--") else string.split("-")[1]
                    #print(setting_name)
                    if setting_name not in pts_settings_names: raise ValueError("Invalid command")
                    if expects_argument[setting_name]: # expects arguments
                        expect_argument = True
                        pts_settings[setting_name] = None
                    else: pts_settings[setting_name] = True # is a flag

                # Argument from previous PTS setting
                elif expect_argument:

                    last_setting_name = pts_settings.keys()[-1]
                    pts_settings[last_setting_name] = strings.unquote(string)

                # Actual command
                else:
                    command_name = string
                    break

            # Get the settings
            settings_list = lst[index+2:]
            settings = OrderedDict()
            last_was_optional = False
            npositional = 0
            #print(settings_list)
            for index, string in enumerate(settings_list):

                if string.startswith("-"):

                    setting_name = string.split("--")[1] if string.startswith("--") else string.split("-")[1]
                    settings[setting_name] = True
                    last_was_optional = True

                elif last_was_optional:

                    last_setting_name = settings.keys()[-1]
                    settings[last_setting_name] = strings.unquote(string)
                    last_was_optional = False

                else:

                    settings[str(npositional)] = strings.unquote(string)
                    last_was_optional = False
                    npositional += 1

            # Get description
            description = self.get_description(command_name)

            # Set empty input dict
            input_dict = dict()

            # Command, description, settings, input_dict, cwd, finish
            # command, description, settings, input_dict
            command = Command(command_name, description, settings, input_dict, pts_settings=pts_settings)

            # Add the command
            commands.append(command)

        # Return
        return commands

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        if self.path is None: raise ValueError("The path is undefined")

        # Save
        self.saveto(self.path)

# -----------------------------------------------------------------
