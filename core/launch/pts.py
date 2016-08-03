#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.pts Contains the PTSRemoteLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import importlib

# Import the relevant PTS classes and modules
from ..basics.remote import Remote
from ...core.tools import introspection
from ...core.tools.logging import log
from ...core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter, InteractiveConfigurationSetter, FileConfigurationSetter

# -----------------------------------------------------------------

class PTSRemoteLauncher(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # -- Attributes --

        # Create the remote execution context
        self.remote = Remote()

    # -----------------------------------------------------------------

    def setup(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Setup the remote
        self.remote.setup(host_id)

    # -----------------------------------------------------------------

    def run(self, pts_command, input):

        """
        This function ...
        :param pts_command:
        :param input:
        :return:
        """

        # Resolve the PTS command
        subproject, command_name, description, class_name, configuration_module_path, configuration_method = find_match(pts_command)


        ## GET THE CONFIGURATION DEFINITION

        try:
            configuration_module = importlib.import_module(configuration_module_path)
            # has_configuration = True
            definition = getattr(configuration_module, "definition")
        except ImportError:
            log.warning("No configuration definition found for the " + class_name + " class")
            # has_configuration = False
            definition = ConfigurationDefinition()  # Create new configuration definition

        ## CREATE THE CONFIGURATION

        # Create the configuration setter
        if configuration_method == "interactive":
            setter = InteractiveConfigurationSetter(command_name, description, log_path="log")
        elif configuration_method == "arguments":
            setter = ArgumentConfigurationSetter(command_name, description, log_path="log")
        elif configuration_method.startswith("file"):
            configuration_filepath = configuration_method.split(":")[1]
            setter = FileConfigurationSetter(configuration_filepath, command_name, description, log_path="log")
        else:
            raise ValueError("Invalid configuration method: " + configuration_method)

        # Create the configuration from the definition and from reading the command line arguments
        config = setter.run(definition)


        ###


        # Exact command name
        exact_command_name = subproject + "/" + command_name

        # Inform the user about starting the PTS command remotely
        log.start("Starting " + exact_command_name + " on remote host " + self.remote.host_id + " ...")

        # Run PTS remotely
        #task = remote.run_pts(exact_command_name, config, keep_remote_temp=True)

        # Return the output dictionary
        return output_dict

# -----------------------------------------------------------------

def find_match(pts_command):

    """
    This function ...
    :param pts_command:
    :return:
    """

    # Find possible PTS commands
    scripts = introspection.get_scripts()
    tables = introspection.get_arguments_tables()

    # Find matches
    matches = introspection.find_matches_scripts(pts_command, scripts)
    table_matches = introspection.find_matches_tables(pts_command, tables)

    # No match
    if len(matches) + len(table_matches) == 0:

        #log.error("There is no match")
        #show_all_available(scripts, tables)
        #exit()
        raise ValueError("There is no match with a PTS command for '" + pts_command + "'")

    # If there is a unique match in an existing script, raise error
    elif len(matches) == 1 and len(table_matches) == 0: raise ValueError("This do command cannot be executed remotely")

    # If there is an unique match in a table
    elif len(table_matches) == 1 and len(matches) == 0:

        subproject, index = table_matches[0]
        command_name = tables[subproject]["Command"][index]
        description = tables[subproject]["Description"][index]
        class_path_relative = tables[subproject]["Path"][index]
        class_path = "pts." + subproject + "." + class_path_relative
        module_path, class_name = class_path.rsplit('.', 1)

        configuration_method_table = tables[subproject]["Configuration method"][index]

    # Show possible matches if there are more than just one
    else:

        #show_possible_matches(matches, table_matches, tables)
        raise ValueError("The PTS command '" + pts_command + "' is ambigious")

    configuration_name = tables[subproject]["Configuration"][index]
    if configuration_name == "--": configuration_name = command_name
    configuration_module_path = "pts." + subproject + ".config." + configuration_name

    # Return
    return subproject, command_name, description, class_name, configuration_module_path, configuration_method_table

# -----------------------------------------------------------------
