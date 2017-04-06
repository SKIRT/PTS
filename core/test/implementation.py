#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.implementation Contains the TestImplementation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import importlib
from collections import OrderedDict
from abc import abstractmethod, ABCMeta

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import introspection
from ..basics.configurable import Configurable
from ..launch.pts import launch_local, launch_remote

# -----------------------------------------------------------------

tables = introspection.get_arguments_tables()

# -----------------------------------------------------------------

class TestImplementation(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(TestImplementation, self).__init__(config, interactive, unlisted=True)

        # The test path
        self.path = None

        # The runnable components
        self.components = OrderedDict()
        self.input_dicts = dict()

    # -----------------------------------------------------------------

    @abstractmethod
    def run(self, **kwargs):

        """
        This function ...
        """
            
        pass

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        self.path = kwargs.pop("path")

    # -----------------------------------------------------------------

    def run_command(self, command, remote=None):

        """
        This function ...
        :param command:
        :param remote:
        :return:
        """

        # Get command properties
        the_command = command.command
        description = command.description
        settings_dict = command.settings
        input_dict = command.input_dict
        cwd = command.cwd

        # Set the cwd if not specified
        if not command.cwd_specified: cwd = self.path

        # Determine the output path
        output_path = fs.absolute_path(fs.join(self.path, cwd))

        # Launch locally or remotely
        if remote is not None: launch_remote(remote, the_command, settings_dict, input_dict)
        else: launch_local(the_command, settings_dict, input_dict, description=description, cwd=output_path)

        # Find match in the tables of configurable classes
        #match = introspection.resolve_command_tables(the_command, tables)

        # Get info
        #module_path = match.module_path
        #class_name = match.class_name
        #configuration_module_path = match.configuration_module_path

        # Get the class
        #cls = introspection.get_class(module_path, class_name)

        # Create the configuration
        #config = create_configuration_passive(the_command, class_name, configuration_module_path, settings_dict, description)

        # Change working directory
        #fs.change_cwd(output_path)

        # Set working directory (output directory)
        #config.path = output_path

        # Create the class instance, configure it with the configuration settings
        #inst = cls(config)

        # Inform the user
        #log.info("Executing command '" + the_command + "': " + description + " ...")

        # Run with input
        #inst.run(**input_dict)

        # Return the instance
        #return inst

# -----------------------------------------------------------------
