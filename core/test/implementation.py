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
from ..tools.logging import log
from ..basics.configuration import DictConfigurationSetter
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

tables = introspection.get_arguments_tables()

# -----------------------------------------------------------------

class TestImplementation(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(TestImplementation, self).__init__(config, unlisted=True)

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

    def run_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get command properties
        the_command = command.command
        description = command.description
        settings_dict = command.settings
        input_dict = command.input_dict
        cwd = command.cwd

        # Find match in the tables of configurable classes
        match = introspection.resolve_command_tables(the_command, tables)

        # Get info
        module_path = match.module_path
        class_name = match.class_name
        configuration_module_path = match.configuration_module_path

        # Get the class
        cls = introspection.get_class(module_path, class_name)

        # Determine the output path
        output_path = fs.absolute_path(fs.join(self.path, cwd))

        ###

        # Debugging
        log.debug("Setting configuration for command '" + the_command + "' ...")

        # Change working directory
        fs.change_cwd(output_path)

        # Get the configuration definition
        configuration_module = importlib.import_module(configuration_module_path)
        definition = getattr(configuration_module, "definition")

        # Parse the configuration
        setter = DictConfigurationSetter(settings_dict, the_command, description=description)
        config = setter.run(definition)

        # Set working directory (output directory)
        config.path = output_path

        # Create the class instance, configure it with the configuration settings
        inst = cls(config)

        # Inform the user
        log.info("Executing command '" + the_command + "': " + description + " ...")

        # Run with input
        inst.run(**input_dict)

        # Return the instance
        return inst

# -----------------------------------------------------------------
