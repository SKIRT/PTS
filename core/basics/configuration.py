#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.configuration Contains the configuration class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys
import argparse
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .map import Map
from ..tools import parsing
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class ConfigurationDefinition(object):

    """
    This function ...
    """

    def __init__(self, prefix=None):

        """
        This function ...
        :param prefix:
        """

        # Prefix
        self.prefix = prefix

        # Dictionary of sections
        self.sections = OrderedDict()

        # Dictionary of fixed parameters
        self.fixed = OrderedDict()

        # Dictionary of required parameters
        self.required = OrderedDict()

        # Dictionary of positional optional parameters
        self.pos_optional = OrderedDict()

        # Dictionary of optional parameters
        self.optional = OrderedDict()

        # Dictionary of flags
        self.flags = OrderedDict()

    # -----------------------------------------------------------------

    def set_arguments(self, parser):

        """
        This function ...
        :param parser:
        :return:
        """

        # Required
        for name in self.required:

            real_type = self.required[name][0]
            description = self.required[name][1]
            choices = self.required[name][2]

            # Add prefix
            if self.prefix is not None: name = self.prefix + "/" + name

            # Add argument to argument parser
            parser.add_argument(name, type=real_type, help=description, choices=choices)

        # Positional optional
        for name in self.pos_optional:

            real_type = self.pos_optional[name][0]
            description = self.pos_optional[name][1]
            default = self.pos_optional[name][2]
            choices = self.pos_optional[name][3]

            # Add prefix
            if self.prefix is not None: name = self.prefix + "/" + name

            # Add argument to argument parser
            parser.add_argument(name, type=real_type, help=description, default=default, nargs='?', choices=choices)

        # Optional
        for name in self.optional:

            # (real_type, description, default, letter)
            real_type = self.optional[name][0]
            description = self.optional[name][1]
            default = self.optional[name][2]
            letter = self.optional[name][3]

            # Add prefix
            if self.prefix is not None: name = self.prefix + "/" + name

            # Add the argument
            if letter is None: parser.add_argument("--" + name, type=real_type, help=description, default=default)
            else: parser.add_argument("-" + letter, "--" + name, type=real_type, help=description, default=default)

        # Flag
        for name in self.flags:

            # (description, letter)
            description = self.flags[name][0]
            letter = self.flags[name][1]

            # Add prefix
            if self.prefix is not None: name = self.prefix + "/" + name

            # Add the argument
            if letter is None: parser.add_argument("--" + name, action="store_true", help=description)
            else: parser.add_argument("-" + letter, "--" + name, action="store_true", help=description)

        # Add arguments of sections
        for section_name in self.sections: self.sections[section_name].set_arguments(parser)

    # -----------------------------------------------------------------

    def get_settings(self, settings, arguments):

        """
        This function ...
        :param settings:
        :param arguments:
        :return:
        """

        # Add fixed
        for name in self.fixed: settings[name] = self.fixed[name]

        # Add required
        for name in self.required:

            if self.prefix is not None: argument_name = self.prefix + "/" + name
            else: argument_name = name

            settings[name] = getattr(arguments, argument_name)

        # Add positional optional
        for name in self.pos_optional:

            if self.prefix is not None: argument_name = self.prefix + "/" + name
            else: argument_name = name

            settings[name] = getattr(arguments, argument_name)

        # Add optional
        for name in self.optional:

            if self.prefix is not None: argument_name = self.prefix + "/" + name
            else: argument_name = name

            settings[name] = getattr(arguments, argument_name)

        # Add flags
        for name in self.flags:

            if self.prefix is not None: argument_name = self.prefix + "/" + name
            else: argument_name = name

            settings[name] = getattr(arguments, argument_name)

        # Add the configuration settings of the various sections
        for name in self.sections:

            # Create a map for the settings
            settings[name] = Map()

            # Recursively add the settings
            self.sections[name].get_settings(settings[name], arguments)

    # -----------------------------------------------------------------

    def add_section(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Add the section
        self.sections[name] = ConfigurationDefinition(prefix=name)

    # -----------------------------------------------------------------

    def add_fixed(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        self.fixed[name] = value

    # -----------------------------------------------------------------

    def add_required(self, name, user_type, description, choices=None):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param choices:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Add
        self.required[name] = (real_type, description, choices)

    # -----------------------------------------------------------------

    def add_positional_optional(self, name, user_type, description, default=None, choices=None, convert_default=False):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param choices:
        :param convert_default:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if convert_default: default = get_real_value(default, user_type)

        # Add
        self.pos_optional[name] = (real_type, description, default, choices)

    # -----------------------------------------------------------------

    def add_optional(self, name, user_type, description, default=None, letter=None, convert_default=False):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param letter:
        :param convert_default:
        :return:
        """

        if self.prefix is not None and letter is not None: raise ValueError("Cannot assign letter argument for child configuration definition")

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if convert_default: default = get_real_value(default, user_type)

        # Add
        self.optional[name] = (real_type, description, default, letter)

    # -----------------------------------------------------------------

    def add_flag(self, name, description, letter=None):

        """
        This function ...
        :param name:
        :param description:
        :param letter:
        :return:
        """

        # Add
        self.flags[name] = (description, letter)

# -----------------------------------------------------------------

class ConfigurationReader(object):

    """
    This function ...
    """

    def __init__(self, name, description=None, add_logging=True, add_cwd=True, log_path=None):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param log_path:
        """

        # The configuration definition
        self.definition = None

        # Create the command-line parser
        self.parser = argparse.ArgumentParser(prog=name, description=description)

        # The parsed arguments
        self.arguments = None

        # ...
        self.add_logging = add_logging
        self.add_cwd = add_cwd
        self.log_path = log_path

        # The config
        self.config = None

    # -----------------------------------------------------------------

    @staticmethod
    def get_arguments():

        """
        This function ...
        :return:
        """

        return sys.argv[1:]

    # -----------------------------------------------------------------

    def read(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set definition
        self.definition = definition

        # Set logging
        self.set_logging_and_cwd()

        # Parse
        self.parse()

        # Create config
        self.create_config()

        # Return the config
        return self.config

    # -----------------------------------------------------------------

    def set_logging_and_cwd(self):

        """
        This function ...
        :return:
        """

        # Add logging options
        if self.add_logging:

            # Log path to absolute path
            log_path = fs.absolute(self.log_path) if self.log_path is not None else fs.cwd()
            self.definition.add_fixed("log_path", log_path)
            self.definition.add_flag("debug", "enable debug output", False)
            self.definition.add_flag("report", "write a report file", False)

        # Add the path to the current working directory
        if self.add_cwd: self.definition.add_fixed("path", fs.cwd())

    # -----------------------------------------------------------------

    def parse(self):

        """
        This function ...
        :return:
        """

        # Set arguments
        self.definition.set_arguments(self.parser)

        # Let the parser parse
        self.arguments = self.parser.parse_args()

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Create the settings Map
        self.config = Map()

        # Add the settings
        self.definition.get_settings(self.config, self.arguments)

# -----------------------------------------------------------------

def get_real_type(user_type):

    """
    This function ...
    :param user_type:
    :return:
    """

    if isinstance(user_type, basestring): return getattr(parsing, user_type)
    else: return user_type

# -----------------------------------------------------------------

def get_real_value(default, user_type):

    """
    This function ...
    :param default:
    :return:
    """

    return user_type(default)

# -----------------------------------------------------------------
