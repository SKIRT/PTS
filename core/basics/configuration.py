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
import argparse
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .map import Map
from ..tools import parsing
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class Configuration(object):

    """
    This function ...
    """

    def __init__(self, add_logging=True, add_cwd=True, prefix=None, log_path=None):

        """
        The constructor ...
        """

        # Create the command-line parser
        self.parser = argparse.ArgumentParser()

        # Dictionary of sections
        self.sections = OrderedDict()

        # Dictionary of fixed parameters
        self.fixed = dict()

        # List of argument names for the settings of the class instance the configuration is converted into
        self.to_instance = []

        # The prefix (for parent configurations)
        self.prefix = prefix

        # The parsed arguments
        self.arguments = None

        # Add logging options
        if add_logging:

            # Log path to absolute path
            log_path = fs.absolute(log_path) if log_path is not None else fs.cwd()
            self.add_fixed("log_path", log_path)
            self.add_flag("debug", "enable debug output", False)
            self.add_flag("report", "write a report file", False)

        # Add the path to the current working directory
        if add_cwd: self.add_fixed("path", fs.cwd())

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, add_logging=True, add_cwd=True):

        """
        This function ...
        """

        # Create new configuration
        config = cls(add_logging, add_cwd)

        # Add options from file ...

        # Return the configuration object
        return config

    # -----------------------------------------------------------------

    def add_section(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        self.sections[name] = Configuration(add_logging=False, prefix=name)

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

    def add_required(self, name, user_type, description, to_instance=True, choices=None):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param to_instance:
        :param choices:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Add prefix
        if self.prefix is not None: name = self.prefix + "/" + name

        # Add the argument
        self.parser.add_argument(name, type=real_type, help=description, choices=choices)

        # To instance
        if to_instance: self.to_instance.append(name)

    # -----------------------------------------------------------------

    def add_positional_optional(self, name, user_type, description, default=None, to_instance=True, choices=None, convert_default=False):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param to_instance:
        :param convert_default:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if convert_default: default = get_real_value(default, user_type)

        # Add prefix
        if self.prefix is not None: name = self.prefix + "/" + name

        # Add the argument
        self.parser.add_argument(name, type=real_type, help=description, default=default, nargs='?', choices=choices)

        # To instance
        if to_instance: self.to_instance.append(name)

    # -----------------------------------------------------------------

    def add_optional(self, name, user_type, description, default=None, to_instance=True, letter=None, convert_default=False):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param to_instance:
        :param letter:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if convert_default: default = get_real_value(default, user_type)

        # Add prefix
        if self.prefix is not None: name = self.prefix + "/" + name

        # Add the argument
        if letter is None: self.parser.add_argument("--" + name, type=real_type, help=description, default=default)
        else: self.parser.add_argument("-" + letter, "--" + name, type=real_type, help=description, default=default)

        # To instance
        if to_instance: self.to_instance.append(name)

    # -----------------------------------------------------------------

    def add_flag(self, name, description, to_instance=True, letter=None):

        """
        This function ...
        :param name:
        :param description:
        :param to_instance:
        :return:
        """

        # Add prefix
        if self.prefix is not None: name = self.prefix + "/" + name

        # Add the argument
        if letter is None: self.parser.add_argument("--" + name, action="store_true", help=description)
        else: self.parser.add_argument("-" + letter, "--" + name, action="store_true", help=description)

        # To instance
        if to_instance: self.to_instance.append(name)

    # -----------------------------------------------------------------

    def read(self):

        """
        This function ...
        :return:
        """

        # Let the parser parse
        self.arguments = self.parser.parse_args()

        # Recursive for the sections
        for name in self.sections: self.sections[name].read()

    # -----------------------------------------------------------------

    def get_settings(self):

        """
        This function ...
        :return:
        """

        # Create the settings Map
        settings = Map()

        # Add fixed
        for name in self.fixed: settings[name] = self.fixed[name]

        # Add base settings
        for argument in vars(self.arguments): settings[argument] = getattr(self.arguments, argument)

        # Add the configuration settings of the various sections
        for name in self.sections:

            # Create a map for the settings
            settings[name] = Map()

            for argument in vars(self.sections[name].arguments): settings[name][argument] = getattr(self.sections[name].arguments, argument)

        # Return the settings
        return settings

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
