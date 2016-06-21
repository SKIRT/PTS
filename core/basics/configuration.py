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

# Import the relevant PTS classes and modules
from .map import Map
from ..tools import parsing
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class Configuration(object):

    """
    This function ...
    """

    def __init__(self, add_logging=True, add_cwd=True):

        """
        The constructor ...
        """

        # Create the command-line parser
        self.parser = argparse.ArgumentParser()

        # Dictionary of fixed parameters
        self.fixed = dict()

        # List of argument names for the settings of the class instance the configuration is converted into
        self.to_instance = []

        # The parsed arguments
        self.arguments = None

        # Add logging options
        if add_logging:

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

    def add_fixed(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        self.fixed[name] = value

    # -----------------------------------------------------------------

    def add_required(self, name, user_type, description, to_instance=True):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param to_instance:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Add the argument
        self.parser.add_argument(name, type=real_type, help=description)

        # To instance
        if to_instance: self.to_instance.append(name)

    # -----------------------------------------------------------------

    def add_optional(self, name, user_type, description, default, to_instance=True):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param to_instance:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        default = get_real_value(default, user_type)

        # Add the argument
        self.parser.add_argument("--" + name, type=real_type, help=description, default=default)

        # To instance
        if to_instance: self.to_instance.append(name)

    # -----------------------------------------------------------------

    def add_flag(self, name, description, to_instance=True):

        """
        This function ...
        :param name:
        :param description:
        :param to_instance:
        :return:
        """

        # Add the argument
        self.parser.add_argument("--" + name, action="store_true", help=description)

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

    # -----------------------------------------------------------------

    def get_settings(self):

        """
        This function ...
        :return:
        """

        pass

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
