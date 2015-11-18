#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import inspect
from config import Config, Mapping

# -----------------------------------------------------------------

def special(self, item):

    """
    This function ...
    :param item:
    :return:
    """

    try:
        return self.__getitem__(item)
    except AttributeError:
        self[item] = Mapping()
        return self[item]

# Replace the __getattr__ function
Config.__getattr__ = special
Mapping.__getattr__ = special

# -----------------------------------------------------------------

def set(classname, config=None):

    """
    This function ...
    :param classname:
    :param config:
    :return:
    """

    # Determine the path to the default configuration file
    directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
    default_config = os.path.join(directory, "config", classname + ".cfg")

    # Open the default configuration if no configuration file is specified, otherwise adjust the default
    # settings according to the user defined configuration file
    if config is None: return open(default_config)
    else: return open(config, default_config)

# -----------------------------------------------------------------

def open(config, default_config=None):

    """
    This function ...
    :param filepath:
    :param default:
    :return:
    """

    # Open the config file
    if isinstance(config, basestring): config = Config(file(config))

    # If a default configuration file is not given, return the opened config file
    if default_config is None: return config
    else:

        # Open the default config file
        if isinstance(default_config, basestring): default_config = Config(file(default_config))

        # Adjust the default configuration according to the user-defined configuration
        adjust(default_config, config)

        # Return the adjusted default configuration
        return default_config

# -----------------------------------------------------------------

def adjust(config, user_config):

    """
    This function ...
    :param config:
    :param user_config:
    :return:
    """

    # Loop over all keys in the configuration
    for key in user_config:

        # If the property is a mapping (consists of more properties), call this function recursively
        if isinstance(user_config[key], Mapping):

            # If the Mapping with the name key does not exist in the config, add an empty mapping
            if not key in config: config[key] = Mapping()


            adjust(config[key], user_config[key])

        # Adapt the value of the property in the configuration to be equal to the value in the user configuration
        else: config[key] = user_config[key]

# -----------------------------------------------------------------



