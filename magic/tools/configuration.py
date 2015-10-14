#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
from config import Config
from config import Mapping

# *****************************************************************

def open(filepath, default=None):

    """
    This function ...
    :param filepath:
    :param default:
    :return:
    """

    # If a default configuration file is not given, open the filepath
    if default is None: return Config(file(filepath))

    # Else, load the default configuration and adapt the properties according to
    # the specified configuration file
    else:

        # Open the default configuration
        config = Config(file(default))

        # Open the configuration file
        user_config = Config(file(filepath))

        adjust(config, user_config)

        return config

# *****************************************************************

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

            if not key in config: config[key] = Mapping()

            adjust(config[key], user_config[key])

        # Adapt the value of the property in the configuration to be equal to the value in the user configuration
        else: config[key] = user_config[key]

# *****************************************************************



