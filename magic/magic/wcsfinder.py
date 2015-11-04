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

# Import Astromagic modules
from ..tools import configuration

# *****************************************************************

class WCSFinder(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "wcsfinder.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

    # *****************************************************************

    def run(self, frame):

        """
        This function ...
        :return:
        """

        pass
        
# *****************************************************************
