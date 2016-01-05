#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.configurable Contains the Configurable class, a class for representing classes that can be
#  configured with a configuration file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from .loggable import Loggable
from ..tools import configuration

# -----------------------------------------------------------------

class Configurable(Loggable):

    """
    This class ...
    """

    def __init__(self, config, subpackage):

        """
        The constructor ...
        :param config:
        :param subpackage:
        :return:
        """

        # Call the constructor of the base class
        super(Configurable, self).__init__()

        # -- Attributes --

        # Set the configuration object
        self.config = configuration.set(subpackage, self.name, config)

    # -----------------------------------------------------------------

    def setup(self):
        
        """
        This function ...
        """

        # Determine the full path to the log file
        path = self.full_output_path(self.config.logging.path)

        # Call the setup function of the Loggable base class
        super(Configurable, self).setup(self.config.logging.level, path)

    # -----------------------------------------------------------------

    def full_input_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if os.path.isabs(name): return name
        elif "input_path" in self.config and self.config.input_path is not None: return os.path.join(self.config.input_path, name)
        else: return os.path.join(os.getcwd(), name)

    # -----------------------------------------------------------------

    def full_output_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if os.path.isabs(name): return name
        elif "output_path" in self.config and self.config.output_path is not None: return os.path.join(self.config.output_path, name)
        else: return os.path.join(os.getcwd(), name)

# -----------------------------------------------------------------
