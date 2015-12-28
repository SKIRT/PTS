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

        # Call the setup function of the Loggable base class
        super(Configurable, self).setup(self.config.logging.level, self.config.logging.path)

# -----------------------------------------------------------------
