#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module ...
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import logging

# Import the relevant PTS classes and modules
from pts.core.tools import configuration

# -----------------------------------------------------------------

class Configurable(object):

    """
    This class ...
    """

    def __init__(self, config):

        # Set the configuration object
        self.config = configuration.set(self.name, config)
        
        # Set the logger to None initially
        self.log = None

    # -----------------------------------------------------------------

    def setup(self):
        
        """
        This function ...
        """
    
        # Create the logger
        self.log = logging.new_log(self.name, self.config.logging.level)
        if self.config.logging.path is not None: logging.link_file_log(self.log, self.config.logging.path, self.config.logging.level)

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return type(self).__name__.lower()

# -----------------------------------------------------------------
