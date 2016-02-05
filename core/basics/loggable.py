#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.loggable Contains the Loggable class, a class representing all classes that
#  are capable of logging.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import old_logging

# -----------------------------------------------------------------

class Loggable(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :param config:
        :param subpackage:
        :return:
        """

        # -- Attributes --

        # The logger
        self.log = None

    # -----------------------------------------------------------------

    def setup(self, level="INFO", path=None):
        
        """
        This function ...
        :param level:
        :param path:
        """
    
        # Create the logger
        self.log = old_logging.new_log(self.name, level)
        if path is not None: old_logging.link_file_log(self.log, path, level)

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        name = type(self).__name__.lower()
        if "plotter" in name: return "plotter"
        else: return name

# -----------------------------------------------------------------
