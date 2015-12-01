#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# *****************************************************************

"""
This module ...
"""

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.basics import Configurable

# -----------------------------------------------------------------

class WCSFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(WCSFinder, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self, frame):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()
        
# -----------------------------------------------------------------
