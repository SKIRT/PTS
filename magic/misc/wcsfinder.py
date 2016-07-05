#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.wcsfinder Contains the WCSFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import OldConfigurable

# -----------------------------------------------------------------

class WCSFinder(OldConfigurable):

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
        super(WCSFinder, self).__init__(config, "magic")

    # -----------------------------------------------------------------

    def run(self, frame):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()
        
# -----------------------------------------------------------------
