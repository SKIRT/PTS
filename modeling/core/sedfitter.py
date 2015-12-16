#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.sedfitter Contains the SEDFitter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.launch.remotelauncher import SkirtRemoteLauncher

# -----------------------------------------------------------------

class SEDFitter(Configurable):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(SEDFitter, self).__init__(config, "modeling")

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup()

# -----------------------------------------------------------------
