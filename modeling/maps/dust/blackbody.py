#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.blackbody Contains the BlackBodyDustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log

# -----------------------------------------------------------------

class BlackBodyDustMapMaker(object):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(BlackBodyDustMapMaker, self).__init__()

        # -- Attributes --

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # ...
        self.make_map()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
