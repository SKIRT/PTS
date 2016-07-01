#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.buat Contains the BuatDustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log

# -----------------------------------------------------------------

class BuatDustMapMaker(object):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(BuatDustMapMaker, self).__init__()

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

        #
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

        # Buat et al: #a_fuv_buat = (-0.0333*x3) + (0.3522*x2) + (1.1960*x) + 0.4967

        a_fuv = -0.0333 * y**3 + 0.3522 * y**2 + 1.1960 * y + 0.4967

        # where y = F_dust / F_FUV


# -----------------------------------------------------------------
