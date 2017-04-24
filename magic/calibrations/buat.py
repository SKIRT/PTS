#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.calibrations.buat Contains the BuatAttenuationCalibration class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class BuatAttenuationCalibration(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        # Set the buat parameters
        self.data = dict()
        self.data["NUV"] = (-0.0495, 0.4718, 0.8998, 0.2269)
        self.data["FUV"] = (-0.0333, 0.3522, 1.1960, 0.4967)

    # -----------------------------------------------------------------

    def get_nuv_parameters(self):

        """
        This function ...
        :return: 
        """

        return self.data["NUV"]

    # -----------------------------------------------------------------

    def get_fuv_parameters(self):

        """
        This function ...
        :return: 
        """

        return self.data["FUV"]

# -----------------------------------------------------------------
