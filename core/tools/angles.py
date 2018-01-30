#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.angles Provides functions related to angles.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

def quadrant(angle):

    """
    This function ...
    :param angle:
    :return:
    """

    if -180 <= angle.to("deg").value < -90.: return 3
    elif -90. <= angle.to("deg").value < 0.0: return 4
    elif 0.0 <= angle.to("deg").value < 90.: return 1
    elif 90. <= angle.to("deg").value < 180.: return 2
    elif 180. <= angle.to("deg").value < 270.: return 3
    elif 270. <= angle.to("deg").value <= 360.: return 4
    else: raise ValueError("Failed to determine quadrant for " + str(angle))

# -----------------------------------------------------------------
