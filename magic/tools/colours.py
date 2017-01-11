#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.colours Provides functions related to colours.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..core.frame import Frame

# -----------------------------------------------------------------

def calculate_colour(flux_a, flux_b):

    """
    This function ...
    :param flux_a:
    :param flux_b:
    :return:
    """

    return -2.5 * np.log10(flux_a / flux_b)

# -----------------------------------------------------------------

def make_colour_map(frame_a, frame_b):

    """
    This function ...
    :param frame_a:
    :param frame_b:
    :return:
    """

    return Frame(-2.5 * np.log10(frame_a / frame_b), wcs=frame_a.wcs)

# -----------------------------------------------------------------
