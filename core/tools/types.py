#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.types Checking types of objects.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# -----------------------------------------------------------------

def is_boolean_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    return isinstance(value, bool) or isinstance(value, np.bool)

# -----------------------------------------------------------------

def is_integer_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    return isinstance(value, int) or isinstance(value, np.int32) or isinstance(value, np.int64) or isinstance(value, np.uint32) or isinstance(value, np.uint64)

# -----------------------------------------------------------------

def is_real_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    return isinstance(value, float) or isinstance(value, np.float32) or isinstance(value, np.float64)

# -----------------------------------------------------------------
