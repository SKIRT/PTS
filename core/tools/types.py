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
try:
    HAS_NP = True
    import numpy as np
except ImportError: HAS_NP = False

# -----------------------------------------------------------------

def is_boolean_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    if HAS_NP: return isinstance(value, bool) or isinstance(value, np.bool)
    else: return isinstance(value, bool)

# -----------------------------------------------------------------

def is_integer_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    if HAS_NP: return isinstance(value, int) or isinstance(value, np.int32) or isinstance(value, np.int64) or isinstance(value, np.uint32) or isinstance(value, np.uint64)
    else: return isinstance(value, int)

# -----------------------------------------------------------------

def is_real_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    if HAS_NP: return isinstance(value, float) or isinstance(value, np.float32) or isinstance(value, np.float64)
    else: return isinstance(value, float)

# -----------------------------------------------------------------

def is_string_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    if HAS_NP: return isinstance(value, basestring) or isinstance(value, np.string_)
    else: return isinstance(value, basestring)

# -----------------------------------------------------------------
