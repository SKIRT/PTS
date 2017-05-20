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
import collections
try:
    HAS_NP = True
    import numpy as np
except ImportError: HAS_NP = False

# -----------------------------------------------------------------

# IMPLEMENT INTELLIGENT META-PROGRAMMING TYPING THINGS:
# SEE: https://www.youtube.com/watch?v=sPiWg5jSoZI (around 1:10 ->  1:30 etc.) !!

# -----------------------------------------------------------------

if HAS_NP: boolean_types = [bool, np.bool]
else: boolean_types = [bool]

# -----------------------------------------------------------------

def is_boolean_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in boolean_types:
        if type(value) == tp: return True
    return False

# -----------------------------------------------------------------

if HAS_NP: integer_types = [int, np.int32, np.int64, np.uint32, np.uint64]
else: integer_types = [int]

# -----------------------------------------------------------------

def is_integer_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in integer_types:
        if type(value) == tp: return True
    return False

# -----------------------------------------------------------------

if HAS_NP: real_types = [float, np.float32, np.float64]
else: real_types = [float]

# -----------------------------------------------------------------

def is_real_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in real_types:
        if type(value) == tp: return True
    return False

# -----------------------------------------------------------------

if HAS_NP: string_types = [basestring, str, np.string_]
else: string_types = [basestring, str]

# -----------------------------------------------------------------

def is_string_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in string_types:

        # Literal
        if type(value) == tp: return True

        # Derived from
        if isinstance(value, tp): return True

    return False

# -----------------------------------------------------------------

def is_dictionary(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    return isinstance(value, dict)

# -----------------------------------------------------------------

def is_dictionary_of_dictionaries(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    if not is_dictionary(value): raise ValueError("Not a dictionary")
    if len(value) == 0: raise ValueError("Empty dictionary")

    # Loop over the keys
    for key in value:
        if not is_dictionary(value[key]): return False
    return True

# -----------------------------------------------------------------

def is_sequence(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    return isinstance(value, collections.Iterable)

# -----------------------------------------------------------------

def is_string_sequence(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    if not is_sequence(value): return False
    for element in value:
        if not is_string_type(element): return False
    return True

# -----------------------------------------------------------------
