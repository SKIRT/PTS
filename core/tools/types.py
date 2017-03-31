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

class Type(object):

    """
    This class ...
    """

# -----------------------------------------------------------------

class Boolean(Type):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This fucntion ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class Real(Type):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class Integer(Type):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        The constructor ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class NonNegativeInteger(Integer):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        The constructor ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class NonPositiveInteger(Integer):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        The constructor ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class PositiveInteger(NonNegativeInteger):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        The constructor ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class NegativeInteger(NonPositiveInteger):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        The constructor ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class NonNegativeReal(Real):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class NonPositiveReal(Real):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class PositiveReal(NonNegativeReal):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class NegativeReal(NonPositiveReal):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class String(Type):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class NoSpacesString(String):

    """
    This class ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class Path(String):

    """
    This function ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class FilePath(Path):

    """
    This function ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

# -----------------------------------------------------------------

class DirectoryPath(Path):

    """
    This fucntion ...
    """

    def __new__(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

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

if HAS_NP: string_types = [basestring, np.string_]
else: string_types = [basestring]

# -----------------------------------------------------------------

def is_string_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in string_types:
        if type(value) == tp: return True
    return False

# -----------------------------------------------------------------
