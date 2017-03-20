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
