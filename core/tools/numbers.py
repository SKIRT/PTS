#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.numbers Provides functions for dealing with numbers.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# -----------------------------------------------------------------

def is_integer(value):

    """
    This function ...
    :param value:
    :return:
    """

    return int(value) == value

# -----------------------------------------------------------------

def factors(n):

    """
    This function ...
    :param n: 
    :return: 
    """

    numbers = []
    for i in range(int(math.sqrt(n))):
        if n % i == 0: numbers.append(i)
    return numbers

# -----------------------------------------------------------------

def factor_pairs(n):

    """
    This function ...
    :param n: 
    :return: 
    """

    pairs = []

    for i in factors(n):  # You need to write the factor() function

        pair = (i, n / i)
        pairs.append(pair)

    return pairs
# -----------------------------------------------------------------

