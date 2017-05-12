#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

# ratio of transitions to transversions is typically around 2

from ...core.tools import numbers

# -----------------------------------------------------------------

# 00, 01, 10, 11
encoding = ["A", "C", "T", "G"]

#  A -> G (purines) and C -> T (pyrimidines): transitions
#  other: tranversions

# -----------------------------------------------------------------

def bits_to_nucleobase(combination):

    """
    This function ...
    :param combination: 
    :return: 
    """

    if combination == [0, 0]: return encoding[0]
    elif combination == [0, 1]: return encoding[1]
    elif combination == [1, 0]: return encoding[2]
    elif combination == [1, 1]: return encoding[3]
    else: raise ValueError("Invalid combination: " + str(combination))

# -----------------------------------------------------------------

def binary_string_to_nucleobase_string(binary_string):

    """
    This function ...
    :param binary_string: 
    :return: 
    """

    if not numbers.is_even(len(binary_string)): raise ValueError("Binary string length must be even")

    bases = []

    # loop over the pairs
    for index in range(len(binary_string)/2):

        index0 = 2*index
        index1 = 2*index + 1

        combination = [binary_string[index0], binary_string[index1]]

        base = bits_to_nucleobase(combination)
        bases.append(base)

    # Return the base string
    return bases

# -----------------------------------------------------------------