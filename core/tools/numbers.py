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
import random
import numpy as np
from itertools import cycle

# Import the relevant PTS classes and modules
from . import sequences

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

def derivatives(x, y):

    """
    This function calculates the derivative values using finite differences
    :param x: 
    :param y: 
    :return: 
    """

    new_x = []
    for i in range(len(x)-1):
        between = 0.5 * (x[i] + x[i+1])
        new_x.append(between)

    # Calculate derivatives
    new_y = np.diff(y) / np.diff(x)

    # Return x, y finite differences derivative data
    return new_x, new_y

# -----------------------------------------------------------------

def test_division_in_n_dimensions(n):

    """
    This function ...
    :param n: 
    :return: 
    """

    amount = 50

    # Generate random integer numbers between 10 and 1000
    for _ in range(amount):

        number = random.randint(10, 1000)
        factors = divide_in_n_dimensions(number, n)

        result = sequences.multiply_all_integers(factors)

        print(number, result)

# -----------------------------------------------------------------

def divide_in_n_dimensions(number, n, sampled_most=None, weights=None):

    """
    This function ...
    :param number: 
    :param n: 
    :param sampled_most:
    :param weights:
    :return: 
    """

    # Check sampled_most and weights argument
    if sampled_most is not None and weights is not None: raise ValueError("Either define 'sampled_most' or 'weights'")

    from . import types
    if not types.is_integer_type(number): raise ValueError("Number must be integer")

    result = number**(1./n)
    result = int(math.ceil(result))

    the_factors = [result] * n

    # Multiply with weights
    if weights is not None: the_factors = [int(math.ceil(the_factors[index] * weights[index])) for index in range(len(the_factors))]

    #print(factors)

    # Create iterator of indices of labels to decrease the value
    if sampled_most is not None: lst = [index for index in range(n) if not sampled_most[index]]
    else: lst = range(n)
    indices = cycle(lst)

    # Lower some of the factors till the result is as small as possible, but still equal to or greater than the initial number
    previous_factors = None
    while True:

        #print("factors", factors)
        product = sequences.multiply_all_integers(the_factors)

        if product < number: return previous_factors
        else:
            previous_factors = the_factors[:] # copy
            index = indices.next()
            #print("index", index)
            # Lower one of the factors
            the_factors[index] = the_factors[index] - 1

# -----------------------------------------------------------------
