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
from math import log10, floor, ceil
from operator import mul

# Import the relevant PTS classes and modules
from . import sequences

# -----------------------------------------------------------------

zero = 0.
nan = float("nan")
inf = float("inf")
min_inf = float("-inf")

# -----------------------------------------------------------------

nan_values = [float("nan"), np.NaN, np.nan]
inf_values = [float("inf"), float("-inf"), np.Inf, -np.Inf, np.inf, -np.inf]

# -----------------------------------------------------------------

def is_multiple_of(number, other_number):
    return number % other_number == 0

# -----------------------------------------------------------------

def is_divisor_of(number, other_number):
    return other_number % number == 0

# -----------------------------------------------------------------

def is_multiple_or_divisor_of(number, other_number):
    return is_multiple_of(number, other_number) or is_divisor_of(number, other_number)

# -----------------------------------------------------------------

def equal_parts(number):

    """
    This function ...
    :param number:
    :return:
    """

    if not is_integer(number): raise ValueError("Not an integer number")

    if is_even(number): a = int(number/2)
    else: a = int(number/2)

    # Return
    return a, number - a

# -----------------------------------------------------------------

def is_number(value_or_string):

    """
    This function ...
    :param value_or_string:
    :return:
    """

    try:
        #print(value_or_string)
        float(value_or_string)
        return True
    except ValueError: return False

# -----------------------------------------------------------------

def find_numbers_in_string(string, left, right):

    """
    This function
    e.g. find numbers 0, 3, 5 in blabla[0]blabla[3]blabla[5] with left = "[" and right = "]"
    :param string:
    :param left:
    :param right:
    :return:
    """

    numbers = []

    for part in string.split(left):

        test = part.split(right)[0]
        if is_number(test): numbers.append(float(test))

    # Return the list of numbers
    return numbers

# -----------------------------------------------------------------

def is_close(value, other_value, rtol=1.e-5, atol=1.e-8):

    """
    This function ...
    :param value:
    :param other_value:
    :param rtol:
    :param atol:
    :return:
    """

    # Check whether values have unit
    if hasattr(value, "unit"):
        unit = value.unit
        if not hasattr(other_value, "unit"): raise ValueError("Both values must have a unit, or neither")
        value = value.value
        other_value = other_value.to(unit).value

    # Check values
    return np.isclose(value, other_value, rtol=rtol, atol=atol)

# -----------------------------------------------------------------

def is_close_to_zero(value, atol=1.e-8):
    return is_close(value, 0., atol=atol)

# -----------------------------------------------------------------

def is_even(integer):
    return integer % 2 == 0

# -----------------------------------------------------------------

def is_odd(integer):
    return integer % 2 != 0

# -----------------------------------------------------------------

def is_numpy_nan(value):
    try: return np.isnan(value)
    except TypeError: return False

# -----------------------------------------------------------------

def is_numpy_inf(value):
    try: return np.isinf(value)
    except TypeError: return False

# -----------------------------------------------------------------

def is_nan(value):
    return value in nan_values or is_numpy_nan(value)

# -----------------------------------------------------------------

def is_inf(value):
    return value in inf_values or is_numpy_inf(value)

# -----------------------------------------------------------------

def is_invalid(value):
    return is_nan(value) or is_inf(value)

# -----------------------------------------------------------------

def is_integer(value, absolute=True, rtol=1.e-5, atol=1.e-8):

    """
    This function ...
    :param value:
    :param absolute:
    :param rtol:
    :param atol:
    :return:
    """

    if isinstance(value, int): return True

    #print(value, int(round(value)), int(value))
    #print(value - round(value))

    # Special values
    if is_nan(value): return False
    if is_inf(value): return False

    #print(value, type(value))

    # Regular value
    if absolute: return int(value) == value
    else:
        #difference = value - round(value)
        return np.isclose(value, round(value), rtol=rtol, atol=atol)

# -----------------------------------------------------------------

def as_integer_check(value, absolute=True, rtol=1.e-5, atol=1.e-8):

    """
    This function ...
    :param value:
    :param absolute:
    :param rtol:
    :param atol:
    :return:
    """

    if not is_integer(value, absolute=absolute, rtol=rtol, atol=atol): raise ValueError("Not an integer")

    if absolute: return int(value)
    else: return int(round(value))

# -----------------------------------------------------------------

def as_integer_check_division(value_a, value_b, absolute=True, rtol=1.e-5, atol=1.e-8):

    """
    This function ...
    :param value_a:
    :param value_b:
    :param absolute:
    :param rtol:
    :param atol:
    :return:
    """

    return as_integer_check(value_a/value_b, absolute=absolute, rtol=rtol, atol=atol)

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

    # Check sampled_most and weights argument: NO, IT IS POSSIBLE FOR BOTH TO HAVE EFFECT
    #if sampled_most is not None and weights is not None: raise ValueError("Either define 'sampled_most' or 'weights'")

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

def round_to_1_significant_digit(x):

    """
    >>> round_to_1_significant_digit(0.0232)
    0.02
    >>> round_to_1_significant_digit(1234243)
    1000000.0
    >>> round_to_1_significant_digit(13)
    10.0
    >>> round_to_1_significant_digit(4)
    4.0
    >>> round_to_1_significant_digit(19)
    20.0
    :param x: 
    :return: 
    """

    return round(x, -int(floor(log10(abs(x)))))

# -----------------------------------------------------------------

def round_to_n_significant_digits(x, n):

    """
    >>> round_to_n_significant_digits(0.0232)
    0.023
    >>> round_to_n_significant_digits(0.0232, 1)
    0.02
    >>> round_to_n_significant_digits(1234243, 3)
    1230000.0
    :param x: 
    :param n: 
    :return: 
    """

    if x == 0: return 0
    return round(x, n-int(floor(log10(abs(x))))-1)

# -----------------------------------------------------------------

def integer_bit_length(integer):

    """
    This fucntion ...
    :param integer: 
    :return: 
    """

    return integer.bit_length()

# -----------------------------------------------------------------

def nbits_for_integer(integer):

    """
    This function ...
    :param integer: 
    :return: 
    """

    return integer_bit_length(integer)

# -----------------------------------------------------------------

def min_nbits_for_nintegers(nintegers):

    """
    This function ...
    :param nintegers: 
    :return: 
    """

    return nbits_for_integer(nintegers-1)

# -----------------------------------------------------------------

def max_integer_for_nbits(nbits):

    """
    This function ...
    :param nbits: 
    :return: 
    """

    return 2 ** nbits - 1

# -----------------------------------------------------------------

def nintegers_for_nbits(nbits):

    """
    This function ...
    :param nbits: 
    :return: 
    """

    return 2 ** nbits

# -----------------------------------------------------------------

def integer_to_binary(integer):

    """
    This function ...
    :param integer: 
    :return: 
    """

    from numpy import binary_repr
    return binary_repr(integer)

# -----------------------------------------------------------------

def binary_to_integer(binary):

    """
    This function ...
    :param binary: 
    :return: 
    """

    return int(str(binary), 2)

# -----------------------------------------------------------------

def integer_to_quaternary(integer):

    """
    This function ...
    :param integer: 
    :return: 
    """

    from numpy import base_repr
    return base_repr(integer, 4)

# -----------------------------------------------------------------

def quaternary_to_integer(quaternary):

    """
    This fucntion ...
    :param quaternary: 
    :return: 
    """

    return int(str(quaternary), 4)

# -----------------------------------------------------------------

def integer_to_binary_string(integer, nbits=None):

    """
    This function ...
    :param integer: 
    :param nbits:
    :return: 
    """

    binary = integer_to_binary(integer)
    return binary_to_binary_string(binary, nbits)

# -----------------------------------------------------------------

def binary_string_to_integer(binary_string):

    """
    This function ...
    :param binary_string: 
    :return: 
    """

    binary = binary_string_to_binary(binary_string)
    return binary_to_integer(binary)

# -----------------------------------------------------------------

# FROM:
# Traditional Techniques of Genetic Algorithms Applied to Floating-Point Chromosome Representations
# Leo Budin, Marin Golub, Andrea Budin

# -----------------------------------------------------------------

def float_to_binary(value, low, high, nbits):

    """
    This function ...
    :param value: 
    :param high:
    :param low:
    :param nbits:
    :return: 
    """

    # Set to floats
    value = float(value)
    high = float(high)
    low = float(low)

    scaled = (value - low) / (high - low) * 2**nbits
    integer = int(round(scaled))
    return integer_to_binary(integer)

# -----------------------------------------------------------------

def binary_to_float(binary, low, high, nbits):

    """
    This function ...
    :param binary: 
    :param high:
    :param low:
    :param nbits:
    :return: 
    """

    # Set to floats
    high = float(high)
    low = float(low)

    integer = binary_to_integer(binary)
    scaled = float(integer)

    value = low + scaled * 2**(-nbits) * (high - low)

    largest_error = 0.5 * (high - low) / nintegers_for_nbits(nbits)
    #print(largest_error)
    #scaled += largest_error

    value += 0.5 * largest_error

    return value

# -----------------------------------------------------------------

def binary_to_binary_string(binary, nbits=None):

    """
    This function ...
    :param binary: 
    :param nbits:
    :return: 
    """

    if nbits is None: characters = list(str(binary))
    else:
        string = str(binary)
        npadded = nbits - len(string)
        characters = list("0" * npadded + string)

    # Return the binary string as a list of integers
    return [int(character) for character in characters]

# -----------------------------------------------------------------

def binary_string_to_binary(binary_string):

    """
    This function ...
    :param binary_string: 
    :return: 
    """

    return int("".join(str(bit) for bit in binary_string))

# -----------------------------------------------------------------

def binary_string_to_float(binary_string, low, high, nbits):
    
    """
    This function ...
    :param binary_string: 
    :param low:
    :param high:
    :param nbits:
    :return: 
    """

    binary = binary_string_to_binary(binary_string)
    return binary_to_float(binary, low, high, nbits)

# -----------------------------------------------------------------

def float_to_binary_string(value, low, high, nbits):
    
    """
    This function ...
    :param value:
    :param low:
    :param high:
    :param nbits:
    :return: 
    """

    binary = float_to_binary(value, low, high, nbits)
    return binary_to_binary_string(binary, nbits)

# -----------------------------------------------------------------

# https://math.stackexchange.com/questions/1968416/number-of-significant-figures-when-going-from-base-10-to-binary
# First comment: I would think math.floor(log2(10)) = 3 significant figures in binary per significant figure in base 10

# -----------------------------------------------------------------

def binary_digits_for_significant_figures(nfigures):

    """
    This function ...
    :param nfigures: 
    :return: 
    """

    return int(math.floor(np.log2(10) * nfigures)) + 1

# -----------------------------------------------------------------

def nbits_for_ndigits(ndigits):

    """
    This fucntion ...
    :param ndigits: 
    :return: 
    """

    return binary_digits_for_significant_figures(ndigits)

# -----------------------------------------------------------------

def largest_error(nbits, low, high):

    """
    This function ...
    :param nbits: 
    :param low: 
    :param high: 
    :return: 
    """

    return 0.5 * (high - low) / nintegers_for_nbits(nbits)

# -----------------------------------------------------------------

def order_of_magnitude(number):

    """
    This function ...
    :param number: 
    :return: 
    """

    # TODO: doesn't work yet with 0.0XXYY !

    string = str(float(number))

    if "e" in string:

        power = int(string.split("e")[1])

        splitted = string.split("e")[0].split(".")
        first_string = splitted[0]

        return len(first_string) - 1 + power

    else:

        splitted = string.split(".")
        first_string = splitted[0]

        if len(splitted) == 1: second_string = ""
        elif len(splitted) == 2: second_string = splitted[1]
        else: raise RuntimeError("Something went wrong")

        return len(first_string) - 1

# -----------------------------------------------------------------

def rounding_error_order_of_magnitude(number, ndigits):

    """
    This function ...
    :param number: 
    :param ndigits: 
    :return: 
    """

    error = order_of_magnitude(number) - ndigits
    return error

# -----------------------------------------------------------------

def rounding_error(number, ndigits):

    """
    This function ...
    :param number: 
    :param ndigits: 
    :return: 
    """

    return 5. * 10**rounding_error_order_of_magnitude(number, ndigits)

# -----------------------------------------------------------------

def maximal_error(number, ndigits):

    """
    This fucntion ...
    :param number: 
    :param ndigits: 
    :return: 
    """

    return rounding_error(number, ndigits)

# -----------------------------------------------------------------

def minimal_nsteps(number, low, high, ndigits):

    """
    This function ...
    :param number: 
    :param low:
    :param high:
    :param ndigits: 
    :return: 
    """

    maximal_error = rounding_error(number, ndigits)
    minimal_number_of_steps = int(math.ceil((high - low) / maximal_error))
    return minimal_number_of_steps

# -----------------------------------------------------------------

def nbits_for_ndigits_experimental(ndigits, low, high):

    """
    This function ...
    :param ndigits: 
    :param low: 
    :param high: 
    :return: 
    """

    nsteps = minimal_nsteps(low, low, high, ndigits)
    return min_nbits_for_nintegers(nsteps)

# -----------------------------------------------------------------

def binary_gray_code(n):

    """
    This function generates the Gray code for dimension n
    :param n: 
    :return: 
    """

    def gray_code_recurse(g, n):

        k = len(g)

        if n <= 0: return
        else:

            for i in range(k-1, -1, -1):

                char = '1' + g[i]
                g.append(char)

            for i in range(k-1, -1, -1):

                g[i] = '0' + g[i]

            gray_code_recurse(g, n-1)

    g = ['0','1']
    gray_code_recurse(g, n-1)

    result = []
    for entry in g: result.append([int(character) for character in entry])
    return result

# -----------------------------------------------------------------

def quaternary_gray_code(n):

    """
    This fucntion ...
    :param n: 
    :return: 
    """

    start = ['0', '1', '2', '3']

    #ncombinations = min_nbits_for_nintegers()

    #total = []

    previous = start
    #new = None

    # Do n-1 times: e.g. for 3 qits, there's two steps from the start to the
    # http://www.eetimes.com/author.asp?section_id=14&doc_id=1283114
    for _ in range(n-1):

        previous_length = len(previous)

        # Add flipped orders
        new = previous + previous[::-1] + previous + previous[::-1]

        # Add first characters
        for i in range(previous_length):

            new[i] = '0' + new[i]

        for i in range(previous_length, 2*previous_length):

            new[i] = '1' + new[i]

        for i in range(previous_length, 3*previous_length):

            new[i] = '2' + new[i]

        for i in range(previous_length, 4*previous_length):

            new[i] = '3' + new[i]

        previous = new

    #print(previous)

    codes = previous

    result = []
    for entry in codes: result.append([int(character) for character in entry])
    return result

# -----------------------------------------------------------------

def binary_string_to_gray_binary_string(bits):

    """
    This fucntion ...
    :param bits: 
    :return: 
    """

    return bits[:1] + [i ^ ishift for i, ishift in zip(bits[:-1], bits[1:])]

# -----------------------------------------------------------------

def gray_binary_string_to_binary_string(bits):

    """
    This fucntion ...
    :param bits: 
    :return: 
    """
    
    b = [bits[0]]
    for nextb in bits[1:]: b.append(b[-1] ^ nextb)
    return b

# -----------------------------------------------------------------

def gray_binary_string_to_integer(binary_string):

    """
    This function ...
    :param binary_string: 
    :return: 
    """

    normal_binary_string = gray_binary_string_to_binary_string(binary_string)
    return binary_string_to_integer(normal_binary_string)

# -----------------------------------------------------------------

def integer_to_gray_binary_string(integer, nbits=None):

    """
    This fucntion ...
    :param integer: 
    :param nbits:
    :return: 
    """

    normal_binary_string = integer_to_binary_string(integer, nbits=nbits)
    return binary_string_to_gray_binary_string(normal_binary_string)

# -----------------------------------------------------------------

def gray_binary_string_to_float(binary_string, low, high, nbits):

    """
    This function ...
    :param binary_string:
    :param low:
    :param high:
    :param nbits:
    :return: 
    """

    binary_string = gray_binary_string_to_binary_string(binary_string)
    return binary_string_to_float(binary_string, low, high, nbits)

# -----------------------------------------------------------------

def float_to_gray_binary_string(value, low, high, nbits):

    """
    This function ...
    :param value: 
    :param low:
    :param high:
    :param nbits:
    :return: 
    """

    binary_string = float_to_binary_string(value, low, high, nbits)
    return binary_string_to_gray_binary_string(binary_string)

# -----------------------------------------------------------------

def next_binary_string(binary_string):

    """
    This function ...
    :param binary_string: 
    :return: 
    """

    nbits = len(binary_string)
    max_integer = max_integer_for_nbits(nbits)
    integer = binary_string_to_integer(binary_string)
    if integer == max_integer: raise ValueError("Cannot increment: last binary for this number of bits")
    return integer_to_binary_string(integer + 1, nbits=nbits)

# -----------------------------------------------------------------

def next_gray_binary_string(binary_string):

    """
    This function ...
    :param binary_string: 
    :return: 
    """

    nbits = len(binary_string)
    max_integer = max_integer_for_nbits(nbits)
    integer = gray_binary_string_to_integer(binary_string)
    if integer == max_integer: raise ValueError("Cannot increment: last binary for this number of bits")
    return integer_to_gray_binary_string(integer + 1, nbits=nbits)

# -----------------------------------------------------------------

def generate_bit_slices(nbits):

    """
    This function ...
    :return: 
    """

    slices = []

    tempsum = 0
    for index in range(len(nbits)):
        nbits = nbits[index]
        slices.append(slice(tempsum, tempsum+nbits))
        tempsum += nbits

    return slices

# -----------------------------------------------------------------

def random_linear(low, high):

    """
    This function ...
    :param low: 
    :param high: 
    :return: 
    """

    random = np.random.uniform(low, high)
    return random

# -----------------------------------------------------------------

def random_logarithmic(low, high):

    """
    This function ...
    :param low: 
    :param high: 
    :return: 
    """

    # Generate random logarithmic variate
    logmin = np.log10(low)
    logmax = np.log10(high)
    lograndom = np.random.uniform(logmin, logmax)
    random = 10 ** lograndom
    return random

# -----------------------------------------------------------------

def to_scale(value, scale):

    """
    This funciton ...
    :param value: 
    :param scale: 
    :return: 
    """

    if scale == "linear": return value
    elif scale == "logarithmic": return np.log10(value)
    else: raise ValueError("Invalid scale " + scale)

# -----------------------------------------------------------------

def unscale(value, scale):

    """
    This funciton ...
    :param value: 
    :param scale: 
    :return: 
    """

    if scale == "linear": return value
    elif scale == "logarithmic": return 10**value
    else: raise ValueError("Invalid scale " + scale)

# -----------------------------------------------------------------

def round_down_to_int(number):

    """
    This function ...
    :param number:
    :return:
    """

    return int(floor(number))

# -----------------------------------------------------------------

def round_down_to_integer(number):
    return round_down_to_int(number)

# -----------------------------------------------------------------

def round_down_to_even_integer(number):

    """
    This function ...
    :param number:
    :return:
    """

    integer = round_down_to_integer(number)
    if is_even(integer): return integer
    else: return integer - 1

# -----------------------------------------------------------------

def round_down_to_odd_integer(number):

    """
    This function ...
    :param number:
    :return:
    """

    integer = round_down_to_integer(number)
    if is_odd(integer): return integer
    else: return integer - 1

# -----------------------------------------------------------------

def round_up_to_int(number):

    """
    This function ...
    :param number:
    :return:
    """

    return int(ceil(number))

# -----------------------------------------------------------------

def round_up_to_integer(number):
    return round_up_to_int(number)

# -----------------------------------------------------------------

def round_up_to_even_integer(number):

    """
    This function ...
    :param number:
    :return:
    """

    integer = round_up_to_integer(number)
    if is_even(integer): return integer
    else: return integer + 1

# -----------------------------------------------------------------

def round_up_to_odd_integer(number):

    """
    This function ...
    :param number:
    :return:
    """

    integer = round_up_to_integer(number)
    if is_odd(integer): return integer
    else: return integer + 1

# -----------------------------------------------------------------

def round_to_base(number, base):

    """
    This function ...
    :param number:
    :param base:
    :return:
    """

    return int(base * round(float(number) / base))

# -----------------------------------------------------------------

def round_up_to_base(number, base):

    """
    This function ...
    :param number:
    :param base:
    :return:
    """

    return int(base * ceil(float(number) / base))

# -----------------------------------------------------------------

def round_down_to_base(number, base):

    """
    This function ...
    :param number:
    :param base:
    :return:
    """

    return int(base * floor(float(number) / base))

# -----------------------------------------------------------------

def round_to_int(number):

    """
    This function ...
    :param number:
    :return:
    """

    return int(round(number))

# -----------------------------------------------------------------

def root(number, n):

    """
    This function ...
    :param number:
    :param n:
    :return:
    """

    return number**(1./float(n))

# -----------------------------------------------------------------

def arithmetic_mean(*numbers):
    
    """
    This function ...
    :param numbers: 
    :return: 
    """

    nnumbers = len(numbers)
    return sum(numbers) / float(nnumbers)

# -----------------------------------------------------------------

def arithmetic_mean_numpy(array):

    """
    This function ...
    :param array:
    :return:
    """

    return np.mean(array)

# -----------------------------------------------------------------

def variance(*numbers, **kwargs):

    """
    This function ...
    :param numbers:
    :return:
    """

    nnumbers = len(numbers)
    mean = kwargs.pop("mean", arithmetic_mean(*numbers))
    return sum(squared_differences(numbers, mean)) / (nnumbers - 1)

# -----------------------------------------------------------------

def standard_deviation(*numbers, **kwargs):

    """
    This function ...
    :param numbers:
    :param kwargs:
    :return:
    """

    var = variance(*numbers, **kwargs)
    return math.sqrt(var)

# -----------------------------------------------------------------

def standard_deviation_numpy(array):

    """
    This function ...
    :param array:
    :return:
    """

    return np.std(array)

# -----------------------------------------------------------------

def weighed_arithmetic_mean(numbers, weights):

    """
    This function ...
    :param numbers:
    :param weights:
    :return:
    """

    #return np.sum(numbers * weights) / np.sum(weights)
    return sum([number * weight for number, weight in zip(numbers, weights)]) / sum(weights)

# -----------------------------------------------------------------

def weighed_arithmetic_mean_numpy(data, weights=None):

    """
    Calculate the weighted mean of an array/list using numpy
    """

    # Not weighted
    if weights is None: return arithmetic_mean_numpy(data)

    import numpy as np

    # Get the number of dimensions
    ndim_data = len(data.shape)
    ndim_weights = len(weights.shape)

    #weights = np.array(weights).flatten() / float(sum(weights))
    #return np.dot(np.array(data), weights)

    if ndim_weights > 1:

        weights = np.copy(weights)
        divisors = np.sum(weights, axis=-1)
        #norm_weights = weights /
        norm_weights = np.moveaxis(weights, -1, 0) # move last to first axis
        #print("1", norm_weights.shape)
        # Loop over
        for index in range(norm_weights.shape[0]): norm_weights[index] /= divisors
        #print(norm_weights.shape)
        norm_weights = np.moveaxis(norm_weights, 0, 1)
        #print("2", norm_weights.shape)

    else: norm_weights = weights / float(np.sum(weights))

    return np.dot(data, norm_weights)
    #return np.dot(norm_weights, data)

# -----------------------------------------------------------------

def weighted_median(data, weights=None):

    """
    Calculate the weighted median of a list
    """

    # Not weighted
    if weights is None: return median(data)

    midpoint = 0.5 * sum(weights)

    if any([j > midpoint for j in weights]):
        return data[weights.index(max(weights))]

    if any([j > 0 for j in weights]):
        sorted_data, sorted_weights = zip(*sorted(zip(data, weights)))
        cumulative_weight = 0
        below_midpoint_index = 0
        while cumulative_weight <= midpoint:
            below_midpoint_index += 1
            cumulative_weight += sorted_weights[below_midpoint_index-1]
        cumulative_weight -= sorted_weights[below_midpoint_index-1]
        if cumulative_weight == midpoint:
            bounds = sorted_data[below_midpoint_index-2:below_midpoint_index]
            return sum(bounds) / float(len(bounds))
        return sorted_data[below_midpoint_index-1]

# -----------------------------------------------------------------

def weighed_median_numpy(data, weights=None):

    """
    Calculate the weighted median of an array/list using numpy
    """

    import numpy as np

    # Not weighed?
    if weights is None: return np.median(np.array(data).flatten())

    data, weights = np.array(data).flatten(), np.array(weights).flatten()

    if any(weights > 0):

        sorted_data, sorted_weights = map(np.array, zip(*sorted(zip(data, weights))))
        midpoint = 0.5 * sum(sorted_weights)

        if any(weights > midpoint): return (data[weights == np.max(weights)])[0]
        cumulative_weight = np.cumsum(sorted_weights)
        below_midpoint_index = np.where(cumulative_weight <= midpoint)[0][-1]

        if cumulative_weight[below_midpoint_index] == midpoint: return np.mean(sorted_data[below_midpoint_index:below_midpoint_index + 2])
        return sorted_data[below_midpoint_index + 1]

# -----------------------------------------------------------------

def squared_differences(numbers, reference):

    """
    This function ...
    :param numbers:
    :param reference:
    :return:
    """

    squared_differences = []
    for number in numbers:
        difference = number - reference
        squared_diff = difference ** 2
        squared_differences.append(squared_diff)
    return squared_differences

# -----------------------------------------------------------------

def weighed_squared_differences(numbers, reference, weights):

    """
    This function ...
    :param numbers:
    :param reference:
    :param weights:
    :return:
    """

    result = []
    sq_diff = squared_differences(numbers, reference)
    for squared_difference, weight in zip(sq_diff, weights):
        weighed = squared_difference * weight
        result.append(weighed)
    return result

# -----------------------------------------------------------------

def weighed_standard_deviation(numbers, weights, mean=None):

    """
    This function ...
    :param numbers:
    :param weights:
    :param mean:
    :return:
    """

    if mean is None: mean = weighed_arithmetic_mean(numbers, weights)
    return np.sqrt(sum(weighed_squared_differences(numbers, mean, weights)) / sum(weights))

# -----------------------------------------------------------------

def weighed_standard_deviation_numpy(numbers, weights, mean=None):

    """
    This function ...
    :param numbers:
    :param weights:
    :param mean:
    :return:
    """

    if mean is None: mean = weighed_arithmetic_mean_numpy(numbers, weights)

    sq_diffs = (numbers - mean)**2
    norm_weights = weights / float(np.sum(weights))
    return np.sqrt(np.dot(sq_diffs, norm_weights))

# -----------------------------------------------------------------

def geometric_mean(*numbers):

    """
    This function ...
    :param numbers:
    :return:
    """

    nnumbers = len(numbers)
    product = reduce(mul, numbers, 1)
    return root(product, nnumbers)

# -----------------------------------------------------------------

def weighed_geometric_mean(numbers, weights):

    """
    This function ...
    :param numbers:
    :param weights:
    :return:
    """

    log_numbers = [np.log10(number) for number in numbers]
    return 10**(weighed_arithmetic_mean(log_numbers, weights))

# -----------------------------------------------------------------

def median(*numbers):

    """
    This function ...
    :param numbers:
    :return:
    """

    return median_numpy(np.asarray(numbers))

# -----------------------------------------------------------------

def median_numpy(array):

    """
    This function ...
    :param array:
    :return:
    """

    return np.median(array)

# -----------------------------------------------------------------

def text2int(textnum, numwords=None):

    """
    Thisf function ...
    :param textnum:
    :param numwords:
    :return:
    """

    if numwords is None:
      numwords = {}

      units = [
        "zero", "one", "two", "three", "four", "five", "six", "seven", "eight",
        "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
        "sixteen", "seventeen", "eighteen", "nineteen",
      ]

      tens = ["", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"]

      scales = ["hundred", "thousand", "million", "billion", "trillion"]

      numwords["and"] = (1, 0)

      for idx, word in enumerate(units):    numwords[word] = (1, idx)
      for idx, word in enumerate(tens):     numwords[word] = (1, idx * 10)
      for idx, word in enumerate(scales):   numwords[word] = (10 ** (idx * 3 or 2), 0)

    current = result = 0
    for word in textnum.split():
        if word not in numwords:
          raise Exception("Illegal word: " + word)

        scale, increment = numwords[word]
        current = current * scale + increment
        if scale > 100:
            result += current
            current = 0

    return result + current

# -----------------------------------------------------------------

number_list = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]
teen_list = ["ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen"]
decades_list = ["twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"]

# -----------------------------------------------------------------

def int2text(number):

    """
    This function ...
    :param number:
    :return:
    """

    if number <= 9: return number_list[number].capitalize()

    elif number >= 10 and number <= 19:

        tens = int(number % 10)
        return teen_list[tens].capitalize()

    elif number > 19 and number <= 99:

        ones = int(math.floor(number / 10))
        twos = ones - 2
        tens = int(number % 10)

        if tens == 0: return decades_list[twos].capitalize()

        else: return decades_list[twos].capitalize() + " " + number_list[tens]

    else: raise ValueError("Your numbers is not supported")

# -----------------------------------------------------------------

def nearest_integer(number):

    """
    This function ...
    :param number:
    :return:
    """

    return int(round(number))

# -----------------------------------------------------------------

def nearest_even_integer(number):

    """
    This function ...
    :param number:
    :return:
    """

    return nearest_integer(number/2)*2

# -----------------------------------------------------------------

def nearest_odd_integer(number):

    """
    This function ...
    :param number:
    :return:
    """

    if is_integer(number):
        if is_odd(number): return number
        else: return number + 1

    new = nearest_integer(number)
    if is_odd(new): return new
    elif new < number: return new + 1
    elif new > number: return new - 1
    else: raise RuntimeError("We shouldn't get here")

# -----------------------------------------------------------------

def nearest_integer_below(number, below, including=False):

    """
    This function ...
    :param number:
    :param below:
    :param including:
    :return:
    """

    integer = nearest_integer(number)

    if not is_integer(below):
        below = int(below)
        including = True

    if integer > below:

        if including: integer = below
        else: integer = below - 1

    # Return the integer
    return integer

# -----------------------------------------------------------------

def nearest_integer_above(number, above, including=False):

    """
    This function ...
    :param number:
    :param above:
    :param including:
    :return:
    """

    integer = nearest_integer(number)

    if not is_integer(above):
        above = int(ceil(above))
        including = True

    if integer < above:

        if including: integer = above
        else: integer = above + 1

    # Return the integer
    return integer

# -----------------------------------------------------------------

def nearest_even_integer_below(number, below, including=False):

    """
    This function ...
    :param number:
    :param below:
    :param including:
    :return:
    """

    integer = nearest_even_integer(number)

    # Convert float to integer if necessary
    if not is_integer(below):
        below = int(below)
        including = True

    if integer > below:

        if is_odd(below): integer = below - 1
        elif including: integer = below
        else: integer = below - 2

    # Return the integer
    return integer

# -----------------------------------------------------------------

def nearest_even_integer_above(number, above, including=False):

    """
    This function ...
    :param number:
    :param above:
    :param including:
    :return:
    """

    integer = nearest_even_integer(number)

    # Convert float to integer if necessary
    if not is_integer(above):
        above = int(ceil(above))
        including = True

    if integer < above:

        if is_odd(above): integer = above + 1
        elif including: integer = above
        else: integer = above + 2

    # Return the integer
    return integer

# -----------------------------------------------------------------

def nearest_odd_integer_below(number, below, including=False):

    """
    This function ....
    :param number:
    :param below:
    :param including:
    :return:
    """

    integer = nearest_odd_integer(number)

    # Convert float to integer if necessary
    if not is_integer(below):
        below = int(below)
        including = True

    if integer > below:

        if is_even(below): integer = below - 1
        elif including: integer = below
        else: integer = below - 2

    # Return the integer
    return integer

# -----------------------------------------------------------------

def nearest_odd_integer_above(number, above, including=False):

    """
    This function ...
    :param number:
    :param above:
    :param including:
    :return:
    """

    integer = nearest_odd_integer(number)

    # Convert float to integer if necessary
    if not is_integer(above):
        above = int(ceil(above))
        including = True

    if integer < above:

        if is_even(above): integer = above + 1
        elif including: integer = above
        else: integer = above + 2

    # Return the integer
    return integer

# -----------------------------------------------------------------

def lowest_missing_integer(integers):

    """
    This function ...
    :param integers:
    :return:
    """

    if len(integers) == 0: return 0

    integers = sorted(integers)
    lowest_missing = max(integers) + 1

    # Adjust to lower to whichever is missing
    for index in range(max(integers)):
        if integers[index] != index:
            lowest_missing = index
            break

    # Return the lowest missing integer
    return lowest_missing

# -----------------------------------------------------------------

def missing_integers(integers):

    """
    This function ...
    :param integers:
    :return:
    """

    if len(integers) == 0: return []

    missing = []
    for index in range(max(integers)):
        if index not in integers: missing.append(index)

    # Return the missing integers
    return missing

# -----------------------------------------------------------------

def sigma_clip_mask(values, sigma_level=3.0, return_nmasked=False):

    """
    This function ...
    :param values:
    :param sigma_level:
    :param return_nmasked:
    :return:
    """

    from astropy.stats import sigma_clip

    # Sigma-clip
    masked_array = sigma_clip(values, sigma=sigma_level, iters=None)

    # Get the mask
    mask = list(masked_array.mask)

    # Get the number of masked values
    nmasked = np.sum(masked_array.mask)

    # Return the mask as a list
    if return_nmasked: return mask, nmasked
    else: return mask

# -----------------------------------------------------------------

def sigma_clip(values, sigma_level=3.0, return_nmasked=False, logarithmic=False):

    """
    This function ...
    :param values:
    :param sigma_level:
    :param return_nmasked:
    :param logarithmic: the values are logarithmically distributed
    :return:
    """

    from astropy.stats import sigma_clip

    # Sigma-clip
    if logarithmic: values = np.log10(logarithmic)
    masked_array = sigma_clip(values, sigma=sigma_level, iters=None)

    # Get the number of masked values
    nmasked = np.sum(masked_array.mask)

    # Get the clipped list of values
    clipped = list(masked_array.compressed())

    # Return as list
    if return_nmasked: return clipped, nmasked
    else: return clipped

# -----------------------------------------------------------------

def fwhm(values, probabilities):

    """
    This function ...
    :param values:
    :param probabilities: frequencies/probabilities/counts
    :return:
    """

    # Check if sorted
    if not sequences.is_sorted(values): raise ValueError("Values must be sorted")

    # Get the height of the half maximum
    max_probability_index = np.argmax(probabilities)
    max_probability = probabilities[max_probability_index]
    half_max_probability = 0.5 * max_probability

    # Initialize variables for the minimum and maximum value of the FWHM span
    min_value_fwhm = None
    max_value_fwhm = None

    # Find the minimum value
    for value, probability in zip(values, probabilities):
        if probability > half_max_probability:
            min_value_fwhm = value
            break
    if min_value_fwhm is None: raise RuntimeError("Something went wrong")

    # Find the maximum value
    for value, probability in zip(reversed(values), reversed(probabilities)):
        if probability > half_max_probability:
            max_value_fwhm = value
            break
    if max_value_fwhm is None: raise RuntimeError("Something went wrong")

    # Return the FWHm
    return max_value_fwhm - min_value_fwhm

# -----------------------------------------------------------------

def same_sign(number_a, number_b):

    """
    This function ...
    :param number_a:
    :param number_b:
    :return:
    """

    return np.sign(number_a) == np.sign(number_b)

# -----------------------------------------------------------------

def different_sign(number_a, number_b):

    """
    This function ...
    :param number_a:
    :param number_b:
    :return:
    """

    return np.sign(number_a) != np.sign(number_b)

# -----------------------------------------------------------------

def differences(numbers):

    """
    This function ...
    """

    result = []
    nnumbers = len(numbers)
    for index in range(1, nnumbers):
        difference = float(numbers[index]) - float(numbers[index-1])
        result.append(difference)
    return result

# -----------------------------------------------------------------

def quotients(numbers):

    """
    This function ...
    :param numbers:
    :return:
    """

    result = []
    nnumbers = len(numbers)
    for index in range(1, nnumbers):
        quotient = float(numbers[index]) / float(numbers[index-1])
        result.append(quotient)
    return result

# -----------------------------------------------------------------

def get_linear_series(npoints, start=1, step=1):

    """
    This function ...
    :param npoints:
    :param start:
    :param step:
    :return:
    """

    stop = start + step * npoints
    return list(range(start, stop, step))

# -----------------------------------------------------------------

def get_alternating_series(npoints, start=1, step=1):

    """
    This function ...
    :param npoints:
    :param start:
    :param step:
    :return:
    """

    numbers = []
    #signs = [1, -1]
    stop = start + npoints
    for i in range(start, stop, step):
        numbers.append(i)
        if len(numbers) == npoints: return numbers
        numbers.append(-i)
        if len(numbers) == npoints: return numbers

    #return numbers
    raise RuntimeError("We shouldn't get here")

# -----------------------------------------------------------------

def closest_half_integer(number):

    """
    Round a number to the closest half integer.
    >>> closest_half_integer(1.3)
    1.5
    >>> closest_half_integer(2.6)
    2.5
    >>> closest_half_integer(3.0)
    3.0
    >>> closest_half_integer(4.1)
    4.0
    """

    return round(number * 2) / 2.0

# -----------------------------------------------------------------

def nearest_half_integer(number):
    return closest_half_integer(number)

# -----------------------------------------------------------------
