#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.strings Provides functions for dealing with strings.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
from textwrap import wrap
from string import ascii_lowercase

# Import the relevant PTS classes and modules
from .sequences import interleave
from . import types

# -----------------------------------------------------------------

alphabet = list(ascii_lowercase)

# -----------------------------------------------------------------

def lowercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.lower()

# -----------------------------------------------------------------

def uppercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.upper()

# -----------------------------------------------------------------

def capitalize(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.lower().title()

# -----------------------------------------------------------------

def is_lowercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.islower()

# -----------------------------------------------------------------

def is_uppercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.isupper()

# -----------------------------------------------------------------

def is_capitalized(string):

    """
    This function ...
    :param string:
    :return:
    """

    return capitalize(string) == string

# -----------------------------------------------------------------

def case_combinations(string, also_one_letter=True):

    """
    This function ...
    :param string:
    :param also_one_letter:
    :return:
    """

    if len(string) == 1 and not also_one_letter: return [string]

    result = []

    if is_lowercase(string):

        result.append(string)
        result.append(uppercase(string))
        result.append(capitalize(string))

    elif is_uppercase(string):

        result.append(string)
        result.append(lowercase(string))
        result.append(capitalize(string))

    elif is_capitalized(string):

        result.append(string)
        result.append(lowercase(string))
        result.append(uppercase(string))

    else:

        result.append(string)
        result.append(lowercase(string))
        result.append(uppercase(string))
        result.append(capitalize(string))

    # Return the resulting list
    return result

# -----------------------------------------------------------------

def case_combinations_list(string_list, also_one_letter=True):

    """
    This function ...
    :param string_list:
    :param also_one_letter:
    :return:
    """

    result = []
    for string in string_list: result += case_combinations(string, also_one_letter)
    return result

# -----------------------------------------------------------------

def quantity_combinations(quantity):

    """
    This function ...
    :param quantity:
    :return:
    """

    result = []

    for value_string in number_representations(quantity.value):

        result.append(value_string + " " + str(quantity.unit))
        result.append(value_string + str(quantity.unit))
        result.append(value_string)

        if quantity.unit == "micron":

            result.append(value_string + "mu")
            result.append(value_string + "um")
            result.append(value_string + " mu")
            result.append(value_string + " um")

    return result

# -----------------------------------------------------------------

def iterate_alphabet():

    """
    This function ...
    :return:
    """

    for letter in alphabet: yield letter

# -----------------------------------------------------------------

def split_except_within_single_quotes(text):

    """
    This function strips the whitespace from a string, except when it is within quotes
    :param text:
    :return:
    """

    parts = []

    lst = text.split("'")

    for i, item in enumerate(lst):

        if i % 2: parts.append("'" + item + "'")
        else:

            for a in item.split(): parts.append(a)

    return parts

# -----------------------------------------------------------------

def stripwhite_around(text, around):

    """
    This function ...
    :param text:
    :param around:
    :return:
    """

    return text.replace(" " + around, around).replace(around + " ", around)

# -----------------------------------------------------------------

def stripwhite_except_quotes(text):

    """
    This function strips the whitespace from a string, except when it is within quotes
    :param text:
    :return:
    """

    lst = text.split('"')
    for i, item in enumerate(lst):
        if not i % 2:
            lst[i] = re.sub("\s+", "", item)
    return '"'.join(lst)

# -----------------------------------------------------------------

def stripwhite_except_curlybrackets(text):

    """
    This function strips the whitespace from a string, except when it is within curly brackets
    :param text:
    :return:
    """

    if "{" in text and "}" in text:
        lst = []
        for part in text.split("{"): lst += part.split("}")
        for i, item in enumerate(lst):
            if not i % 2:
                lst[i] = re.sub("\s+", "", item)
        return "".join(list(interleave([lst, ["{", "}"]])))
    else: return text.replace(" ", "")

# -----------------------------------------------------------------

def stripwhite_except_singlequotes(text):

    """
    This function strips the whitespace from a string, except when it is within single quotes
    :param text:
    :return:
    """

    lst = text.split("'")
    for i, item in enumerate(lst):
        if not i % 2:
            lst[i] = re.sub("\s+", "", item)
    return "'".join(lst)

# -----------------------------------------------------------------

def num_to_ith(num):

    """
    1 becomes 1st, 2 becomes 2nd, etc.
    """

    value             = str(num)
    last_digit        = value[-1]
    if len(value) > 1 and value[-2] == '1': return value +'th'
    if last_digit == '1': return value + 'st'
    if last_digit == '2': return value + 'nd'
    if last_digit == '3': return value + 'rd'
    return value + 'th'

# -----------------------------------------------------------------

def number_representations(value):

    """
    This function ...
    :param value:
    :return:
    """

    strings = set()

    strings.add(str(value))
    strings.add(repr(value))
    strings.add(str_from_real_or_integer(value))

    return list(strings)

# -----------------------------------------------------------------

def str_from_real_or_integer(value):

    """
    This function ...
    :param value:
    :return:
    """

    if types.is_integer_type(value): return str(value)
    elif types.is_real_type(value):

        if int(value) == value: return str(int(value))
        else: return repr(value)

    else: raise ValueError("Value must be real or integer")

# -----------------------------------------------------------------

def generate_from_two_parts(part_a, part_b, connectors=(" ", "-", ".", "_"), also_reverse=False):

    """
    This function ...
    :param part_a:
    :param part_b:
    :param connectors:
    :param also_reverse:
    :return:
    """

    for connector in connectors: yield part_a + connector + part_b

    if also_reverse:
        for connector in connectors: yield part_b + connector + part_a

# -----------------------------------------------------------------

def find_first_digit(string):

    """
    This function ...
    :param string:
    :return:
    """

    index = re.search("\d", string)
    return index

# -----------------------------------------------------------------

def find_last_digit(string):

    """
    This function ...
    :param string:
    :return:
    """

    index = re.search("\d", string[::-1])
    return len(string) - 1 - index

# -----------------------------------------------------------------

def split_in_lines(string, length=60, as_list=False):

    """
    This function ...
    :param length:
    :param as_list:
    :return:
    """

    # Split into list
    lst = wrap(string, length)

    if as_list: return lst
    else: return "\n".join(lst)

# -----------------------------------------------------------------

def startswith_any(line, patterns):

    """
    This function ...
    :param patterns:
    :return:
    """

    for pattern in patterns:
        if line.startswith(pattern): return True
    return False

# -----------------------------------------------------------------
