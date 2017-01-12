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
from string import ascii_lowercase
from .lists import interleave

# -----------------------------------------------------------------

alphabet = list(ascii_lowercase)

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
