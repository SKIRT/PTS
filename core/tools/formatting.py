#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.formatting Formatting text in the terminal.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .stringify import stringify_list_fancy, stringify, tostr
from . import filesystem as fs
from .logging import log
from .sequences import equal_sizes

# -----------------------------------------------------------------

# SOURCE: http://misc.flogisoft.com/bash/tip_colors_and_formatting

# -----------------------------------------------------------------

# Set
bold = "\033[1m"
dim = "\033[2m"
underlined = "\033[4m"
blink = "\033[5m"
inverted = "\033[7m"
hidden = "\033[8m"

# -----------------------------------------------------------------

# Reset
reset = "\033[0m"
reset_bold = "\033[21m"
reset_dim = "\033[22m"
reset_underlined = "\033[24m"
reset_blink = "\033[25m"
reset_inverted = "\033[27m"
reset_hidden = "\033[28m"

# -----------------------------------------------------------------

# Text colours
default_text = "\033[39m"
black = "\033[30m"
red = "\033[31m"
green = "\033[32m"
yellow = "\033[33m"
blue = "\033[34m"
magenta = "\033[35m"
cyan = "\033[36m"
lightgray = "\033[37m"
darkgray = "\033[90m"
lightred = "\033[91m"
lightgreen = "\033[92m"
lightyellow = "\033[93m"
lightblue = "\033[94m"
lightmagenta = "\033[95m"
lightcyan = "\033[96m"
white = "\033[97m"

# -----------------------------------------------------------------

# Background colours
default_background = "\033[49m"
black_background = "\033[40m"
red_background = "\033[41m"
green_background = "\033[42m"
yellow_background = "\033[43m"
blue_background = "\033[44m"
magenta_background = "\033[45m"
cyan_background = "\033[46m"
lightgray_background = "\033[47m"
darkgray_background = "\033[100m"
lightred_background = "\033[101m"
lightgreen_background = "\033[102m"
lightyellow_background = "\033[103m"
lightblue_background = "\033[104m"
lightmagenta_background = "\033[105m"
lightcyan_background = "\033[106m"
white_background = "\033[107m"

# -----------------------------------------------------------------

def center_text_around(text, ch, length=50):

    """
    This function ...
    :param text:
    :param ch:
    :param length:
    :return:
    """

    spaced_text = '%s' % text
    result = spaced_text.center(length, ch)
    return result

# -----------------------------------------------------------------

def print_empty():

    """
    This function ...
    :return:
    """

    print("")

# -----------------------------------------------------------------

def print_filled(ch, length=50, prefix=""):

    """
    This function ...
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    result = ch * length
    print(prefix + result)

# -----------------------------------------------------------------

def print_border(ch, length=50, prefix=""):

    """
    This function ...
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    text = " " * (length - 2)
    result = center_text_around(text, ch, length)
    print(prefix + result)

# -----------------------------------------------------------------

def print_centered_around(text, ch, length=50, prefix=""):

    """
    This function ...
    :param text:
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    result = center_text_around(text, ch, length)
    print(prefix + result)

# -----------------------------------------------------------------

def print_centered_around_border(text, ch, length=50, prefix=""):

    """
    This function ...
    :param text:
    :param ch:
    :param length:
    :param prefix:
    :return:
    """

    without_border = center_text_around(text, " ", length=length-2)
    result = center_text_around(without_border, ch, length=length)
    print(prefix + result)

# -----------------------------------------------------------------

def print_files_in_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    if not log.is_debug(): return
    in_path = fs.files_in_path(path, returns="name", extensions=True)
    if len(in_path) == 0: log.debug("No files in '" + path + "'")
    else:
        log.debug(str(len(in_path)) + " files in '" + path + "':")
        print("")
        print(stringify_list_fancy(in_path, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_files_in_list(lst, name, only_name=False):

    """
    This function ...
    :param lst:
    :param name:
    :param only_name:
    :return:
    """

    if not log.is_debug(): return
    if len(lst) == 0: log.debug("No files in '" + name + "'")
    else:
        log.debug(str(len(lst)) + " files in '" + name + "':")
        print("")
        if only_name: strings = [fs.name(path) for path in lst]
        else: strings = lst
        print(stringify_list_fancy(strings, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_directories_in_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    if not log.is_debug(): return
    in_path = fs.directories_in_path(path, returns="name")
    if len(in_path) == 0: log.debug("No directories in '" + path + "'")
    else:
        log.debug(str(len(in_path)) + " directories in '" + path + "':")
        print("")
        print(stringify_list_fancy(in_path, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_directories_in_list(lst, name):

    """
    This function ...
    :param lst:
    :param name:
    :return:
    """

    if not log.is_debug(): return
    if len(lst) == 0: log.debug("No directories in '" + name + "'")
    else:
        log.debug(str(len(lst)) + " directories in '" + name +"':")
        print("")
        print(stringify_list_fancy(lst, lines_prefix="  ")[1])
        print("")

# -----------------------------------------------------------------

def print_dictionary(dictionary):

    """
    This function ...
    :param dictionary: 
    :return: 
    """

    print("")
    for label in dictionary: print(" - " + label + ": " + stringify(dictionary[label])[1])
    print("")

# -----------------------------------------------------------------

def print_columns(*columns, **kwargs):

    """
    This function ...
    :param columns:
    :return:
    """

    delimiter = kwargs.pop("delimiter", "  ")

    # Check sizes
    if not equal_sizes(*columns): raise ValueError("Columns must have equal lengths")

    # Convert all to strings
    string_columns = []

    max_lengths = []

    for column in columns:

        new_column = [tostr(entry) for entry in column]
        lengths = [len(string) for string in new_column]

        max_lengths.append(max(lengths))
        string_columns.append(new_column)

    ncolumns = len(string_columns)

    for i in range(len(columns[0])):

        row = ""

        for j in range(ncolumns): # number of columns

            part = string_columns[j][i]
            row += part

            if j != ncolumns - 1:

                nextra = max_lengths[j] - len(part)
                spaces = " " * nextra + delimiter

                row += spaces

        print(row)

# -----------------------------------------------------------------
