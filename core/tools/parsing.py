#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.parsing Provides useful functions for parsing strings into a variety of types.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Unit

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from ..basics.range import IntegerRange, RealRange, QuantityRange
from ...magic.basics.vector import Vector

# -----------------------------------------------------------------

def boolean(entry):

    """
    This function ...
    :param entry:
    :return:
    """

    return entry == "True"

# -----------------------------------------------------------------

def integer(argument):

    """
    This function ....
    :param argument:
    :return:
    """

    return int(argument)

# -----------------------------------------------------------------

def real(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return float(argument)

# -----------------------------------------------------------------

def string(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return argument

# -----------------------------------------------------------------

def real_range(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    min_value, max_value = real_tuple(argument)
    return RealRange(min_value, max_value)

# -----------------------------------------------------------------

def integer_range(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    min_value, max_value = int_tuple(argument)
    return IntegerRange(min_value, max_value)

# -----------------------------------------------------------------

def quantity_range(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    min_quantity, max_quantity = quantity_tuple(argument)
    return QuantityRange(min_quantity, max_quantity)

# -----------------------------------------------------------------

def absolute_path(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return os.path.abspath(argument)

# -----------------------------------------------------------------

def int_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(int, argument.split(','))
        return a, b
    except: raise argparse.ArgumentTypeError("Tuple must be of format a,b")

# -----------------------------------------------------------------

def real_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(float, argument.split(","))
        return a, b
    except: raise argparse.ArgumentTypeError("Tuple must be of format a,b")

# -----------------------------------------------------------------

def quantity_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(quantity, argument.split(","))
        return a, b
    except: raise argparse.ArgumentTypeError("Tuple must be of format a unit_a, b unit_b")

# -----------------------------------------------------------------

def quantity_vector(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    tuple_string = argument.split("(")[1].split(")")[0]
    x, y = tuple_string.split(", ")

    x = quantity(x[2:])
    y = quantity(y[2:])

    return Vector(x, y)

# -----------------------------------------------------------------

def string_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return argument.split(",")

# -----------------------------------------------------------------

def duration(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    # Calculate the walltime in seconds
    hours, minutes, seconds = argument.split(':')
    duration = int(hours)*3600 + int(minutes)*60 + int(seconds)

    # Return the duration in seconds
    return duration

# -----------------------------------------------------------------

def int_list(string, name="ids"):

    """
    This function returns a list of integer values, based on a string denoting a certain range (e.g. '3-9') or a
    set of integer values seperated by commas ('2,14,20')
    :param string:
    :param name:
    :return:
    """

    if "-" in string and "," in string:

        parts = string.split(",")
        total_int_list = []
        for part in parts: total_int_list += int_list(part)
        return total_int_list

    # Split the string
    splitted = string.split('-')

    if len(splitted) == 0: raise argparse.ArgumentError(name, "No range given")
    elif len(splitted) == 1:

        splitted = splitted[0].split(",")

        # Check if the values are valid
        for value in splitted:
            if not value.isdigit(): raise argparse.ArgumentError(name, "Argument contains unvalid characters")

        # Only leave unique values
        return list(set([int(value) for value in splitted]))

    elif len(splitted) == 2:

        if not (splitted[0].isdigit() and splitted[1].isdigit()): raise argparse.ArgumentError(name, "Not a valid integer range")
        return range(int(splitted[0]), int(splitted[1])+1)

    else: raise argparse.ArgumentError(name, "Values must be seperated by commas or by a '-' in the case of a range")

# -----------------------------------------------------------------

def simulation_ids(string):

    """
    This function ...
    :param string:
    :return:
    """

    # Initialize a dictionary
    delete = dict()

    # If the string is empty, raise an error
    if not string.strip(): raise argparse.ArgumentError("ids", "No input for argument")

    # Split the string by the ';' character, so that each part represents a different remote host
    for entry in string.split(";"):

        # Split again to get the host ID
        splitted = entry.split(":")

        # Get the host ID
        host_id = splitted[0]

        # Get the simulation ID's
        values = int_list(splitted[1])

        # Add the simulation ID's to the dictionary for the correspoding host ID
        delete[host_id] = values

    # Return the dictionary with ID's of simulations that should be deleted
    return delete

# -----------------------------------------------------------------

def quantity(argument):

    """
    >>> quantity("2GB")
    (2.0, 'GB')
    >>> quantity("17 ft")
    (17.0, 'ft')
    >>> quantity("   3.4e-27 frobnitzem ")
    (3.4e-27, 'frobnitzem')
    >>> quantity("9001")
    (9001.0, '')
    >>> quantity("spam sandwhiches")
    (1.0, 'spam sandwhiches')
    >>> quantity("")
    (1.0, '')
    """

    # NEW IMPLEMENTATION
    units = ""
    number = 1.0
    while argument:
        try:
            number = float(argument)
            break
        except ValueError:
            units = argument[-1:] + units
            argument = argument[:-1]
    return number, Unit(units.strip())

    # FIRST IMPLEMENTATION
    #splitted = argument.split()
    #value = float(splitted[0])
    #unit = Unit(splitted[1])
    #return value * unit

    # http://stackoverflow.com/questions/2240303/separate-number-from-unit-in-a-string-in-python

    # SECOND IMPLEMENTATION
    #numeric = '0123456789-.'
    #for i, c in enumerate(argument + " "):
    #    if c not in numeric:
    #        break
    #value = argument[:i]
    #unit = Unit(argument[i:].lstrip())
    #return value * unit

# -----------------------------------------------------------------

def get_quantity(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create a quantity object and return it
    if unit is not None: value = value * Unit(unit)
    return value

# -----------------------------------------------------------------

def get_angle(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create an Angle object and return it
    if unit is not None: value = Angle(value, unit)
    return value

# -----------------------------------------------------------------
