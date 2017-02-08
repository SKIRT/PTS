#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.arrays Provides useful functions for dealing with arrays and derived classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import numpy as np

# Import the relevant PTS classes and modules
from ..basics.unit import parse_unit

# -----------------------------------------------------------------

def find_closest_index(array, value, array_unit=None):

    """
    This function ...
    :param array:
    :param value:
    :param array_unit:
    :return:
    """

    closest_delta = None
    #closest_delta = float("inf")
    closest_index = None

    column_unit = array_unit
    value_unit = value.unit if hasattr(value, "unit") else None

    # Check units
    if value_unit is not None:
        if column_unit is None: raise ValueError("Value has a unit but column has not: cannot compare these values")
        else: value = value.to(column_unit).value # for correct comparison inside loop
    elif column_unit is not None: raise ValueError("Value has no unit but the column has: cannot compare these values")

    # Loop over all entries in the column
    for i in range(len(array)):

        delta = abs(array[i] - value)

        if closest_delta is None or delta < closest_delta:
            closest_delta = delta
            closest_index = i

    # Return index of the closest
    return closest_index

# -----------------------------------------------------------------

def find_closest_above_index(array, value, array_unit=None):

    """
    This function ...
    :param array:
    :param value:
    :param array_unit:
    :return:
    """

    column_unit = array_unit
    value_unit = value.unit if hasattr(value, "unit") else None

    # Check units
    if value_unit is not None:
        if column_unit is None: raise ValueError("Value has a unit but column has not: cannot compare these values")
        else: value = value.to(column_unit).value  # for correct comparison inside loop
    elif column_unit is not None: raise ValueError("Value has no unit but the column has: cannot compare these values")

    # Loop from beginning, return index of first value that is greater
    for i in range(len(array)):

        table_value = array[i]
        if table_value > value: return i

    # Return None
    return None

# -----------------------------------------------------------------

def find_closest_below_index(array, value, array_unit=None):

    """
    This function ...
    :param array:
    :param value:
    :param array_unit:
    :return:
    """

    column_unit = array_unit
    value_unit = value.unit if hasattr(value, "unit") else None

    # Check units
    if value_unit is not None:
        if column_unit is None: raise ValueError("Value has a unit but column has not: cannot compare these values")
        else: value = value.to(column_unit).value  # for correct comparison inside loop
    elif column_unit is not None: raise ValueError("Value has no unit but the column has: cannot compare these values")

    # Loop from end (reversed loop), return index of first value that is smaller
    for i in reversed(range(len(array))):

        table_value = array[i]
        if table_value < value: return i

    # Return None
    return None

# -----------------------------------------------------------------

def array_as_list(array, add_unit=True, unit=None, masked_value=None, array_unit=None, conversion_info=None):

    """
    This function ...
    :param array:
    :param add_unit:
    :param unit:
    :param masked_value:
    :param array_unit:
    :param conversion_info:
    :return:
    """

    if conversion_info is None: conversion_info = dict()

    # Initialize a list to contain the column values
    result = []

    #if array_unit is None:
    #    if hasattr(array, "unit"): array_unit = array.unit
    has_unit = array.unit is not None
    has_mask = hasattr(array, "mask")

    # If unit has to be converted, check whether the original unit is specified
    if not has_unit and unit is not None: raise ValueError("Cannot determine the unit of the column so values cannot be converted to " + str(unit))

    # Loop over the entries in the column
    for i in range(len(array)):

        if has_mask and array.mask[i]: result.append(masked_value)
        else:

            if has_unit:

                # Add the unit initially to be able to convert
                value = array[i] * array_unit

                # Set the wavelength
                conversion_info = copy.deepcopy(conversion_info)
                if "wavelengths" in conversion_info:
                    conversion_info["wavelength"] = conversion_info["wavelengths"][i]
                    del conversion_info["wavelengths"]

                # If a target unit is specified, convert
                if unit is not None: value = value.to(unit, **conversion_info).value * parse_unit(unit) # If converted, do not add any unit

                if not add_unit: value = value.value # If not converted and add_unit is enabled, add the unit

            else: value = array[i]

            # Add the value to the list
            result.append(value)

    # Return the list
    return result

# -----------------------------------------------------------------

def plain_array(column, unit=None, array_unit=None, conversion_info=None):

    """
    This function ...
    :param column:
    :param unit:
    :param array_unit:
    :param conversion_info:
    :return:
    """

    return np.array(array_as_list(column, unit=unit, add_unit=False, masked_value=float('nan'), array_unit=array_unit, conversion_info=conversion_info))

# -----------------------------------------------------------------
