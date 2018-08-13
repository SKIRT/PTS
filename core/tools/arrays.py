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
from ..units.parsing import parse_unit
from ..units.quantity import PhotometricQuantity

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

def find_closest_above_index(array, value, array_unit=None, inclusive=False):

    """
    This function ...
    :param array:
    :param value:
    :param array_unit:
    :param inclusive:
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

        if inclusive:
            if table_value >= value: return i
        else:
            if table_value > value: return i

    # Return None
    return None

# -----------------------------------------------------------------

def find_closest_below_index(array, value, array_unit=None, inclusive=False):

    """
    This function ...
    :param array:
    :param value:
    :param array_unit:
    :param inclusive:
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

        if inclusive:
            if table_value <= value: return i
        else:
            if table_value < value: return i

    # Return None
    return None

# -----------------------------------------------------------------

def array_as_list(array, add_unit=True, unit=None, masked_value=None, array_unit=None, conversion_info=None,
                  density=False, brightness=False, equivalencies=None, mask=None):

    """
    This function ...
    :param array:
    :param add_unit:
    :param unit:
    :param masked_value:
    :param array_unit:
    :param conversion_info:
    :param density:
    :param brightness:
    :param equivalencies:
    :param mask:
    :return:
    """

    # If conversion info is not defined
    wavelengths = None
    if conversion_info is None:
        conversion_info = dict()

    # Wavelengths are specified in the conversion info
    elif "wavelengths" in conversion_info:
        wavelengths = conversion_info["wavelengths"]
        conversion_info = copy.deepcopy(conversion_info)
        del conversion_info["wavelengths"]

    #print(wavelengths)

    # Initialize a list to contain the column values
    result = []

    # Get array unit if not passed to this function
    if array_unit is None and hasattr(array, "unit"): array_unit = array.unit

    # Parse the unit
    if unit is not None: unit = parse_unit(unit, density=density, brightness=brightness)

    #print("array unit", array_unit)
    #print("unit", unit)

    has_unit = array_unit is not None
    has_mask = hasattr(array, "mask")

    # If unit has to be converted, check whether the original unit is specified
    if not has_unit and unit is not None: raise ValueError("Cannot determine the unit of the column so values cannot be converted to " + str(unit))

    # Loop over the entries in the column
    for i in range(len(array)):

        # Masked value: skip
        if mask is not None and mask[i]: continue

        # Masked array value
        if has_mask and array.mask[i]: result.append(masked_value)

        # Not masked value
        else:

            # If the array has a unit
            if has_unit:

                # Add the unit initially to be able to convert
                value = array[i] * array_unit

                # Check if a target unit is specified
                if unit is not None:

                    # Check if the requested unit is different from the array unit
                    if unit != array_unit:

                        # Set the wavelength
                        if isinstance(value, PhotometricQuantity):
                            if wavelengths is not None: conversion_info["wavelength"] = wavelengths[i]
                        # Not a photometric quantity
                        else: conversion_info = dict()

                        # Add equivalencies
                        if equivalencies is not None: conversion_info["equivalencies"] = equivalencies

                        #print(value, value.physical_type)
                        #print(unit, unit.physical_type)
                        #print(conversion_info)

                        # If a target unit is specified, convert
                        value = value.to(unit, **conversion_info).value * unit  # If converted, do not add any unit

                    # Don't add the unit: get the value
                    if not add_unit: value = value.value # If not converted and add_unit is enabled, add the unit

                # If adding the unit is not requested but no unit was specified for the output, then the user wouldn't have any information on what the units are of the resulting list, so throw an error
                elif not add_unit: raise ValueError("You cannot know which units the values are going to be if you don't specifiy the target unit and you put add_unit to False")

            # The array doesn't have a unit
            else:

                if unit is not None: raise ValueError("Cannot convert the unit because the original unit is unknown")
                value = array[i]

            # Add the value to the list
            result.append(value)

    # Return the list
    return result

# -----------------------------------------------------------------

def plain_array(column, unit=None, array_unit=None, conversion_info=None, density=False, brightness=False,
                equivalencies=None, mask=None):

    """
    This function ...
    :param column:
    :param unit:
    :param array_unit:
    :param conversion_info:
    :param density:
    :param brightness:
    :param equivalencies:
    :param mask:
    :return:
    """

    # Array has unit and target unit is specified?
    if unit is not None and array_unit is not None:

        # Parse the target unit
        unit = parse_unit(unit, density=density, brightness=brightness)

        # Same unit as column?
        if unit == array_unit: return _column_to_array(column, mask=mask)

        #elif TODO: CAN THERE BE ONE CONVERSION FACTOR FOR THE ENTIRE ARRAY?

        else: return _plain_array_from_column(column, unit=unit, array_unit=array_unit, conversion_info=conversion_info,
                                              equivalencies=equivalencies, mask=mask)

    # Array has no unit and target unit is also None
    elif unit is None and array_unit is None: return _column_to_array(column, mask=mask)

    # All other cases
    else: return _plain_array_from_column(column, unit=unit, array_unit=array_unit, conversion_info=conversion_info,
                                          equivalencies=equivalencies, mask=mask)

# -----------------------------------------------------------------

def _column_to_array(column, mask=None):

    """
    This function ...
    :param column:
    :param mask:
    :return:
    """

    array = np.asarray(column)
    if mask is not None: array = array[np.logical_not(mask)]
    return array

# -----------------------------------------------------------------

def _plain_array_from_column(column, unit=None, array_unit=None, conversion_info=None, density=False, brightness=False,
                            equivalencies=None, mask=None):

    """
    Thisnf unction ...
    :param column:
    :param unit:
    :param array_unit:
    :param conversion_info:
    :param density:
    :param brightness:
    :param equivalencies:
    :param mask:
    :return:
    """

    return np.array(array_as_list(column, unit=unit, add_unit=False, masked_value=float('nan'), array_unit=array_unit,
                                  conversion_info=conversion_info, density=density, brightness=brightness,
                                  equivalencies=equivalencies, mask=mask))

# -----------------------------------------------------------------
