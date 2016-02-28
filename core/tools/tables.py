#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.tables Provides useful functions for dealing with tables

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.table import Table, Column

# -----------------------------------------------------------------

def find_index(table, key):

    """
    This function ...
    :param key:
    :return:
    """

    first_column_name = table.colnames[0]

    # Loop over all entries in the column
    for i in range(len(table)):

        if table[first_column_name][i] == key: return i

    return None

# -----------------------------------------------------------------

def write(table, path):

    """
    This function ...
    :param catalog:
    :param path:
    :return:
    """

    # TODO: add metadata ?

    # Write the table
    table.write(path, format="ascii.commented_header")

# -----------------------------------------------------------------

def from_file(path, format="ascii.commented_header"):

    """
    This function ...
    :param path:
    :param format:
    :return:
    """

    # Read the table from file
    fill_values = [('--', '0')]
    table = Table.read(path, fill_values=fill_values, format=format)

    # Fix boolean values
    fix_logical(table)

    # Return the new table
    return table

# -----------------------------------------------------------------

def new(data, names, meta=None):

    """
    This function ...
    :param names:
    :param meta:
    :return:
    """

    # Create a new table from the data
    table = Table(data, names=names, meta=meta, masked=True)

    # Set mask for each column from None values
    for column_index in range(len(names)):
        table[names[column_index]].mask = [value is None for value in data[column_index]]

    # Return the new table
    return table

# -----------------------------------------------------------------

def fix_logical(table):
    
    """
    This function ...
    """

    for column in table.columns.values():
        if column.dtype.str.endswith('S5'):
            falses = column == 'False'
            trues = column == 'True'
            if np.all(falses | trues):

                bool_column = Column(trues)
                table.replace_column(column.name, bool_column)

# -----------------------------------------------------------------

def column_as_list(column, add_unit=True, unit=None):

    """
    This function ...
    :param column:
    :param add_unit:
    :param unit:
    :return:
    """

    # Initialize a list to contain the column values
    result = []

    has_unit = column.unit is not None
    has_mask = hasattr(column, "mask")

    # If unit has to be converted, check whether the original unit is specified
    if not has_unit and unit is not None: raise ValueError("Cannot determine the unit of the column so values cannot be converted to " + str(unit))

    # Loop over the entries in the column
    for i in range(len(column)):

        if has_mask and column.mask[i]: result.append(None)
        else:

            if has_unit:

                # Add the unit initially to be able to convert
                value = column[i] * column.unit

                # If a target unit is specified, convert
                if unit is not None: value = value.to(unit).value # If converted, do not add any unit
                elif not add_unit: value = column[i] # If not converted and add_unit is enabled, add the unit

            else: value = column[i]

            # Add the value to the list
            result.append(value)

    # Return the list
    return result

# -----------------------------------------------------------------

def columns_as_objects(columns, cls, add_unit=True, unit=None):

    """
    This function ...
    :param columns:
    :param cls:
    :param add_unit:
    :param unit:
    :return:
    """

    # Initialize a list to contain the objects
    result = []

    # Number of columns
    ncols = len(columns)

    # Initialize a list to contain the column values
    column_lists = []

    # Loop over the columns, get list of values
    for column in columns: column_lists.append(column_as_list(column, add_unit=add_unit, unit=unit))

    # Loop over the column values for each entry
    for i in range(len(column_lists[0])):

        arguments = []
        for j in range(ncols): arguments.append(column_lists[j][i])

        if all(argument is None for argument in arguments): obj = None
        else: obj = cls(*arguments)

        result.append(obj)

    # Return the resulting list of objects
    return result

# -----------------------------------------------------------------
