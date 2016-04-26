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
from astropy.units import Unit

# -----------------------------------------------------------------

def find_index(table, key, column_name=None):

    """
    This function ...
    :param table:
    :param key:
    :param column_name:
    :return:
    """

    if isinstance(key, list):

        if column_name is None: raise ValueError("Column names must be specified when specifying multiple keys")
        if not isinstance(column_name, list): raise ValueError("If key(s) is a list, column_name(s) must also be a list")

        # Loop over all entries in the table
        for i in range(len(table)):

            found_mismatch = False

            for k, c in zip(key, column_name):

                if not (table[c][i] == k):
                    found_mismatch = True
                    break

            if not found_mismatch: return i

        return None

    #elif isinstance(key, basestring):
    else:

        # Get first column name if none is given
        if column_name is None: column_name = table.colnames[0]

        # Loop over all entries in the column
        for i in range(len(table)):

            if table[column_name][i] == key: return i

        return None

    #else: raise ValueError("Invalid key: must be a list (of strings) or a string")

# -----------------------------------------------------------------

def find_indices(table, key, column_name=None):

    """
    This function ...
    :param table:
    :param key:
    :param column_name:
    :return:
    """

    if isinstance(key, list):

        if column_name is None: raise ValueError("Column names must be specified when specifying multiple keys")
        if not isinstance(column_name, list): raise ValueError("If key(s) is a list, column_name(s) must also be a list")

        indices = []

        # Loop over all entries in the table
        for i in range(len(table)):

            found_mismatch = False

            for k, c in zip(key, column_name):

                if not (table[c][i] == k):
                    found_mismatch = True
                    break

            if not found_mismatch: indices.append(i)

        return indices

    #elif isinstance(key, basestring):
    else:

        # Get first column name is none is given
        if column_name is None: column_name = table.colnames[0]

        indices = []

        # Loop over all entries in the column
        for i in range(len(table)):

            if table[column_name][i] == key: indices.append(i)

        return indices

    #else: raise ValueError("Invalid key: must be a list (of strings) or a string")

# -----------------------------------------------------------------

def write(table, path, format="ascii.commented_header"):

    """
    This function ...
    :param table
    :param path:
    :param format:
    :return:
    """

    # TODO: add metadata ?

    # Write the table
    table.write(path, format=format)

# -----------------------------------------------------------------

def from_file(path, format="ascii.commented_header", fix_floats=False, fix_string_length=False):

    """
    This function ...
    :param path:
    :param format:
    :param fix_float:
    :return:
    """

    # Read the table from file
    fill_values = [('--', '0')]
    table = Table.read(path, fill_values=fill_values, format=format)

    # Fix boolean values
    fix_logical(table)

    if fix_floats: fix_float(table)
    # Sometimes, a column of floats is parsed as a column of strings ... But then importing this function from the python command line and loading the same table does work ... straaange..

    if fix_string_length: fix_string_length_column(table, fix_string_length[0], fix_string_length[1])

    # Return the new table
    return table

# -----------------------------------------------------------------

def new(data, names, meta=None, dtypes=None):

    """
    This function ...
    :param data:
    :param names:
    :param meta:
    :param dtypes:
    :return:
    """

    # Create a new table from the data
    table = Table(data, names=names, meta=meta, masked=True, dtype=dtypes)

    # Set mask for each column from None values
    for column_index in range(len(names)):
        table[names[column_index]].mask = [value is None for value in data[column_index]]

    # Return the new table
    return table

# -----------------------------------------------------------------

def fix_logical(table):
    
    """
    This function ...
    :param table:
    """

    for column in table.columns.values():
        if column.dtype.str.endswith('S5'):
            falses = column == 'False'
            trues = column == 'True'
            if np.all(falses | trues):

                bool_column = Column(trues)
                table.replace_column(column.name, bool_column)

# -----------------------------------------------------------------

def fix_float(table):

    """
    This function ...
    :param table:
    :return:
    """

    for column in table.columns.values():

        if str(column.dtype).startswith("|S"):

            failed = False

            values = []

            for value in column:

                try:
                    value = float(value)
                    values.append(value)
                except ValueError:
                    failed = True
                    break # break the loop over the values in the column

            if not failed:

                float_column = Column(values)
                table.replace_column(column.name, float_column)

# -----------------------------------------------------------------

def fix_string_length_column(table, column_name, length):

    """
    This function ...
    :param table:
    :param column_name:
    :param length:
    :return:
    """

    length_str = "S"+str(length)
    table[column_name].dtype = length_str

# -----------------------------------------------------------------

def column_as_list(column, add_unit=True, unit=None, masked_value=None):

    """
    This function ...
    :param column:
    :param add_unit:
    :param unit:
    :param masked_value:
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

        if has_mask and column.mask[i]: result.append(masked_value)
        else:

            if has_unit:

                # Add the unit initially to be able to convert
                value = column[i] * column.unit

                # If a target unit is specified, convert
                if unit is not None: value = value.to(unit).value * Unit(unit) # If converted, do not add any unit

                if not add_unit: value = value.value # If not converted and add_unit is enabled, add the unit

            else: value = column[i]

            # Add the value to the list
            result.append(value)

    # Return the list
    return result

# -----------------------------------------------------------------

def column_as_array(column, unit=None):

    """
    This function ...
    :param column:
    :param unit:
    :return:
    """

    return np.array(column_as_list(column, unit=unit, add_unit=False, masked_value=float('nan')))

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
