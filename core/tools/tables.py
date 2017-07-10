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

# Import the relevant PTS classes and modules
from . import arrays

# -----------------------------------------------------------------

def find_index(table, key, column_name=None, where=None):

    """
    This function ...
    :param table:
    :param key:
    :param column_name:
    :param where:
    :return:
    """

    # If multiple keys are passed
    if isinstance(key, list):

        if column_name is None: raise ValueError("Column names must be specified when specifying multiple keys")
        if not isinstance(column_name, list): raise ValueError("If key(s) is a list, column_name(s) must also be a list")

        # Loop over all entries in the table
        for i in range(len(table)):

            if where is not None:
                if skip_entry_based_on_where(table, i, where): continue

            found_mismatch = False

            for k, c in zip(key, column_name):

                if not (table[c][i] == k):
                    found_mismatch = True
                    break

            if not found_mismatch: return i

        return None

    # Just a single key is passed
    else:

        # Get first column name if none is given
        if column_name is None: column_name = table.colnames[0]

        # Loop over all entries in the column
        for i in range(len(table)):

            # Skip
            if where is not None:
                if skip_entry_based_on_where(table, i, where): continue

            if table[column_name][i] == key: return i

        return None

# -----------------------------------------------------------------

def skip_entry_based_on_where(table, index, where):

    """
    This function ...
    :param table: 
    :param index: 
    :param where: 
    :return: 
    """

    skip_entry = False
    for cname in where:
        #print(table[cname][index], where[cname])
        if table[cname][index] != where[cname]:
            skip_entry = True
            break

    return skip_entry

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

def equal_columns(columns):

    """
    This function ...
    :param columns:
    :return:
    """

    first_column = columns[0]
    for column in columns[1:]:
        if not np.array_equal(column, first_column): return False
    return True

# -----------------------------------------------------------------

def write(table, path, format="ascii.ecsv"):

    """
    This function ...
    :param table
    :param path:
    :param format:
    :return:
    """

    # Write the table
    table.write(path, format=format)

# -----------------------------------------------------------------

def from_file(path, format="ascii.ecsv", fix_floats=False, fix_string_length=False):

    """
    This function ...
    :param path:
    :param format:
    :param fix_floats:
    :param fix_string_length:
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

def new(data=None, names=None, meta=None, dtypes=None, copy=True):

    """
    This function ...
    :param data:
    :param names:
    :param meta:
    :param dtypes:
    :param copy:
    :return:
    """

    if data is None: return Table()
    else:

        # Create a new table from the data
        table = Table(data=data, names=names, meta=meta, masked=True, dtype=dtypes, copy=copy)

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

def columns_as_objects(columns, cls, add_unit=True, unit=None, column_units=None, density=False, brightness=False):

    """
    This function ...
    :param columns:
    :param cls:
    :param add_unit:
    :param unit:
    :param column_units:
    :param density:
    :param brightness:
    :return:
    """

    # Initialize a list to contain the objects
    result = []

    # Number of columns
    ncols = len(columns)

    # Initialize a list to contain the column values
    column_lists = []

    if column_units is None: column_units = [None] * ncols

    # Loop over the columns, get list of values
    for column, column_unit in zip(columns, column_units):
        column_lists.append(arrays.array_as_list(column, add_unit=add_unit, unit=unit, array_unit=column_unit, density=density, brightness=brightness))

    # Loop over the column values for each entry
    for i in range(len(column_lists[0])):

        arguments = []
        for j in range(ncols): arguments.append(column_lists[j][i])

        if all(argument is None for argument in arguments): obj = None
        else:

            try: obj = cls(*arguments)
            except Exception: obj = None

        result.append(obj)

    # Return the resulting list of objects
    return result

# -----------------------------------------------------------------

def column_names(table):

    """
    This function ...
    :param table:
    :return:
    """

    return table.colnames

# -----------------------------------------------------------------

def as_tuples(table):

    """
    This function ...
    :param table:
    :return:
    """

    tuples = []
    for i in range(len(table)):
        tuples.append(tuple(list(table[i])))
    return tuples

# -----------------------------------------------------------------
