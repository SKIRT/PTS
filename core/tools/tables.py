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
from collections import defaultdict

# Import astronomical modules
from astropy.table import Table, Column

# Import the relevant PTS classes and modules
from . import arrays
from . import types
from . import numbers

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
    if types.is_sequence(key):

        # Checks
        if column_name is None: raise ValueError("Column names must be specified when specifying multiple keys")
        if not types.is_sequence(column_name): raise ValueError("If key(s) is a list, column_name(s) must also be a list")

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

    # Sequence
    if types.is_sequence(key):

        # Checks
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

    # Single key
    else:

        # Get first column name is none is given
        if column_name is None: column_name = table.colnames[0]

        indices = []

        # Loop over all entries in the column
        for i in range(len(table)):

            if table[column_name][i] == key: indices.append(i)

        return indices

# -----------------------------------------------------------------

def find_closest_indices(table, key, column_name=None, nindices=5, rel_or_abs="rel", quadratic=True):

    """
    This function ...
    :param table:
    :param key:
    :param column_name:
    :param nindices:
    :param rel_or_abs:
    :param quadratic:
    :return:
    """

    # Keep the differences
    diffs = []

    # If multiple keys are passed
    if types.is_sequence(key):

        # Checks
        if column_name is None: raise ValueError("Column names must be specified when specifying multiple keys")
        if not types.is_sequence(column_name): raise ValueError("If key(s) is a list, column_name(s) must also be a list")

        # Set rel or abs dict
        if types.is_dictionary(rel_or_abs): relabs = rel_or_abs
        elif types.is_sequence(rel_or_abs) or types.is_tuple(rel_or_abs):
            relabs = dict()
            for k, ra in zip(key, rel_or_abs): relabs[k] = ra
        elif types.is_string_type(rel_or_abs):
            relabs = dict()
            for k in key: relabs[k] = rel_or_abs
        else: raise ValueError("Invalid value for 'rel_or_abs'")

        # Loop over all entries in the table
        for i in range(len(table)):

            diff_keys = []

            # Loop over the columns and reference value
            for k, c in zip(key, column_name):

                value = table[c][i]

                # Calculate the difference
                if relabs[k] == "abs": diff = abs(value - k)
                elif relabs[k] == "rel":
                    if numbers.different_sign(value, k): diff = float("inf")
                    else: diff = np.exp(abs(np.log(value/k)))
                else: raise ValueError("Invalid option: '" + relabs[k] + "'")

                # Add the difference
                diff_keys.append(diff)

            # Calculate difference
            if quadratic: diff = np.sqrt(np.sum(np.power(diff_keys, 2)))
            else: diff = sum(diff_keys)

            # Add difference for this entry
            diffs.append(diff)

    # Single key
    else:

        # Get first column name is none is given
        if column_name is None: column_name = table.colnames[0]

        # Loop over all entries in the table
        for i in range(len(table)):

            # Get the value
            value = table[column_name][i]

            # Calculate the difference
            if rel_or_abs == "abs": diff = abs(value - key)
            elif rel_or_abs == "rel":
                if numbers.different_sign(value, key): diff = float("inf")
                else: diff = np.exp(abs(np.log(value/key)))
            else: raise ValueError("Invalid option: '" + rel_or_abs + "'")

            # Add the difference
            diffs.append(diff)

    # Return the indices
    indices = np.argsort(diffs)
    return indices[:nindices]

# -----------------------------------------------------------------

def find_one_index(table, key, column_name=None):

    """
    This function ...
    :param table:
    :param key:
    :param column_name:
    :return:
    """

    indices = find_indices(table, key, column_name=column_name)
    if len(indices) == 0: raise ValueError("Not found")
    elif len(indices) > 1: raise ValueError("Multiple matches found")
    else: return indices[0]

# -----------------------------------------------------------------

def find_closest_index(table, key, column_name=None, rel_or_abs="rel"):

    """
    This function ...
    :param table:
    :param key:
    :param column_name:
    :param rel_or_abs:
    :return:
    """

    indices = find_closest_indices(table, key, column_name=column_name, rel_or_abs=rel_or_abs)
    if len(indices) == 0: raise ValueError("Not found")
    else: return indices[0] # return closest

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

def from_remote_file(path, remote, format="ascii.ecsv", fix_floats=False, fix_string_length=False):

    """
    This function ...
    :param path:
    :param remote:
    :param format:
    :param fix_floats:
    :param fix_string_length:
    :return:
    """

    contents = remote.get_text(path)

    # Read the table from file
    fill_values = [('--', '0')]
    table = Table.read(contents, fill_values=fill_values, format=format)

    # Fix boolean values
    fix_logical(table)

    if fix_floats: fix_float(table)
    # Sometimes, a column of floats is parsed as a column of strings ... But then importing this function from the python command line and loading the same table does work ... straaange..

    if fix_string_length: fix_string_length_column(table, fix_string_length[0], fix_string_length[1])

    # Return the new table
    return table

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

    #if data is None: return Table()
    #else:

    # Create a new table from the data
    table = Table(data=data, names=names, meta=meta, masked=True, dtype=dtypes, copy=copy)

    if data is not None:
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

def columns_as_objects(columns, cls, add_unit=True, unit=None, column_units=None, density=False, brightness=False, conversion_info=None):

    """
    This function ...
    :param columns:
    :param cls:
    :param add_unit:
    :param unit:
    :param column_units:
    :param density:
    :param brightness:
    :param conversion_info:
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
        column_lists.append(arrays.array_as_list(column, add_unit=add_unit, unit=unit, array_unit=column_unit, density=density, brightness=brightness, conversion_info=conversion_info))

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

def filtered(table, column_name, value):

    """
    This function ...
    :param table:
    :param key:
    :param value:
    :return:
    """

    names = table.colnames

    data = []

    indices = find_indices(table, value, column_name=column_name)

    # Fill the columns
    for name in names:
        column = [value for index, value in enumerate(table[name]) if index in indices]
        data.append(column)

    # Create the filtered table
    newtable = new(data, names)
    return newtable

# -----------------------------------------------------------------
