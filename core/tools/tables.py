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

def from_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Read the table from file
    fill_values = [('--', '0')]
    table = Table.read(path, fill_values=fill_values, format="ascii.commented_header")

    # Fix boolean values
    fix_logical(table)

    # Return the new table
    return table

# -----------------------------------------------------------------

def new(data, names, meta=None):

    """
    This function ...
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
