#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.table Contains the SmartTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.table import Table

# -----------------------------------------------------------------

class SmartTable(Table):

    """
    This class ...
    """

    # column_info defined in the sublasses

    @classmethod
    def initialize(cls):

        """
        This function ...
        :return:
        """

        # Create the table names and types lists
        names = []
        dtypes = []
        for entry in cls.column_info:
            name = entry[0]
            dtype = entry[1]
            names.append(name)
            dtypes.append(dtype)

        # Call the constructor of the base class
        table = cls(names=names, dtype=dtypes, masked=True)

        # Set the column units
        for entry in cls.column_info:

            name = entry[0]
            unit = entry[2]

            if unit is None: continue

            # Set the column unit
            table[name].unit = unit

        # Add the path attribute
        table.path = None

        # Return the smart table instance
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        # Open the table
        table = super(SmartTable, cls).read(path, format="ascii.ecsv")

        # Set the path
        table.path = path

        # Return the table
        return table

    # -----------------------------------------------------------------

    def _resize_string_columns(self, values):

        """
        This function ...
        :return:
        """

        new_sizes = dict()

        colnames = self.colnames
        for index in range(len(colnames)):

            colname = colnames[index]

            dtype_str = str(self[colname].dtype)

            if not dtype_str.startswith("|S"): continue

            current_string_length = int(dtype_str.split("S")[1])

            new_string_length = len(values[index])

            if new_string_length > current_string_length: new_sizes[colname] = new_string_length

        # Resize columns
        for colname in new_sizes:

            current_resized = self[colname].astype("S" + str(new_sizes[colname]))
            self.replace_column(colname, current_resized)

    # -----------------------------------------------------------------

    def add_row(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        mask = [value is None for value in values]
        values = [self.column_info[i][1]("0") if values[i] is None else values[i] for i in range(len(values))]
        super(SmartTable, self).add_row(values, mask=mask)

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        if self.path is None: raise RuntimeError("Path has not been set yet")

        # Save to the current path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Write the table in ECSV format
        self.write(path, format="ascii.ecsv")

        # Set the path
        self.path = path

# -----------------------------------------------------------------
