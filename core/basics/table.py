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

# Import standard modules
import numpy as np

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

        fill_values = [('--', '0')]

        # Open the table
        table = super(SmartTable, cls).read(path, format="ascii.ecsv", fill_values=fill_values)

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

            if values[index] is None: new_string_length = 0
            else: new_string_length = len(values[index])

            if new_string_length > current_string_length: new_sizes[colname] = new_string_length

        # Resize columns
        for colname in new_sizes:

            current_resized = self[colname].astype("S" + str(new_sizes[colname]))
            self.replace_column(colname, current_resized)

    # -----------------------------------------------------------------

    def _resize_string_column(self, colname, value):

        """
        This function ...
        :param colname:
        :param value:
        :return:
        """

        if value is None: return

        dtype_str = str(self[colname].dtype)

        if not dtype_str.startswith("|S"): raise ValueError("Column " + colname + " is not a column of strings")

        current_string_length = int(dtype_str.split("S")[1])

        new_string_length = len(value)

        if new_string_length > current_string_length:

            new_size = new_string_length

            # Replace the column by a resized one
            current_resized = self[colname].astype("S" + str(new_size))
            self.replace_column(colname, current_resized)

    # -----------------------------------------------------------------

    def _strip_units(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        scalar_values = []

        for i, value in enumerate(values):

            # If this value has a unit, we have to make sure it is converted into the proper column unit
            if hasattr(value, "unit"):

                column_unit = self.column_info[i][2]
                assert column_unit is not None

                scalar_value = value.to(column_unit).value

                scalar_values.append(scalar_value)

            # A scalar value (or string, int, ...)
            else: scalar_values.append(value)

        # Return the values without unit
        return scalar_values

    # -----------------------------------------------------------------

    def add_row(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Strip units
        values = self._strip_units(values)

        # Create mask
        mask = [value is None for value in values]

        # Set masked values to have a default value (None will not work for Astropy)
        new_values = []
        for i in range(len(values)):

            if values[i] is not None: new_values.append(values[i])
            else:
                coltype = self[self.colnames[i]].dtype.name
                if coltype.startswith("string"): new_values.append("")
                elif coltype.startswith("float"): new_values.append(0.)
                elif coltype.startswith("int"): new_values.append(0)
                else: raise ValueError("Unknown column type: " + coltype)

        # Add the row
        super(SmartTable, self).add_row(new_values, mask=mask)

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
