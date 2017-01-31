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
from collections import OrderedDict

# Import astronomical modules
from astropy.table import Table, MaskedColumn

# Import the relevant PTS classes and modules
from .unit import PhotometricUnit
from .unit import parse_unit as u

# -----------------------------------------------------------------

class SmartTable(Table):

    """
    This class ...
    """

    # default extension
    default_extension = "dat"

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Always used masked tables
        kwargs["masked"] = True

        # Call the constructor of the base class
        super(SmartTable, self).__init__(*args, **kwargs)

        # Column info
        self.column_info = []

        # Path
        self.path = None

    # -----------------------------------------------------------------

    def add_column_info(self, name, dtype, unit, description):

        """
        This function ...
        :param name:
        :param dtype:
        :param unit:
        :param description:
        :return:
        """

        self.column_info.append((name, dtype, unit, description))

    # -----------------------------------------------------------------

    def column_unit(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        if self[column_name].unit is None: return None

        # Construct unit
        if column_name in self.meta["density"]: density = True
        else: density = False
        return u(self[column_name].unit, density=density)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Setup has already been called
        if len(self.colnames) != 0: return

        # Create the table names and types lists
        for entry in self.column_info:

            name = entry[0]
            dtype = entry[1]
            unit = entry[2]

            data = []

            # Add column
            col = MaskedColumn(data=data, name=name, dtype=dtype, unit=unit)
            self.add_column(col)

            # Set whether this column is a spectral density
            if isinstance(unit, PhotometricUnit) and unit.density:
                if "density" not in self.meta: self.meta["density"] = []
                self.meta["density"].append(name)

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

        # Clear the column info so that we can rebuild it
        table.column_info = []

        # Set the column info
        # Loop over the columns
        for name in table.colnames:

            # Get the type
            dtype = table[name].dtype
            if np.issubdtype(dtype, np.string_): simple_dtype = str
            elif np.issubdtype(dtype, np.float): simple_dtype = float
            elif np.issubdtype(dtype, np.int): simple_dtype = int
            elif np.issubdtype(dtype, np.bool): simple_dtype = bool
            else: raise ValueError("Did not recognize the dtype of column '" + name + "'")

            # Get unit of the column
            unit = table.column_unit(name)

            # Add column info
            table.add_column_info(name, simple_dtype, unit, None)

            # Initialize "density" meta
            if "density" not in table.meta: table.meta["density"] = []

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

                # Quantity with photometric unit
                if isinstance(value.unit, PhotometricUnit):

                    factor = value.unit.conversion_factor(column_unit)
                    scalar_value = value.value * factor

                # Quantity with regular Astropy Unit
                else: scalar_value = value.to(column_unit).value

                # Add the value
                scalar_values.append(scalar_value)

            # A scalar value (or string, int, ...)
            else: scalar_values.append(value)

        # Return the values without unit
        return scalar_values

    # -----------------------------------------------------------------

    def _convert_lists(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        converted_values = []

        for i, value in enumerate(values):

            if isinstance(value, list):

                column_type = self.column_info[i][1]

                if len(value) == 0: converted_value = None
                elif column_type == str: converted_value = ",".join(map(str, value))
                else: raise ValueError("Cannot have a list element in the row at the column that is not of string type")

            else: converted_value = value

            converted_values.append(converted_value)

        # Return the converted values
        return converted_values

    # -----------------------------------------------------------------

    def get_quantity(self, colname, index):

        """
        This function ...
        :param colname:
        :param index:
        :return:
        """

        value = self[colname][index]

        if self[colname].mask[index]: return None
        elif self[colname].unit is not None:
            quantity = value * self[colname].unit
        else: quantity = value

        # Return the quantity
        return quantity

    # -----------------------------------------------------------------

    def get_row(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        row = OrderedDict()

        for name in self.colnames:

            value = self[name][index]

            if self[name].mask[index]: value = None
            elif self[name].unit is not None:

                # Add unit
                value = value * self[name].unit

            # Add the value
            row[name] = value

        # Return the row
        return row

    # -----------------------------------------------------------------

    def add_row(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        # Setup if necessary
        if len(self.colnames) == 0: self.setup()

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Strip units
        values = self._strip_units(values)

        # Convert lists to string
        values = self._convert_lists(values)

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

    def all_equal(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        if len(self[column_name]) == 0: return True

        # Doesn't work with strings I think ...
        #return (self[column_name] == self[column_name][0]).all() # all equal to the first element

        first = self[column_name][0]
        for i in range(len(self[column_name])):
            if self[column_name][i] != first: return False
        return True

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

        # Setup if necessary
        if len(self.colnames) == 0: self.setup()

        # Write the table in ECSV format
        self.write(path, format="ascii.ecsv")

        # Set the path
        self.path = path

# -----------------------------------------------------------------
