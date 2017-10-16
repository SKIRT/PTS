#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.table Contains the SkirtTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..tools import tables

# -----------------------------------------------------------------

class TruncatedSKIRTTableError(Exception):

    """
    This class ...
    """

    def __init__(self, message, path=None):

        """
        Thisf unction ...
        :param message:
        :param path:
        """

        # Call the base class constructor with the parameters it needs
        super(TruncatedSKIRTTableError, self).__init__(message)

        # The FITS file path
        self.path = path

# -----------------------------------------------------------------

def is_valid(path):

    """
    This function ...
    :param path:
    :return:
    """

    columns = np.loadtxt(path, unpack=True)

    number_of_columns = len(columns)

    # Try to interpret the number of rows
    try: number_of_rows = len(columns[0])
    except TypeError:  # object of type 'numpy.float64' has no len()
        #raise TruncatedSKIRTTableError("The file only contains one line", path=path)
        return False

    from ..tools import filesystem as fs

    # Get first ncolumn+1 rows of the file = the header
    header = fs.get_first_lines(path, number_of_columns + 1)

    # Loop over the lines
    has_title = False
    for index, line in enumerate(header):

        # Set flags
        first = index == 0
        last = index == number_of_columns

        if last and not has_title: break # one row less

        # Should start as comment
        if not line.startswith("#"): return False

        # Should state column index, no not always!
        if not line.startswith("# column "):
            if not first: return False
            else: has_title = True

    # Every check passed
    return True

# -----------------------------------------------------------------

class SkirtTable(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        pass

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, expected_nrows=None):

        """
        This function ...
        :param path:
        :param expected_nrows:
        :return:
        """

        columns = np.loadtxt(path, unpack=True)

        number_of_columns = len(columns)

        try: number_of_rows = len(columns[0])
        except TypeError: # object of type 'numpy.float64' has no len()
            raise TruncatedSKIRTTableError("The file only contains one line", path=path)

        # Check expected number of rows
        if expected_nrows is not None and number_of_rows != expected_nrows:
            raise IOError("Expected " + str(expected_nrows) + " rows but only found " + str(number_of_rows))

        data = []
        names = []
        units = dict()

        # Get names and units
        with open(path) as table_file:

            for i in range(number_of_columns):

                line = table_file.next()
                column_number = i + 1

                if "column " + str(column_number) not in line:

                    if i != 0: raise IOError("Column name and unit of column " + str(i+1) + " not found in file header")
                    else:

                        # Determine offset
                        while "column" not in line:
                            line = table_file.next()

                name_and_unit = line.split(": ")[1].split("\n")[0]
                if "(" in name_and_unit and ")" in name_and_unit:
                    name = name_and_unit.split(" (")[0].capitalize()
                    unit = name_and_unit.split(" (")[1].split(")")[0]
                else:
                    name = name_and_unit.capitalize()
                    unit = None

                data.append(columns[i])
                names.append(name)

                if unit is not None: units[name] = unit

        # Construct the table
        table = tables.new(data, names)

        # Set the column units
        for column_name in units:
            table[column_name].unit = units[column_name]

        # Return the table
        return table

# -----------------------------------------------------------------
