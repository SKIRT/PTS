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
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        columns = np.loadtxt(path, unpack=True)

        number_of_columns = len(columns)

        data = []
        names = []
        units = dict()

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

        table = tables.new(data, names)

        for column_name in units:
            table[column_name].unit = units[column_name]

        return table

# -----------------------------------------------------------------
