#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.generations Contains the GenerationsTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable

# -----------------------------------------------------------------

class GenerationsTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Generation name", str, None, "Name for the generation"),
                   ("Generation index", int, None, "Index of the generation"),
                   ("Wavelength grid level", int, None, "level of the wavelength gid"),
                   ("Dust grid level", int, None, "level of the dust grid"),
                   ("Number of simulations", int, None, "number of simulations (individuals) in the generation"),
                   ("Self-absorption", bool, None, "dust self-absorption enabled")]

    @classmethod
    def initialize(cls, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        for label in parameters:

            cls.column_info.append(("Minimum value for " + label, float, None, "minimum value for " + label))
            cls.column_info.append(("Maximum value for " + label, float, None, "Maximum value for " + label))

        # Call the initialize function of the generations table function
        return super(GenerationsTable, cls).initialize()

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        # Open the table
        table = super(GenerationsTable, cls).read(path, format="ascii.ecsv")

        # Set the path
        table.path = path

        # Return the table
        return table

    # -----------------------------------------------------------------

    @property
    def parameter_names(self):

        """
        This function ...
        :return:
        """

        names = []

        for name in self.colnames:

            if name.startswith("Minimum value for "):
                names.append(name.split("value for ")[1])

        return names

    # -----------------------------------------------------------------

    def add_entry(self, name, index, wavelength_grid_level, dust_grid_level, nsimulations, selfabsorption, ranges):

        """
        This function ...
        :param name:
        :param index:
        :param wavelength_grid_level:
        :param dust_grid_level:
        :param nsimulations:
        :param selfabsorption:
        :param ranges:
        :return:
        """

        values = [name, index, wavelength_grid_level, dust_grid_level, nsimulations, selfabsorption]

        # Add the boundaries (min, max) of the parameter ranges as seperate column values
        for name in self.parameter_names:
            values.append(ranges[name].min)
            values.append(ranges[name].max)

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
