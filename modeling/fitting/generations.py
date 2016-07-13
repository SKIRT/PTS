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

# Import astronomical modules
from astropy.table import Table

# -----------------------------------------------------------------

class GenerationsTable(Table):

    """
    This class ...
    """

    @classmethod
    def initialize(cls, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        # Create the table
        names = ["Generation name", "Generation index", "Wavelength grid level", "Dust grid level",
                 "Number of simulations", "Self-absorption"]
        dtypes = [str, int, int, int, int, bool]

        for label in parameters:
            names.append("Minimum value for " + label)
            names.append("Maximum value for " + label)
            dtypes.append(float)
            dtypes.append(float)

        # Call the constructor of the base class
        table = cls(names=names, masked=True)

        # Set the column units
        #table["a"].unit = "a_unit"
        #table["b"].unit = "b_unit"

        # Set the path
        table.path = None

        # Return the generations table instance
        return table

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

    def add_entry(self, name, index, wavelength_grid_level, dust_grid_level, nsimuations, selfabsorption, ranges):

        """
        This function ...
        :param name:
        :param index:
        :param wavelength_grid_level:
        :param dust_grid_level:
        :param nsimuations:
        :param selfabsorption:
        :param ranges:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

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
