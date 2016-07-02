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

# Import the relevant PTS classes and modules
from ...core.tools import tables
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class GenerationsTable(Table):
    
    """
    This class ...
    """

    def __init__(self, path):

        """
        This function ...
        :param path:
        """

        # Set the path of the generations table
        self.path = path

        # If the file does not exist yet, create it
        if not fs.is_file(self.path): self.initialize()

    # -----------------------------------------------------------------

    def initialize(self):

        """
        This function ...
        :return:
        """

        # Create the table
        names = ["Generation name", "Generation index", "Wavelength grid level", "Dust grid level",
                 "Number of simulations", "FUV young range"]
        data = [[] for _ in names]
        dtypes = ["S24", "S23", "S15", "S15", "int64", "int64", "int64", "int64", "int64", "int64", "bool", "bool",
                  "bool", "float64", "float64", "float64", "float64", "float64", "float64", "float64", "float64",
                  "float64"]
        table = tables.new(data, names, dtypes=dtypes)

        # Set the column units
        table["Total runtime"] = "s"  # runtimes are expressed in seconds
        table["Setup time"] = "s"
        table["Stellar emission time"] = "s"
        table["Spectra calculation time"] = "s"
        table["Dust emission time"] = "s"
        table["Writing time"] = "s"
        table["Waiting time"] = "s"
        table["Communication time"] = "s"
        table["Intermediate time"] = "s"

        # Write the (empty) table
        tables.write(table, self.path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    @classmethod
    def read(cls, path):

        """
        This function ...
        :return:
        """

        # Load the generations table
        table = tables.from_file(path, format="ascii.ecsv")

        # Return the table
        return table

    # -----------------------------------------------------------------

    def add_entry(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Open the timing file in 'append' mode
        timing_file = open(self.path, 'a')

        # Initialize a list to contain the values of the row
        row = []

# -----------------------------------------------------------------
