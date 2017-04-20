#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.table import SmartTable
from pts.core.tools import tables
from pts.core.tools import sequences

# -----------------------------------------------------------------

class ScoresTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ScoresTable, self).__init__(*args, **kwargs)

        # Add column info
        self.column_info.append(("Individual name", str, None, "name of the individual"))
        self.column_info.append(("Score", float, None, "individual's score"))

        # Initialize 'min_or_max' meta attribute
        self.meta["min_or_max"] = kwargs.pop("min_or_max")

    # -----------------------------------------------------------------

    @property
    def min_or_max(self):

        """
        This function ...
        :return: 
        """

        return self.meta["min_or_max"]

    # -----------------------------------------------------------------

    @property
    def minimization(self):

        """
        This function ...
        :return: 
        """

        return self.min_or_max == "min"

    # -----------------------------------------------------------------

    @property
    def maximization(self):

        """
        This function ...
        :return: 
        """

        return self.min_or_max == "max"

    # -----------------------------------------------------------------

    @property
    def best_individual_name(self):

        """
        This function ...
        :return:
        """

        # Get index
        if self.maximization: index = np.argmax(self["Score"])
        else: index = np.argmin(self["Score"])

        # Return individual name
        return self["Individual name"][index]

    # -----------------------------------------------------------------

    @property
    def best_score(self):

        """
        This function ...
        :return: 
        """

        if self.maximization: return np.max(self["Score"])
        else: return np.min(self["Score"])

    # -----------------------------------------------------------------

    @property
    def individual_names(self):

        """
        This function ...
        :return: 
        """

        return list(self["Individual name"])

    # -----------------------------------------------------------------

    @property
    def scores(self):

        """
        This function ...
        :return: 
        """

        return list(self["Score"])

    # -----------------------------------------------------------------

    def score_for(self, individual_name):

        """
        This function ...
        :param individual_name: 
        :return: 
        """

        index = tables.find_index(self, individual_name, "Individual name")
        return self["Score"][index]

    # -----------------------------------------------------------------

    def sort_as(self, individual_names):

        """
        This function ...
        :param individual_names: 
        :return: 
        """

        if len(individual_names) != len(self): raise ValueError("Number of individuals does not have the same length as the table")
        if sequences.contains_duplicates(individual_names): raise ValueError("Invalid list of individual names: multiple occurences of the same name")

        new_column = []

        # Loop over the individual names
        for name in individual_names: new_column.append(self.score_for(name))

        # Replace the column
        self["Score"] = new_column

    # -----------------------------------------------------------------

    def add_entry(self, name, score):

        """
        This function ...
        :param name:
        :param score:
        :return:
        """

        # Set the values
        values = [name, score]

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
