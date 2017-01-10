#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.curve Contains the Curve class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .table import SmartTable

# -----------------------------------------------------------------

class Curve(SmartTable):

    """
    This class ...
    """

    column_info = []

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, x_unit=None, y_unit=None, x_name="x", y_name="y", x_description="x values", y_description="y values"):

        """
        This function ...
        :param x_unit:
        :param y_unit:
        :param x_name:
        :param y_name:
        :param x_description:
        :param y_description:
        :return:
        """

        # Add columns
        cls.column_info.append((x_name, float, str(x_unit), x_description))
        cls.column_info.append((y_name, float, str(y_unit), y_description))

        # Call the initialize function of the SmartTable table function
        return super(Curve, cls).initialize()

    # -----------------------------------------------------------------

    def add_entry(self, x_value, y_value):

        """
        This function ...
        :param x_value:
        :param y_value:
        :return:
        """

        # Set values
        values = [x_value, y_value]

        # Add a row
        self.add_row(values)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Sort the table by the x values
        self.sort(self.colnames[0])

        # Call the saveto function of the base class
        super(Curve, self).saveto(path)

# -----------------------------------------------------------------
