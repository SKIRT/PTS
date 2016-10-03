#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.history Contains the ModelingHistory class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.tools import time

# -----------------------------------------------------------------

class ModelingHistory(SmartTable):
    
    """
    This class...
    """

    column_info = [("Command", str, None, "name of the modeling command"),
                   ("Start time", str, None, "timestamp for start of command"),
                   ("End time", str, None, "timestamp for end of command")]

    # -----------------------------------------------------------------

    def add_entry(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Set the values
        values = [command, time.timestamp(), None]

        # Add a row to the table
        self.add_row(values)

    # -----------------------------------------------------------------

    def mark_end(self):

        """
        This function ...
        :return:
        """

        timestamp = time.timestamp()

        self._resize_string_column("End time", timestamp)

        # Set the value
        self["End time"].mask[-1] = False
        self["End time"][-1] = timestamp

# -----------------------------------------------------------------
