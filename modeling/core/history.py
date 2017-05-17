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
from ...core.tools import time, tables

# -----------------------------------------------------------------

class ModelingHistory(SmartTable):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelingHistory, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Command", str, None, "name of the modeling command")
        self.add_column_info("Start time", str, None, "timestamp for start of command")
        self.add_column_info("End time", str, None, "timestamp for end of command")

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

    def remove_entry(self, command):

        """
        This function removes
        :param command: 
        :return: 
        """

        # Find the index of the corresponding row
        index = tables.find_index(self, command)
        self.remove_row(index)

    # -----------------------------------------------------------------

    def __contains__(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        return command_name in list(self["Command"])

    # -----------------------------------------------------------------

    def finished(self, command_name):

        """
        This function ...
        :param command_name:
        :return:
        """

        if command_name not in self: return False
        else:
            index = tables.find_index(command_name, self)
            return not self["End time"][index] # not masked

    # -----------------------------------------------------------------

    @property
    def has_configured_fit(self):

        """
        This function ...
        :return:
        """

        return "configure_fit" in self

    # -----------------------------------------------------------------

    @property
    def has_initialized_fit(self):

        """
        This function ...
        :return:
        """

        if "initialize_fit_sed" in self: return True
        elif "initialize_fit_galaxy" in self: return True
        else: return False

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
