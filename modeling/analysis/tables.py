#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.tables Contains different table classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable

# -----------------------------------------------------------------

class CachedRunsTable(SmartTable):

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
        super(CachedRunsTable, self).__init__(*args, **kwargs)

        # Add columns
        self.add_column_info("Run name", str, None, "Name of the analysis run")
        self.add_column_info("Host id", str, None, "Cache remote host ID")

    # -----------------------------------------------------------------

    @property
    def run_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Run name"])

    # -----------------------------------------------------------------

    def add_entry(self, run_name, host_id):

        """
        This function ...
        :param run_name:
        :param host_id:
        :return:
        """

        values = [run_name, host_id]
        self.add_row(values)

# -----------------------------------------------------------------
