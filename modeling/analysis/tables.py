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

    @property
    def cache_host_ids(self):

        """
        This function ...
        :return:
        """

        return list(set(self["Host id"]))

    # -----------------------------------------------------------------

    def run_names_for_host_id(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        names = []
        for index in range(len(self)):
            host_id_index = self["Host id"][index]
            if host_id_index == host_id: names.append(self["Run name"][index])
        return names

    # -----------------------------------------------------------------

    def index_for_run_name(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.run_names.index(run_name)

    # -----------------------------------------------------------------

    def host_id_for_run_name(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        index = self.index_for_run_name(run_name)
        return self["Host id"][index]

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
