#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.context Contains the AnalysisContext class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.utils import lazyproperty
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .run import AnalysisRuns
from .run import AnalysisRunInfo, AnalysisRun
from .tables import CachedRunsTable

# -----------------------------------------------------------------

timing_filename = "timing.dat"
memory_filename = "memory.dat"
cached_filename = "cached.dat"

# -----------------------------------------------------------------

class AnalysisContext(object):

    """
    This class ...
    """

    def __init__(self, analysis_path):

        """
        The constructor ...
        :param analysis_path:
        """

        # Set the analysis path
        self.path = analysis_path

        # Timing table --

        # Set the path to the timing table
        self.timing_table_path = fs.join(self.analysis_path, timing_filename)

        # Initialize the timing table if necessary
        if not fs.is_file(self.timing_table_path):

            # Create the table and save it
            timing_table = TimingTable()
            timing_table.saveto(self.timing_table_path)

        # Memory table --

        # Set the path to the memory table
        self.memory_table_path = fs.join(self.analysis_path, memory_filename)

        # Initialize the memory table if necessary
        if not fs.is_file(self.memory_table_path):

            # Create the table and save it
            memory_table = MemoryTable()
            memory_table.saveto(self.memory_table_path)

        # Cached table --

        # Set the path to the cached runs table
        self.cached_table_path = fs.join(self.analysis_path, cached_filename)

        # Initialize the cached table if necessary
        if not fs.is_file(self.cached_table_path):

            # Create the table and save it
            cached_table = CachedRunsTable()
            cached_table.saveto(self.cached_table_path)

        # Load the analysis runs object
        self.runs = AnalysisRuns(self.modeling_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_modeling_path(cls, modeling_path):

        """
        This function ...
        :param modeling_path:
        :return:
        """

        fit_path = fs.join(modeling_path, "analysis")
        return cls(fit_path)

    # -----------------------------------------------------------------

    @property
    def analysis_path(self):

        """
        This function ...
        :return:
        """

        return self.path

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This funtion ...
        :return:
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @lazyproperty
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return TimingTable.from_file(self.timing_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return MemoryTable.from_file(self.memory_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cached_table(self):

        """
        This function ...
        :return:
        """

        return CachedRunsTable.from_file(self.cached_table_path)

    # -----------------------------------------------------------------

    @property
    def cached_analysis_run_names(self):

        """
        This function ...
        :return:
        """

        return self.cached_table.run_names

    # -----------------------------------------------------------------

    @property
    def cache_host_ids(self):

        """
        This function ...
        :return:
        """

        return self.cached_table.cache_host_ids

    # -----------------------------------------------------------------

    def get_run_names_for_host_id(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.cached_table.run_names_for_host_id(host_id)

    # -----------------------------------------------------------------

    @property
    def analysis_run_names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.analysis_path, returns="name")

    # -----------------------------------------------------------------

    def get_run_path(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        path = fs.join(self.analysis_path, run_name)
        return path

    # -----------------------------------------------------------------

    def get_run_info(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        path = fs.join(self.get_run_path(run_name), "info.dat")
        return AnalysisRunInfo.from_file(path)

    # -----------------------------------------------------------------

    def get_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        info_path = fs.join(self.get_run_path(run_name), "info.dat")
        return AnalysisRun.from_info(info_path)

# -----------------------------------------------------------------
