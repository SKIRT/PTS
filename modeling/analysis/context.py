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

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.utils import lazyproperty
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .run import AnalysisRuns
from .run import AnalysisRunInfo, AnalysisRun, CachedAnalysisRun, CachedAnalysisRuns
from .tables import CachedRunsTable
from ...core.remote.remote import Remote
from ...core.tools import tables

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

    @property
    def modeling_data_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.modeling_path, "data")

    # -----------------------------------------------------------------

    @property
    def galaxy_info_path(self):

        """
        This function ...
        :return:
        """

        #  Set the path to the galaxy info file
        return fs.join(self.modeling_data_path, "info.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_info(self):

        """
        This function ...
        :return:
        """

        # Load the info table
        table = tables.from_file(self.galaxy_info_path)

        # To ordered dict
        info = OrderedDict()
        for name in table.colnames: info[name] = table[name][0]

        # Return the info
        return info

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_type(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_info["Hubble Type"]

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_stage(self):

        """
        Thisf unction ...
        :return:
        """

        return self.galaxy_info["Hubble Stage"]

    # -----------------------------------------------------------------

    @property
    def galaxy_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.modeling_path)

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

    def has_cached_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return host_id in self.cache_host_ids

    # -----------------------------------------------------------------

    def get_run_names_for_host_id(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.cached_table.run_names_for_host_id(host_id)

    # -----------------------------------------------------------------

    def get_host_id_for_run_name(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.cached_table.host_id_for_run_name(run_name)

    # -----------------------------------------------------------------

    @property
    def analysis_run_names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.analysis_path, returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def all_analysis_run_names(self):

        """
        This function ...
        :return:
        """

        return self.cached_analysis_run_names + self.analysis_run_names

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

    def get_run_info_path(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        path = fs.join(self.get_run_path(run_name), "info.dat")
        return path

    # -----------------------------------------------------------------

    def get_run_info(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        path = self.get_run_info_path(run_name)
        return AnalysisRunInfo.from_file(path)

    # -----------------------------------------------------------------

    def get_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        info_path = fs.join(self.get_run_path(run_name), "info.dat")
        return AnalysisRun.from_info(info_path, hubble_stage=self.hubble_stage)

    # -----------------------------------------------------------------

    def get_remote_for_run_name(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        # Get host ID
        host_id = self.get_host_id_for_run_name(run_name)

        # Connect and return
        remote = Remote(host_id=host_id, silent=True)
        return remote

    # -----------------------------------------------------------------

    def get_cached_run_info(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        # Get the remote
        remote = self.get_remote_for_run_name(run_name)

        # Determine name of the remote analysis cache directory
        cache_directory_name = self.galaxy_name + "_analysis"

        # Get the cache path
        cache_path = fs.join(remote.home_directory, cache_directory_name)

        # Get the analysis run path on the remote
        run_path = fs.join(cache_path, run_name)

        # Get the info file path
        info_path = fs.join(run_path, "info.dat")

        # Return the info
        return AnalysisRunInfo.from_remote_file(info_path, remote)

    # -----------------------------------------------------------------

    def get_cached_runs_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Return the analysis runs object
        return CachedAnalysisRuns(self.modeling_path, host_id)

    # -----------------------------------------------------------------

    @lazyproperty
    def cached_runs(self):

        """
        Thins function ...
        :return:
        """

        # Initialize a dictionary to contain the runs per remote
        runs = dict()

        # Loop over the host IDs
        for host_id in self.cached_table.cache_host_ids:

            # Load the runs
            runs[host_id] = self.get_cached_runs_for_remote(host_id)

        # Return the runs dict
        return runs

    # -----------------------------------------------------------------

    def get_cached_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        # Get the remote
        remote = self.get_remote_for_run_name(run_name)

        # Determine name of the remote analysis cache directory
        cache_directory_name = self.galaxy_name + "_analysis"

        # Get the cache path
        cache_path = fs.join(remote.home_directory, cache_directory_name)

        # Get the analysis run path on the remote
        run_path = fs.join(cache_path, run_name)

        # Create the cached analysis run and return it
        return CachedAnalysisRun(run_path, remote)

# -----------------------------------------------------------------
