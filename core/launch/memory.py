#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.memory Contains the MemoryTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.table import SmartTable

# -----------------------------------------------------------------

class MemoryTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Simulation name", str, None, "name of the simulation"),
                   ("Timestamp", str, None, "timestamp"),
                   ("Host id", str, None, "remote host ID"),
                   ("Cluster name", str, None, "remote cluster name"),
                   ("Cores", int, None, "number of cores"),
                   ("Threads per core", int, None, "number of threads per core"),
                   ("Processes", int, None, "number of processes"),
                   ("Wavelengths", int, None, "number of wavelengths"),
                   ("Dust cells", int, None, "number of dust cells"),
                   ("Self-absorption", bool, None, "self-absorption enabled"),
                   ("Transient heating", bool, None, "transient (non-LTE) heating enabled"),
                   ("Data-parallel", bool, None, "data parallelization enabled"),
                   ("Number of pixels", int, None, "total number of spatial pixels for all instruments"),
                   ("Peak memory usage", float, "GB", "peak memory usage")]

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  dust_cells, selfabsorption, transient_heating, data_parallel, npixels, peak_memory_usage):

        """
        This function ...
        :param name:
        :param timestamp:
        :param host_id:
        :param cluster_name:
        :param cores:
        :param threads_per_core:
        :param processes:
        :param wavelengths:
        :param dust_cells:
        :param selfabsorption:
        :param transient_heating:
        :param data_parallel:
        :param npixels
        :param peak_memory_usage:
        :return:
        """

        values = [name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  dust_cells, selfabsorption, transient_heating, data_parallel, npixels, peak_memory_usage]

        # Resize string columns for longer entries
        self._resize_string_columns(values)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
