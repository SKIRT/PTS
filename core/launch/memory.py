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
from ..tools import tables
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class MemoryTable(object):

    """
    This class ...
    """

    def __init__(self, path):

        """
        This function ...
        """

        # Set the path of the memory table
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
        names = ["Simulation name", "Timestamp", "Host id", "Cluster name", "Cores", "Threads per core",
                 "Processes", "Wavelengths", "Dust cells", "Self-absorption", "Transient heating", "Data-parallel",
                 "Peak memory usage"]
        data = [[] for _ in names]
        dtypes = ["S24", "S23", "S15", "S15", "int64", "int64", "int64", "int64", "int64", "bool", "bool", "bool", "float64"]
        table = tables.new(data, names, dtypes=dtypes)

        # Set the column units
        table["Peak memory usage"] = "GB"  # memory usage is expressed in gigabytes

        # Write the (empty) table
        tables.write(table, self.path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    @classmethod
    def read(cls, path):

        """
        This function ...
        :return:
        """

        # Load the memory table
        memory_table = tables.from_file(path, format="ascii.ecsv")

        # Return the table
        return memory_table

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  dust_cells, selfabsorption, transient_heating, data_parallel, peak_memory_usage):

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
        :param peak_memory_usage:
        :return:
        """

        if cluster_name is None: cluster_name = "--"

        # Open the memory file in 'append' mode
        memory_file = open(self.path, 'a')

        # Initialize a list to contain the values of the row
        row = []

        # (12) Columns:
        # "Simulation name"
        # "Timestamp"
        # "Host id"
        # "Cluster name"
        # "Cores"
        # "Hyperthreads per core"
        # "Processes"
        # "Wavelengths"
        # "Dust cells"
        # "Self-absorption"
        # "Transient heating"
        # "Data-parallel"
        # "Peak memory usage"
        row.append(name)
        row.append(timestamp)
        row.append(host_id)
        row.append(cluster_name)
        row.append(str(cores))
        row.append(str(threads_per_core))
        row.append(str(processes))
        row.append(str(wavelengths))
        row.append(str(dust_cells))
        row.append(str(selfabsorption))
        row.append(str(transient_heating))
        row.append(str(data_parallel))
        row.append(str(peak_memory_usage))

        # Add the row to the runtime file
        memory_file.write(" ".join(row) + "\n")

        # Close the file
        memory_file.close()

# -----------------------------------------------------------------
