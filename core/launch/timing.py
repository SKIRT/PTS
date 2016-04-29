#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.parallelization Contains the Parallelization class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import tables
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class TimingTable(object):

    """
    This class ...
    """

    def __init__(self, path):

        """
        This function ...
        """

        # Set the path of the timing table
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
                 "Processes", "Packages", "Total runtime", "Serial runtime", "Parallel runtime", "Runtime overhead"]
        data = [[], [], [], [], [], [], [], [], [], [], [], []]
        dtypes = ["S24", "S23" "S15", "S15", "int64", "int64", "int64", "int64", "float64", "float64", "float64", "float64"]
        table = tables.new(data, names, dtypes=dtypes)

        # Set the column units
        table["Total runtime"] = "s"  # runtimes are expressed in seconds
        table["Serial runtime"] = "s"
        table["Parallel runtime"] = "s"
        table["Runtime overhead"] = "s"

        # Write the (empty) table
        tables.write(table, self.path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    @classmethod
    def read(cls, path):

        """
        This function ...
        :return:
        """

        # Load the timing table
        timing_table = tables.from_file(path, format="ascii.ecsv")

        # Return the table
        return timing_table

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, packages,
                  total_runtime, serial_runtime, parallel_runtime, runtime_overhead):

        """
        This function ...
        :param name:
        :param timestamp:
        :param host_id:
        :param cluster_name:
        :param cores:
        :param threads_per_core:
        :param processes:
        :param packages:
        :param total_runtime:
        :param serial_runtime:
        :param parallel_runtime:
        :param runtime_overhead:
        :return:
        """

        if cluster_name is None: cluster_name = "--"

        # Open the timing file in 'append' mode
        timing_file = open(self.path, 'a')

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
        # "Packages"
        # "Total runtime"
        # "Serial runtime"
        # "Parallel runtime"
        # "Runtime overhead"
        row.append(name)
        row.append(timestamp)
        row.append(host_id)
        row.append(cluster_name)
        row.append(str(cores))
        row.append(str(threads_per_core))
        row.append(str(processes))
        row.append(str(packages))
        row.append(str(total_runtime))
        row.append(str(serial_runtime))
        row.append(str(parallel_runtime))
        row.append(str(runtime_overhead))

        # Add the row to the runtime file
        timing_file.write(" ".join(row) + "\n")

        # Close the file
        timing_file.close()

# -----------------------------------------------------------------
