#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.timing Contains the TimingTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.table import SmartTable

# -----------------------------------------------------------------

class TimingTable(SmartTable):

    """
    This class ...
    """

    @classmethod
    def initialize(cls):

        """
        This function ...
        :return:
        """

        # Total time: Total simulation time (s)
        # Setup time: Time spent in simulation setup (s)
        # Stellar emission time: Time spent shooting stellar photon packages (s)
        # Spectra calculation time: Time spent in calculation of dust emission spectra (s)
        # Dust emission time: Time spent shooting dust emission photon packages (s)
        # Writing time: Time spent writing to disk (s)
        # Waiting time: Time spent waiting for other processes (s)
        # Communication time: Time spent in inter-process communication (s)
        # Intermediate time: Time spent in between other phases (s)

        # Create the table
        names = ["Simulation name", "Timestamp", "Host id", "Cluster name", "Cores", "Threads per core",
                 "Processes", "Wavelengths", "Packages", "Dust cells", "Self-absorption", "Transient heating",
                 "Data-parallel", "Total runtime", "Setup time", "Stellar emission time", "Spectra calculation time",
                 "Dust emission time", "Writing time", "Waiting time", "Communication time", "Intermediate time"]
        #data = [[] for _ in names]
        #dtypes = ["S24", "S23", "S15", "S15", "int64", "int64", "int64", "int64", "int64", "int64", "bool", "bool",
        #          "bool", "float64", "float64", "float64", "float64", "float64", "float64", "float64", "float64",
        #          "float64"]
        dtypes = [str, str, str, str, int, int, int, int, int, int, bool, bool, bool, float, float, float, float, float, float, float, float]

        # Call the constructor of the base class
        table = cls(names=names, dtype=dtypes, masked=True)

        # Set the column units
        table["Total runtime"] = "s"  # runtimes are expressed in seconds
        table["Setup time"] = "s"
        table["Stellar emission time"] = "s"
        table["Spectra calculation time"] = "s"
        table["Dust emission time"] = "s"
        table["Writing time"] = "s"
        table["Waiting time"] = "s"
        table["Communication time"] = "s"
        table["Intermediate time"] = "s"

        # Add the path attribute
        table.path = None

        # Return the timing table instance
        return table

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        # Open the table
        table = super(TimingTable, cls).read(path, format="ascii.ecsv")

        # Set the path
        table.path = path

        # Return the table
        return table

    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths,
                  packages, cells, selfabsorption, transient_heating, data_parallel, total_runtime, setup_time,
                  stellar_time, spectra_time, dust_time, writing_time, waiting_time, communication_time,
                  intermediate_time):

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
        :param packages:
        :param cells:
        :param selfabsorption:
        :param transient_heating:
        :param data_parallel:
        :param total_runtime:
        :param setup_time:
        :param stellar_time:
        :param spectra_time:
        :param dust_time:
        :param writing_time:
        :param waiting_time:
        :param communication_time:
        :param intermediate_time:
        :return:
        """

        #if cluster_name is None: cluster_name = "--"

        # COLUMNS:
        #names = ["Simulation name", "Timestamp", "Host id", "Cluster name", "Cores", "Threads per core",
        #         "Processes", "Wavelengths", "Packages", "Dust cells", "Self-absorption", "Transient heating",
        #         "Data-parallel", "Total runtime", "Setup time", "Stellar emission time", "Spectra calculation time",
        #         "Dust emission time", "Writing time", "Waiting time", "Communication time", "Intermediate time"]

        values = [name, timestamp, host_id, cluster_name, cores, threads_per_core, processes, wavelengths, packages,
                  cells, selfabsorption, transient_heating, data_parallel, total_runtime, setup_time, stellar_time,
                  spectra_time, dust_time, writing_time, waiting_time, communication_time, intermediate_time]

        # Resize string columns for longer entries
        self._resize_string_columns(values)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
