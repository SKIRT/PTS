#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.batchanalyser Contains the BatchAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable

# -----------------------------------------------------------------

class BatchAnalyser(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(BatchAnalyser, self).__init__(config, "core")

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The paths to the timing and memory table files
        self.timing_table_path = None
        self.memory_table_path = None

        # The timeline and memory extractors
        self.te = None
        self.me = None

        # The log file of the simulation
        self.log_file = None

        # The ski file corresponding to the simulation
        self.ski = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline_extractor, memory_extractor):

        """
        This function ...
        :param simulation:
        :param timeline_extractor:
        :param memory_extractor:
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation, timeline_extractor, memory_extractor)

        # 2. Load the log file of the simulation
        self.load_log_file()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the batch analyser ...")

        # Set the attributes to default values
        self.simulation = None
        self.timing_table_path = None
        self.memory_table_path = None
        self.te = None
        self.me = None
        self.log_file = None
        self.ski = None

    # -----------------------------------------------------------------

    def setup(self, simulation, timeline_extractor, memory_extractor):

        """
        This function ...
        :param simulation:
        :param timeline_extractor:
        :param memory_extractor:
        :return:
        """

        # Call the setup function of the base class
        super(BatchAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

        # Also make references to the timing and memory table files (for shortness of notation)
        self.timing_table_path = self.simulation.analysis.timing_table_path
        self.memory_table_path = self.simulation.analysis.memory_table_path

        # Make a reference to the timeline extractor
        self.te = timeline_extractor

        # Load the ski file
        self.ski = self.simulation.parameters()

    # -----------------------------------------------------------------

    def load_log_file(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulation log file ...")

        # Get the log file produced by the simulation
        self.log_file = self.simulation.log_file

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the timing information
        if self.timing_table_path is not None: self.write_timing()

        # Write the memory information
        if self.memory_table_path is not None: self.write_memory()

    # -----------------------------------------------------------------

    def write_timing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the timing information of the simulation ...")

        # Get the name of the host on which the simulation was run
        host_id = self.simulation.host_id
        cluster_name = self.simulation.cluster_name

        # Get the parallelization object from the simulation
        parallelization = self.simulation.parallelization

        # Get the paralleliation properties
        cores = parallelization.cores
        hyperthreads = parallelization.threads_per_core
        processes = parallelization.processes

        # Get the total runtime (in seconds)
        total_runtime = self.log_file.total_runtime

        # Get the number of photon packages
        packages = self.ski.packages()

        # Get the serial, parallel runtime and runtime overhead (in seconds)
        serial_runtime = self.te.serial
        parallel_runtime = self.te.parallel
        runtime_overhead = self.te.overhead

        # Open the timing table
        timing_table = TimingTable(self.timing_table_path)

        # Add an entry to the timing table
        timing_table.add_entry(self.simulation.name, self.simulation.submitted_at, host_id, cluster_name, cores,
                               hyperthreads, processes, packages, total_runtime, serial_runtime, parallel_runtime,
                               runtime_overhead)

    # -----------------------------------------------------------------

    def write_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the memory usage information of the simulation ...")

        # Get the name of the host on which the simulation was run
        host_id = self.simulation.host_id
        cluster_name = self.simulation.cluster_name

        # Get the parallelization object from the simulation
        parallelization = self.simulation.parallelization

        # Get the paralleliation properties
        cores = parallelization.cores
        hyperthreads = parallelization.threads_per_core
        processes = parallelization.processes

        # Get the peak memory usage
        peak_memory_usage = self.log_file.peak_memory

        # Get the number of wavelengths
        input_path = self.simulation.input_path
        wavelengths = self.ski.nwavelengthsfile(input_path)

        # Open the memory table
        memory_table = MemoryTable(self.memory_table_path)

        # Add an entry to the memory table
        memory_table.add_entry(self.simulation.name, self.simulation.submitted_at, host_id, cluster_name, cores,
                               hyperthreads, processes, wavelengths, peak_memory_usage)

# -----------------------------------------------------------------
