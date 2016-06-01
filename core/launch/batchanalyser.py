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

        # Get the number of wavelengths
        wavelengths = self.log_file.wavelengths

        # Get the number of photon packages
        packages = self.ski.packages()

        # Get the number of dust cells
        cells = self.log_file.dustcells()

        # Check whether dust self-absorption was enabled for the simulation
        selfabsorption = self.ski.dustselfabsorption()

        # Check whether transient heating was enabled for the simulation
        transient_heating = self.ski.transientheating()

        # Get the different contributions to the simulation's runtime
        setup_time = self.te.setup
        stellar_time = self.te.stellar
        spectra_time = self.te.spectra
        dust_time = self.te.dust
        writing_time = self.te.writing
        waiting_time = self.te.waiting
        communication_time = self.te.communication
        intermediate_time = self.te.other

        # Open the timing table
        timing_table = TimingTable(self.timing_table_path)

        # Add an entry to the timing table
        # Simulation name, Timestamp, Host id, Cluster name, Cores, Hyperthreads per core, Processes, Wavelengths,
        # Packages, cells, Self-absorption, Total runtime, Setup time, Stellar emission time, Spectra calculation time,
        # Dust emission time, Writing time, Waiting time, Communication time, Intermediate time
        timing_table.add_entry(self.simulation.name, self.simulation.submitted_at, host_id, cluster_name, cores,
                               hyperthreads, processes, wavelengths, packages, cells, selfabsorption, transient_heating,
                               total_runtime, setup_time, stellar_time, spectra_time, dust_time, writing_time,
                               waiting_time, communication_time, intermediate_time)

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
        peak_memory_usage = None
        try: peak_memory_usage = self.log_file.peak_memory
        except RuntimeError:
            for log_file in self.simulation.logfiles():
                try:
                    peak_memory_usage = log_file.peak_memory
                    break
                except RuntimeError: pass
        if peak_memory_usage is None: raise RuntimeError("All log files were aborted")

        # Get the number of wavelengths
        input_path = self.simulation.input_path
        wavelengths = self.ski.nwavelengthsfile(input_path)

        # Get the number of dust cells
        dust_cells = self.log_file.dustcells()

        # Check whether dust self-absorption was enabled for the simulation
        selfabsorption = self.ski.dustselfabsorption()

        # Check whether transient heating was enabled for the simulation
        transient_heating = self.ski.transientheating()

        # Open the memory table
        memory_table = MemoryTable(self.memory_table_path)

        # Add an entry to the memory table
        memory_table.add_entry(self.simulation.name, self.simulation.submitted_at, host_id, cluster_name, cores,
                               hyperthreads, processes, wavelengths, dust_cells, selfabsorption, transient_heating,
                               peak_memory_usage)

# -----------------------------------------------------------------
