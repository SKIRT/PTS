#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.bestmodelanalyser Contains the BestModelAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools.logging import log
from ...core.launch.timing import TimingTable

# -----------------------------------------------------------------

class BestModelAnalyser(AnalysisComponent):

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
        super(AnalysisComponent, self).__init__(config)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The timeline extractor
        self.te = None

        # The log file of the simulation
        self.log_file = None

        # The ski file corresponding to the simulation
        self.ski = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline_extractor):

        """
        This function ...
        :param simulation:
        :param timeline_extractor:
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation, timeline_extractor)

        # 2. Load the log file of the simulation
        self.load_log_file()

        # 6. Write
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the best model analyser ...")

        # Set the attributes to default values
        self.simulation = None
        self.te = None
        self.log_file = None
        self.ski = None

    # -----------------------------------------------------------------

    def setup(self, simulation, timeline_extractor):

        """
        This function ...
        :param simulation:
        :param timeline_extractor:
        :return:
        """

        # Call the setup function of the base class
        super(BestModelAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

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
        self.write_timing()

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
        if cluster_name is None: cluster_name = "--"

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
