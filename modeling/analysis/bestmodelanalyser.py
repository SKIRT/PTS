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

        # Get the runtime in seconds
        runtime = self.log_file.total_runtime

        # Get the number of photon packages
        packages = self.ski.packages()

        # Open the timing file in 'append' mode
        timing_file = open(self.timing_table_path, 'a')

        # Initialize a list to contain the values of the row
        row = []

        # Columns:
        # "Submission time"
        # "Host id"
        # "Cluster name"
        # "Cores"
        # "Hyperthreads per core"
        # "Processes"
        # "Packages"
        # "Total runtime"
        # "Serial runtime"
        # "Parallel runtime"
        # "Overhead runtime"
        row.append(self.simulation.submitted_at)
        row.append(host_id)
        row.append(cluster_name)
        row.append(str(cores))
        row.append(str(hyperthreads))
        row.append(str(processes))
        row.append(str(packages))
        row.append(str(runtime))
        row.append(str(self.te.serial))
        row.append(str(self.te.parallel))
        row.append(str(self.te.overhead))

        # Add the row to the table file
        timing_file.write(" ".join(row) + "\n")

        # Close the file
        timing_file.close()

# -----------------------------------------------------------------
