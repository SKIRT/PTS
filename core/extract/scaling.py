#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.table import Table

# -----------------------------------------------------------------

class ScalingExtractor():

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        ## Attributes

        # The scaling table
        self.table = None

        # The timeline and memory extractors
        self.te = None
        self.me = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline_extractor, memory_extractor):

        """
        This function ...
        :return:
        """

        # Set the parallelization mode
        self.mode = simulation.scaling_run_name.split("__")[4]

        # Set the number of processes and threads
        self.processes = simulation.processes
        self.threads = simulation.threads

        # Set the scaling file path
        self.output_path = simulation.scaling_file_path

        # Cache local references to the timeline and memory extractors
        self.te = timeline_extractor
        self.me = memory_extractor

        # Write the relevant of the current simulation
        self.write()

        # Read in the extracted scaling table
        self.read()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Open the output file
        resultfile = open(self.output_path, 'a')

        # Add a line to the output file containing the runtimes for the current simulation
        resultfile.write(self.mode + ' ' + str(self.processes) + ' ' + str(self.threads) + ' ' + str(self.te.setup)
                         + ' ' + str(self.te.stellar) + ' ' + str(self.te.spectra) + ' ' + str(self.te.dust)
                         + ' ' + str(self.te.writing) + ' ' + str(self.te.waiting) + ' ' + str(self.te.communication)
                         + ' ' + str(self.te.total) + ' ' + str(self.me.peak) + '\n')

        # Close the output file
        resultfile.close()

    # -----------------------------------------------------------------

    def read(self):

        """
        This function ...
        :return:
        """

        # Read in the scaling data file
        self.table = Table.read(self.output_path)

# -----------------------------------------------------------------
