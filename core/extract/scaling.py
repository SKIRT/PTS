#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.extract.scaling Contains the ScalingExtractor class, used for extracting scaling information
#  from a simulation's log files.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.table import Table

# -----------------------------------------------------------------

class ScalingExtractor(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # -- Attributes --

        # The parallelization mode
        self.mode = None

        # The number of processes and threads
        self.processes = None
        self.threads = None

        # The path to the scaling file
        self.scaling_file_path = None

        # The scaling table
        self.table = None

        # The timeline and memory usage tables
        self.timeline = None
        self.memory = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline, memory):

        """
        This function ...
        :return:
        :param simulation:
        :param timeline:
        :param memory:
        """

        # Set the parallelization mode
        self.mode = simulation.analysis.scaling_run_name.split("__")[4]

        # Set the number of processes and threads
        self.processes = simulation.processes()
        self.threads = simulation.threads()

        # Set the path to the scaling file
        self.scaling_file_path = simulation.analysis.scaling_data_file

        # Cache local references to the timeline and memory usage tables
        self.timeline = timeline
        self.memory = memory

        # Write the relevant of the current simulation
        self.write()

        # Read in the extracted scaling table
        self.read()

        # Return the scaling table
        return self.table

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Open the output file
        resultfile = open(self.scaling_file_path, 'a')

        # Add a line to the output file containing the runtimes for the current simulation
        resultfile.write(self.mode + ' ' + str(self.processes) + ' ' + str(self.threads) + ' ' + str(self.timeline.setup)
                         + ' ' + str(self.timeline.stellar) + ' ' + str(self.timeline.spectra) + ' ' + str(self.timeline.dust)
                         + ' ' + str(self.timeline.writing) + ' ' + str(self.timeline.waiting) + ' ' + str(self.timeline.communication)
                         + ' ' + str(self.timeline.total) + ' ' + str(self.memory.peak) + '\n')

        # Close the output file
        resultfile.close()

    # -----------------------------------------------------------------

    def read(self):

        """
        This function ...
        :return:
        """

        # Read in the scaling data file
        self.table = Table.read(self.scaling_file_path, format="ascii.ecsv")

# -----------------------------------------------------------------
