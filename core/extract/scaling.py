#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

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

        self.te = None
        self.me = None

    # -----------------------------------------------------------------

    def run(self, simulation, timeline_extractor, memory_extractor, output_path):

        """
        This function ...
        :return:
        """

        # Set the number of processes and threads
        self.processes = simulation.processes
        self.threads = simulation.threads

        # Cache local references to the timeline and memory extractors
        self.te = timeline_extractor
        self.me = memory_extractor

        # Write the results
        self.write(output_path)

    # -----------------------------------------------------------------

    def write(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Open the output file
        resultfile = open(output_path, 'a')

        # Add a line to the output file containing the runtimes for the current simulation
        resultfile.write(str(self.processes) + ' ' + str(self.threads) + ' ' + str(self.processes*self.threads) + ' '
                         + str(self.te.setup) + ' ' + str(self.te.stellar) + ' ' + str(self.te.spectra)
                         + ' ' + str(self.te.dust) + ' ' + str(self.te.writing) + ' ' + str(self.te.waiting) + ' '
                         + str(self.te.communication) + ' ' + str(self.te.total) + ' ' + str(self.me.peak) + '\n')

# -----------------------------------------------------------------
