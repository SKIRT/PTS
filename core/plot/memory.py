#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module ...
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

# Import the relevant PTS classes and modules


# -----------------------------------------------------------------

class MemoryPlotter(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Set the table to None initially
        self.table = None

    # -----------------------------------------------------------------

    def run(self, input, output_path):

        """
        This function ...
        :return:
        """

        # If the input is a Table object
        if isinstance(input, Table): self.table = input

        # If the input is a string
        elif isinstance(input, basestring):

            fill_values = [('--', '0', 'Simulation phase'), ('--', '0', 'Array (de)allocation'), ('--', '0', 'Array ID')]
            self.table = ascii.read(input, fill_values=fill_values)

        # Invalid input
        else: raise ValueError("Input must be either an Astropy Table object or a filename (e.g. memory.dat)")

        # Create the plots
        self.plot(output_path)

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :return:
        """

        # Create a PDF Pages object
        pp = PdfPages(path)

        # Initialize figure
        plt.figure()
        plt.clf()

        # Plot the memory usage
        plt.plot(self.table["Simulation time"], self.table["Memory usage"])

        if "Array (de)allocation" in self.table.colnames:

            # Calculate the cumulative allocated memory
            totals = np.cumsum(self.table["Array (de)allocation"])

            #plt.step(times, totals)
            plt.fill_between(self.table["Simulation time"], totals, color='green')
            #plt.bar(times, totals, color='r')

        # Save the figure
        pp.savefig()
        pp.close()

# -----------------------------------------------------------------

