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
import os
import numpy as np
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

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

        # Set the path to the output directory
        self.output_path = output_path

        # Create the plots
        self.plot()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the plot file
        plot_path = os.path.join(self.output_path, "memory.pdf")

        # Initialize figure
        plt.figure()
        plt.clf()

        # Plot the memory usage
        plt.plot(self.table["Simulation time"], self.table["Memory usage"])

        # Check whether (de)allocation information is present in the memory table
        if "Array (de)allocation" in self.table.colnames:

            # Get the mask covering entries that do not contain array (de)allocation information
            mask = self.table["Array (de)allocation"].mask

            # Get the compressed (de)allocation data
            allocation = self.table["Array (de)allocation"].compressed()

            # Get the compressed time data
            times = np.ma.masked_array(self.table["Simulation time"], mask=mask).compressed()

            # Calculate the cumulative allocated memory
            totals = np.cumsum(allocation)

            plt.step(times, totals, where='post', linestyle='--')  # or linestyle=':'
            #plt.fill_between(self.table["Simulation time"], 0, totals, color='green')
            #plt.bar(times, totals, color='r')

        # Save the figure
        plt.savefig(plot_path, bbox_inches="tight")
        plt.close()

# -----------------------------------------------------------------
