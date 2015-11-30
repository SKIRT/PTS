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
import matplotlib

# Import astronomical modules
from astropy.table import Table


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
            #self.table = ascii.read("memory.dat", fill_values=fill_values)

        # Invalid input
        else: raise ValueError("Input must be either an Astropy Table object or a filename (e.g. memory.dat)")

        # Calculate the number of wavelengths and dust cells
        #Nlambda = skifile.nwavelengths()
        #Ncells = skifile.ncells()
        #Ndoubles = Nlambda * Ncells

        #plt.plot(seconds_list, memory_list)

        #totals = np.cumsum(deltas)

        #for i in range(len(times)):
        #    print times[i], totals[i]

        #plt.step(times, totals)
        #plt.fill_between(times, totals, color='green')
        #plt.bar(times, totals, color='r')
        #plt.show()

        pass

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Use a non-interactive back-end to generate high-quality vector graphics
        if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
        import matplotlib.pyplot as plt

# -----------------------------------------------------------------

