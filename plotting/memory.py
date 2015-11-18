#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plottimeline Plot a timeline of a simulation, for the different parallel processes, from
#  a text file containing the timeline information extracted from the simulation's log files.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os.path
import numpy as np

# Use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------

# Ignore warnings, otherwise Canopy would give a UserWarning on top of the error encountered when a progress
# file does not contain any data (an error which is catched an produces an error message).
import warnings
warnings.filterwarnings("ignore")

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

        pass

    # -----------------------------------------------------------------

    def run(self, input_path, output_path):

        """
        This function ...
        :return:
        """

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

