#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package plotting.plotprogress

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

# Import the relevant PTS classes and modules
from ..basics.map import Map

# -----------------------------------------------------------------

class ProgressPlotter(object):
    
    """
    This class ...
    """
    
    def __init__(self):
        
        """
        The constructor ...
        """

        # Set the table to None initially
        self.table = None

    # -----------------------------------------------------------------

    def run(self, input, output_path):
        
        """
        This function ...
        """

        # If the input is a Table object
        if isinstance(input, Table): self.table = input

        # If the input is a string
        elif isinstance(input, basestring):

            fill_values = ('--', '0', 'Simulation phase')
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

        # Get the number of processes
        ranks = np.unique(self.table["Process rank"])

        assert len(ranks) == max(ranks) + 1
        processes = len(ranks)

        # Loop over the different phases
        for phase in "stellar", "spectra", "dust":

            # Initialize a data structure to contain the progress information in plottable format
            data = [Map({"times": [], "progress": []}) for i in range(processes)]

            # Determine the path to the plot file for this phase
            plot_path = os.path.join(self.output_path, "progress_" + phase + ".pdf")

            # Determine the title for the plot
            title = phase + " progress"

            # Loop over the different entries in the progress table
            for i in range(len(self.table)):

                # Skip entries that do not belong to the current simulation phase
                if not self.table["Simulation phase"][i] == phase: continue

                # Get the process rank
                rank = self.table["Process rank"][i]

                # Get the time and progress
                time = self.table["Time"][i]
                progress = self.table["Progress"][i]

                # Add the data point to the data structure
                data[rank].times.append(time)
                data[rank].progress.append(progress)

            # Create the plot for this simulation phase
            self.create_plot(data, plot_path, title)

    # -----------------------------------------------------------------

    def create_plot(self, data, path, title):

        """
        This function ...
        :return:
        """

        # Create a PDF Pages object
        pp = PdfPages(path)

        # Initialize figure
        plt.figure()
        plt.clf()

        # Loop over all the different process ranks for which we have data
        for rank in data:

            # Name of the current process
            process = "P" + str(rank)

            # Add the progress of the current process to the figure
            plt.plot(data[rank].times, data[rank].progress, label=process)

        plt.xlim(0)
        plt.grid('on')
        plt.xlabel("Time (s)", fontsize='large')
        plt.ylabel("Progress (%)", fontsize='large')
        plt.title(title)

        if len(data) > 16: plt.legend(loc='upper center', ncol=8, bbox_to_anchor=(0.5,-0.1), prop={'size':8})
        else: plt.legend(loc='lower right', ncol=4, prop={'size':8})

        # Save the figure
        pp.savefig(bbox_inches="tight", pad_inches=0.25)
        pp.close()
        
# -----------------------------------------------------------------
