#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.progress Contains the ProgressPlotter class, used for creating plots of the progress
#  of different phases of a SKIRT simulation as a function of time.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .plotter import Plotter
from ..basics.map import Map
from ..basics.log import log
from ..tools import filesystem as fs

# -----------------------------------------------------------------

full_phase_names = {"stellar": "stellar emission phase",
                    "spectra": "calculation of dust emission spectra",
                    "dust": "dust emission phase"}

# -----------------------------------------------------------------

class ProgressPlotter(Plotter):
    
    """
    This class ...
    """
    
    def __init__(self):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ProgressPlotter, self).__init__()

    # -----------------------------------------------------------------

    @staticmethod
    def default_input():

        """
        This function ...
        :return:
        """

        return "progress.dat"

    # -----------------------------------------------------------------

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the input data into plottable format ...")

        # Get the number of processes
        ranks = np.unique(self.table["Process rank"])

        assert len(ranks) == max(ranks) + 1
        processes = len(ranks)

        # Initialize the data structure to contain the progress information in plottable format
        self.data = defaultdict(lambda: [Map({"times": [], "progress": []}) for i in range(processes)])

        # Loop over the different phases
        for phase in "stellar", "spectra", "dust":

            # Loop over the different entries in the progress table
            for i in range(len(self.table)):

                # Skip entries that do not belong to the current simulation phase
                if not self.table["Phase"][i] == phase: continue

                # Get the process rank
                rank = self.table["Process rank"][i]

                # Get the time and progress
                time = self.table["Time"][i]
                progress = self.table["Progress"][i]

                # Add the data point to the data structure
                self.data[phase][rank].times.append(time)
                self.data[phase][rank].progress.append(progress)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the plots ...")

        # Loop over the different phases in the data structure
        for phase in self.data:

            # Determine the path to the plot file for this phase
            plot_path = fs.join(self.output_path, "progress_" + phase + ".pdf")

            # Determine the title for the plot
            title = "Progress of " + full_phase_names[phase]

            # Create the plot for this simulation phase
            create_progress_plot(self.data[phase], plot_path, title)

# -----------------------------------------------------------------

def create_progress_plot(data, path, title):

    """
    This function ...
    :return:
    """

    # Initialize figure
    plt.figure()
    plt.clf()

    # Loop over all the different process ranks for which we have data
    for rank in range(len(data)):

        # Name of the current process
        process = "P" + str(rank)

        # Add the progress of the current process to the figure
        plt.plot(data[rank].times, data[rank].progress, label=process)

    plt.xlim(0)
    plt.grid('on')

    # Set the axis labels
    plt.xlabel("Time (s)", fontsize='large')
    plt.ylabel("Progress (\%)", fontsize='large')

    # Set the plot title
    plt.title(title)

    # Set the legend
    if len(data) > 16: plt.legend(loc='upper center', ncol=8, bbox_to_anchor=(0.5, -0.1), prop={'size': 8})
    elif len(data) > 1: plt.legend(loc='lower right', ncol=4, prop={'size': 8})
    else: pass

    # Save the figure
    plt.savefig(path, bbox_inches="tight", pad_inches=0.25)
    plt.close()
        
# -----------------------------------------------------------------
