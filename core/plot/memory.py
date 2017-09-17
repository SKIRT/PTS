#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.memory Contains the MemoryPlotter class, used for creating plots of the memory consumption
#  of a SKIRT simulation as a function of time.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ..basics.map import Map
from .plotter import Plotter
from ..basics.log import log
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class MemoryPlotter(Plotter):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(MemoryPlotter, self).__init__()

        # -- Attributes --

        # A data structure to store the memory (de)allocation information
        self.allocation = None

    # -----------------------------------------------------------------

    @staticmethod
    def default_input():

        """
        This function ...
        :return:
        """

        return "memory.dat"

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

        # Initialize the data structure to contain the memory usage information in plottable format
        self.data = [Map({"times": [], "memory": []}) for i in range(processes)]

        # Loop over the different entries in the memory table
        for i in range(len(self.table)):

            # Get the process rank
            rank = self.table["Process rank"][i]

            # Get the time and memory usage
            time = self.table["Simulation time"][i]
            memory = self.table["Memory usage"][i]

            # Add the data point to the data structure
            self.data[rank].times.append(time)
            self.data[rank].memory.append(memory)

        # Check whether (de)allocation information is present in the memory table
        if "Array (de)allocation" in self.table.colnames:

            # Initialize the data structure for plotting the memory usage of the root process and the memory
            # allocation curve
            self.allocation = Map({"times": [], "allocation": [], "cumulative": []})

            # Get the mask covering entries that do not contain array (de)allocation information
            mask = self.table["Array (de)allocation"].mask

            # Check whether the first entry of the table corresponds to the root process
            assert self.table["Process rank"][0] == 0

            # Create a variable to store the cumulative sum of allocated memory
            cumulative_sum = 0.0

            # Loop over the different entries in the memory table
            for i in range(len(self.table)):

                # Get the process rank
                rank = self.table["Process rank"][i]

                # Only add the contributions from the root process
                if rank > 0: break

                # If the entry is masked because it does not contain memory allocation information, skip it
                if mask[i]: continue

                # Get the time and the amount of (de)allocated memory
                time = self.table["Simulation time"][i]
                allocation = self.table["Array (de)allocation"][i]

                # Add the allocated memory to the sum
                cumulative_sum += allocation

                # Add the data point to the data structure
                self.allocation.times.append(time)
                self.allocation.allocation.append(allocation)
                self.allocation.cumulative.append(cumulative_sum)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the plots ...")

        # Make a plot of the memory usage as a function of time
        self.plot_memory()

        # Make a plot of the memory (de)allocation information, if present
        if self.allocation is not None: self.plot_allocation()

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the plot file
        plot_path = fs.join(self.output_path, "memory.pdf")

        # Initialize figure
        plt.figure()
        plt.clf()

        # Loop over the different processes
        for rank in range(len(self.data)):

            # Name of the current process
            process = "P" + str(rank)

            # Plot the memory usage
            plt.plot(self.data[rank].times, self.data[rank].memory, label=process)

        # Set the axis labels
        plt.xlabel("Time (s)", fontsize='large')
        plt.ylabel("Memory usage (GB)", fontsize='large')

        # Set the plot title
        plt.title("Memory consumption")

        # Set the legend
        if len(self.data) > 16: plt.legend(loc='upper center', ncol=8, bbox_to_anchor=(0.5, -0.1), prop={'size': 8})
        else: plt.legend(loc='lower right', ncol=4, prop={'size': 8})

        # Save the figure
        plt.savefig(plot_path, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    # -----------------------------------------------------------------

    def plot_allocation(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the plot file
        plot_path = fs.join(self.output_path, "allocation.pdf")

        # Initialize figure
        plt.figure()
        plt.clf()

        # Plot the memory usage of the root process
        plt.plot(self.data[0].times, self.data[0].memory, label="total memory usage")

        # Plot the memory allocation of the root process
        plt.step(self.allocation.times, self.allocation.cumulative, where="post", linestyle="--", label="allocated array memory")

        # Set the axis labels
        plt.xlabel("Time (s)", fontsize='large')
        plt.ylabel("Memory usage (GB)", fontsize='large')

        # Set the plot title
        plt.title("Memory (de)allocation")

        # Set the legend
        plt.legend(loc='lower right', prop={'size': 8})

        # Save the figure
        plt.savefig(plot_path, bbox_inches='tight', pad_inches=0.25)
        plt.close()

# -----------------------------------------------------------------
