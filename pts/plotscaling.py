#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# This script creates a PDF plot showing the results of the parallel scaling benchmark on one or more target computers.

# -----------------------------------------------------------------

# Import standard modules
import numpy as np
import os.path
import multiprocessing

# -----------------------------------------------------------------

# Import the relevant PTS classes
from pts.log import Log

# -----------------------------------------------------------------

# Use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------

## An instance of the ScalingTest class represents a SKIRT scaling benchmark test for a particular ski file.
#
class ScalingPlotter:

    ## The constructor accepts the following arguments:
    #
    #  - directory: the path of the directory where the result files are located and where the plots should be created
    #  - system: the system for which to plot the statistics. If no system is given, plots are created with different
    #            curves for all the different systems (and modes) for which a results file is found.
    #
    def __init__(self, directory, system=""):

        # Create a logger
        self._log = Log()

        # Set the directory path
        self._directory = directory

        # Set the system name
        self._system = system

        # Create a list of result files (i.e. all *.dat files in the current directory),
        # ordered in a dictionary keyed on a combination of the system and mode in which the scaling test was run.
        # Each key in the dictionary will correspond to a different curve in the plots.
        filenames = dict()
        for filename in os.listdir(self._directory):

            # Check whether this file is a data file and not hidden
            if filename.endswith(".dat") and not filename.startswith("."):

                # Split the file name into its segments
                segments = filename.split("_")

                # Get the system name in which the scaling test was run for this results file
                systemname = segments[0]

                # If no system is specified, or the system corresponds with the system of the file, we proceed
                # In the former case, all files will pass through
                if not self._system or systemname == self._system:

                    # Get the mode in which the scaling test was run for this results file
                    mode = segments[1]

                    # Add the filename to the library at the correct key (the name of the system + the mode),
                    # if necessary, this will create a new key
                    filenames.setdefault((systemname,mode),[]).append(filename)

        # Show an error if no files were found (for the specified system)
        if len(filenames) == 0:
            self._log.error("No results file found for system " + self._system)

        # Create a dictionary to contain the statistics (timings with error bars for different thread counts)
        # with the keys equal to the keys of the filenames dictionary (systemname, mode).
        self._statistics = dict.fromkeys(filenames.keys(), 0)

        # For each different system and mode in which the scaling test has been performed
        for (systemname, mode), filelist in filenames.items():

            # A dictionary containing the runtimes for different amount of threads for this system and mode
            runtimes = dict()

            # For every results file in this run
            for file in filelist:

                # Get the values from this results file
                filepath = os.path.join(directory, file)
                threadcounts, times = np.loadtxt(filepath, usecols=(2,4), unpack=True)

                # Put these in a dictionary, with the keys being the number of threads
                for threadcount, time in zip(threadcounts, times):

                    runtimes.setdefault(threadcount,[]).append(time)

            # Now we want a list of the mean runtimes and their sample standard deviations
            meantimes = []
            errortimes = []

            # For each number of threads
            for threadcount, times in runtimes.items():

                # Calculate the mean runtime and add it to the list
                meantimes.append(np.mean(times))

                # Check whether we have more than one timing for this number of threads
                if len(times) > 1:

                    # If this is the case, calculate the sample standard deviation and add it to the list
                    errortimes.append(np.std(times, ddof=1))

                else:

                    # If only one timing was available, set the standard deviation to infinity.
                    # No error bars will be plotted for an data point with infinite standard deviation.
                    errortimes.append(float("inf"))

            # Make a list of the different threadcount
            nthreads = runtimes.keys()

            # Add the statistics for this system and mode
            self._statistics[(systemname,mode)] = [nthreads, meantimes, errortimes]

    # -----------------------------------------------------------------

    ## This function creates a PDF plot showing the execution time in function of thread count.
    # The function takes the following (optional) arguments:
    #
    # - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
    # - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
    # - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
    #
    def plottimes(self, figsize=(12,8), xlim=None, ylim=None):

        # Initialize a figure with the appropriate size
        plt.figure(figsize=figsize)

        # Loop over the different combinations of systemname and mode
        for (systemname, mode), stats in self._statistics.items():

            # Determine a label to identify this curve
            label = mode if self._system else systemname + " (" + mode + ")"

            # Plot the data points for this curve
            plt.errorbar(stats[0], stats[1], stats[2], marker='.', label=label)

        # Set axis limits if requested
        plt.grid(True)
        if xlim != None: plt.xlim(xlim)
        if ylim != None: plt.ylim(ylim)

        # Add axis labels and a legend
        plt.xlabel("Total number of threads", fontsize='large')
        plt.ylabel("Stellar emission time (s)", fontsize='large')
        plt.legend(title="Modes") if self._system else plt.legend(title="Systems")

        # Save the figure
        filename = "scaling_" + self._system + "_times.pdf" if self._system else "scaling_times.pdf"
        plotfilepath = os.path.join(self._directory, filename)
        plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    # -----------------------------------------------------------------

    ## This function creates a PDF plot showing the speedup in function of thread count.
    #  The speedup is defined as T(1)/T(N). It is a dimensionless quantity, theoretically it is >= 1.
    #  The function takes the following (optional) arguments:
    #
    # - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
    # - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
    # - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
    #
    def plotspeedups(self, figsize=(12,8), xlim=None, ylim=None):

        # Initialize a figure with the appropriate size
        plt.figure(figsize=figsize)

        # Loop over the different combinations of systemname and mode
        for (systemname, mode), stats in self._statistics.items():

            # Determine a label to identify this curve
            label = mode if self._system else systemname + " (" + mode + ")"

            # Calculate the speedups

            # TODO: check whether the first value is indeed for one thread
            serialtime = stats[1][0]
            speedups = serialtime / stats[1]

            # TODO: calculate the errors on the speedups

            # Plot the data points for this curve
            plt.errorbar(stats[0], speedups, stats[2], marker='.', label=label)

        # Set axis limits if requested
        plt.grid(True)
        if xlim != None: plt.xlim(xlim)
        if ylim != None: plt.ylim(ylim)

        # Add axis labels and a legend
        plt.xlabel("Total number of threads", fontsize='large')
        plt.ylabel("Speedup (dimensionless)", fontsize='large')
        plt.legend(title="Modes") if self._system else plt.legend(title="Systems")

        # Save the figure
        filename = "scaling_" + self._system + "_speedups.pdf" if self._system else "scaling_speedups.pdf"
        plotfilepath = os.path.join(self._directory, filename)
        plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    ## This function creates a PDF plot showing the efficiency in function of thread count.
    # Efficiency is defined as T(1)/T(N)/N. It is a dimensionless quantity <= 1.
    # The function takes the following (optional) arguments:
    #
    # - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
    # - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
    # - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
    #
    def ploteffs(self, figsize=(12,8), xlim=None, ylim=None):

        # Initialize a figure with the appropriate size
        plt.figure(figsize=figsize)

        # Loop over the different combinations of systemname and mode
        for (systemname, mode), stats in self._statistics.items():

            # Determine a label to identify this curve
            label = mode if self._system else systemname + " (" + mode + ")"

            # Calculate the efficiencies

            # TODO: check whether the first value is indeed for one thread
            serialtime = stats[1][0]
            efficiencies = (serialtime / stats[1]) / stats[0]

            # TODO: calculate the errors on the efficiencies

            # Plot the data points for this curve
            plt.errorbar(stats[0], efficiencies, stats[2], marker='.', label=label)

        # Set axis limits if requested
        plt.grid(True)
        if xlim != None: plt.xlim(xlim)
        if ylim != None: plt.ylim(ylim)

        # Add axis labels and a legend
        plt.xlabel("Total number of threads", fontsize='large')
        plt.ylabel("Efficiency (dimensionless)", fontsize='large')
        plt.legend(title="Modes") if self._system else plt.legend(title="Systems")

        # Save the figure
        filename = "scaling_" + self._system + "_efficiencies.pdf" if self._system else "scaling_efficiencies.pdf"
        plotfilepath = os.path.join(self._directory, filename)
        plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.25)
        plt.close()

# -----------------------------------------------------------------
