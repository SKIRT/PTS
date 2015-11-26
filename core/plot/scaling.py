#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plotscaling The class ScalingPlotter in this module makes plots of the results SKIRT
#  scaling benchmark tests performed with the scalingtest module.

# -----------------------------------------------------------------

# Import standard modules
import math
import numpy as np
import os.path
from collections import defaultdict
from scipy.optimize import curve_fit

# Import the relevant PTS classes and modules
from ..basics import Log

# -----------------------------------------------------------------

# Use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

import matplotlib.ticker

# -----------------------------------------------------------------

# Ignore warnings, otherwise Canopy would give a UserWarning on top of the error encountered when a scaling
# test results file does not contain any data (an error which is catched and produces an error message).
import warnings
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------

# Define which column in the scaling test results file defines which quantity
columns = {'threads': 2, 'setup': 3, 'stellar': 4, 'dustselfabs': 5, 'dustem': 6, 'writing': 7, 'total': 8}

# For each different phase, define its full name
fullnames = {'setup': 'Setup', 'stellar': 'Stellar emission', 'dustselfabs': 'Dust self-absorption',
             'dustem': 'Dust emission', 'writing': 'Writing', 'total': 'Total simulation'}

# -----------------------------------------------------------------

## An instance of the ScalingPlotter is used to create plots of the runtimes, speedups and efficiencies as a function
#  of the number of threads, based on the results of (multiple) SKIRT scaling benchmark tests.
#
class ScalingPlotter(object):

    ## The constructor accepts the following arguments:
    #
    #  - directory: the path of the directory where the result files are located and where the plots should be created
    #  - phase: the phase of the simulation for which the scaling should be plotted. This argument can take the
    #           following values:
    #               * 'setup': for the setup of the simulation
    #               * 'stellar': for the stellar emission phase
    #               * 'dustselfabs': for the dust selfabsorption phase
    #               * 'dustem': for the dust emission phase
    #               * 'writing': for the writing phase
    #               * 'total': for the total simulation
    #  - system: the system for which to plot the statistics. If no system is given, plots are created with different
    #            curves for all the different systems (and modes) for which a results file is found.
    #
    def __init__(self, directory, phase, system=""):

        # Create a logger
        self._log = Log()

        # Set the results and visualization paths
        self._respath = os.path.join(directory, "res")
        self._vispath = os.path.join(directory, "vis")

        # Set the phase
        self._phase = phase

        # Set the system name
        self._system = system

        # Create a list of all scaling results files ('scaling.dat') found in a subdirectory of the results path,
        # ordered in a dictionary keyed on a combination of the system and mode in which the scaling test was run.
        # Each key in the dictionary will correspond to a different curve in the plots.
        filepaths = defaultdict(list)
        self._log.info("Gathering the data files...")
        for itemname in os.listdir(self._respath):

            # Define the full path to this item
            itempath = os.path.join(self._respath, itemname)

            # Check whether this item is a directory and it is not hidden
            if os.path.isdir(itempath) and not itemname.startswith("."):

                # Split the file name into its segments
                segments = itemname.split("_")

                # Get the system name in which the scaling test was run for this results file
                systemname = segments[0]

                # If no system is specified, or the system corresponds with the system of the file, we proceed
                # In the former case, all files will pass through
                if not self._system or systemname == self._system:

                    # Get the mode in which the scaling test was run for this results file
                    mode = segments[1]

                    # Define the name of the scaling results file inside this directory
                    filepath = os.path.join(itempath, "scaling.dat")

                    # If this file exists, add its path to the library at the correct key (the name of the
                    # system + the mode), if necessary, this will create a new key
                    if os.path.isfile(filepath): filepaths[(systemname,mode)].append(filepath)

        # Show an error if no files were found (for the specified system)
        if len(filepaths) == 0: self._log.error("No results file found for system " + self._system)

        # Create a dictionary to contain the statistics (timings with error bars for different thread counts)
        # with the keys equal to the keys of the filenames dictionary (systemname, mode).
        self._statistics = dict.fromkeys(filepaths.keys(), 0)

        # For each different system and mode in which the scaling test has been performed
        for (systemname, mode), filelist in filepaths.items():

            # A dictionary containing the runtimes for different amount of threads for this system and mode
            runtimes = defaultdict(list)

            # For every scaling results file in this run
            for filepath in filelist:

                # Get the values from this results file. If no data could be found in the file, we skip it.
                try:
                    threadcounts, times = np.loadtxt(filepath, usecols=(columns['threads'],columns[phase]), unpack=True)
                except ValueError:
                    self._log.warning("The file " + filepath + " does not contain any data")
                    continue

                # If there was only one entry in the results file, make lists of the number of threads and the runtime
                if isinstance(threadcounts, float):
                    threadcounts = [threadcounts]
                    times = [times]

                # Put these in a dictionary, with the keys being the number of threads
                for threadcount, time in zip(threadcounts, times):

                    runtimes[threadcount].append(time)

            # Check whether any valid data file was found for this system and mode. Otherwise, show a warning
            # and remove the corresponding key from the dictionary.
            if not runtimes:
                self._log.warning("No valid data was found for " + systemname + " in " + mode + " mode")
                self._statistics.pop((systemname, mode))
                continue

            # Make a list of the different threadcounts, order them from lowest to highest
            nthreads = sorted(runtimes.keys())

            # Now we want a list of the mean runtimes and a list of their standard deviations,
            # corresponding with the threadcounts in the 'nthreads' list.
            meantimes = []
            errortimes = []

            self._log.info("Calculating the average runtimes and standard deviations for " + systemname + " in " + mode + " mode...")

            # For each number of threads (lowest to highest)
            for threadcount in nthreads:

                # Get the timings for this threadcount from the dictionary
                times = runtimes[threadcount]

                # Calculate the mean runtime for this number of threads
                meantime = np.mean(times)

                # Add the mean runtime to the list
                meantimes.append(meantime)

                # Check whether we have more than one timing for this number of threads
                if len(times) > 1:

                    # If this is the case, calculate the sample standard deviation and add it to the list
                    errortimes.append(np.std(times, ddof=1))

                else:

                    # If only one timing was available, set the standard deviation to infinity.
                    # No error bars will be plotted for an data point with infinite standard deviation.
                    errortimes.append(float("inf"))

            # Add the statistics for this system and mode
            self._statistics[(systemname,mode)] = [nthreads, meantimes, errortimes]

        # Make dictionaries to store the average serial runtime and its standard deviation for each system
        self._serialruntime = dict()
        self._serialerror = dict()

        # Use a dictionary to temporarily store the serial runtimes for each different system
        serialruntimes = defaultdict(list)

        # For each system, average the serial (1 process, 1 threads) runtimes over the different modes
        for (systemname, mode), [nthreads, times, _] in self._statistics.items():

            try:

                # Try to find the index of the entry corresponding with a total number of threads equal to 1
                serialindex = nthreads.index(1)

                # For this serial run, get the runtime
                serialtime = times[serialindex]

                # Add this serial runtime to the dictionary, corresponding to the appropriate system name
                serialruntimes[systemname].append(serialtime)

            except ValueError:
                    # Try the next (systemname, mode): no serial runtime could be found for this combination
                    pass

        # If the dictionary with serial runtimes is empty, exit with an error
        if not serialruntimes:

            self._log.error("No serial runtimes could be found for the simulation; required to normalize the speedups"
                            " and the efficiencies")
            exit()

        # For each system, calculate the average serial runtime and the standard deviation
        for systemname, serialtimes in serialruntimes.items():

            # Try to calculate the mean
            mean = np.mean(serialtimes)

            # If the result is NaN, the simulation was never run on one thread for this system
            if np.isnan(mean):

                # Exit with an error message
                self._log.error("Could not find the " + phase + " runtime for one thread on this system")
                exit()

            # Check whether we have more than one serial timing for this system
            if len(serialtimes) > 1:

                # If this is the case, calculate the sample standard deviation
                error = np.std(serialtimes, ddof=1)

            else:

                # If only one timing was available, set the standard deviation to infinity.
                error = float("inf")

            # Add the mean value and the standard deviation to the dictionary
            self._serialruntime[systemname] = mean
            self._serialerror[systemname] = error

        # For each system, replace the serial runtime in every mode by the mean average runtime for that system
        for (systemname, mode), [nthreads, _, _] in self._statistics.items():

            try:

                # Try to find the index of the entry corresponding with a total number of threads equal to 1
                serialindex = nthreads.index(1)

                # Replace the serial runtimes for this system and this mode
                self._statistics[(systemname, mode)][1][serialindex] = self._serialruntime[systemname]
                self._statistics[(systemname, mode)][2][serialindex] = self._serialerror[systemname]

            except ValueError:

                # For this configuration no runtime was recorded with only one thread, so we add an entry
                # with the appropriate serial runtime and its error, if the mode is 'threads' or 'mpi' and not hybrid
                # (because the curve for hybrid should only start for a number of processors = threads per process)
                # Basically, this procedure lets us use the average runtime on 1 processor as part of an 'mpi' mode
                # scaling test as the runtime on 1 processor for the 'threads' mode scaling relation, if we didn't
                # have this runtime available from the 'threads' mode scaling test results file
                if mode == "threads" or mode == "mpi":

                    self._statistics[(systemname, mode)][0].insert(0, 1)
                    self._statistics[(systemname, mode)][1].insert(0, self._serialruntime[systemname])
                    self._statistics[(systemname, mode)][2].insert(0, self._serialerror[systemname])

    # -----------------------------------------------------------------

    ## This function creates a PDF plot showing the execution time in function of thread count.
    #  The function takes the following (optional) arguments:
    #
    #  - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
    #  - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
    #  - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
    #
    def plottimes(self, figsize=(12,8), xlim=None, ylim=None):

        # Inform the user of the fact that the runtimes are being plotted
        self._log.info("Plotting the runtimes...")

        # Initialize a figure with the appropriate size
        plt.figure(figsize=figsize)

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different combinations of systemname and mode (i.e. different curves)
        for (systemname, mode), [nthreads, times, errors] in self._statistics.items():

            # Determine a label to identify this curve
            if "hybrid" in mode: mode = "hybrid " + mode.strip('hybrid') + ":1"
            label = mode if self._system else systemname + " (" + mode + ")"

            # Plot the data points for this curve
            plt.errorbar(nthreads, times, errors, marker='.', label=label)

            # Add the appropriate ticks
            ticks |= set(nthreads)

        # Use a logarithmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1]*2)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        plt.grid(True)

        # Set custom axis limits if requested
        if xlim != None: plt.xlim(xlim)
        if ylim != None: plt.ylim(ylim)

        # Add axis labels and a legend
        plt.xlabel("Total number of threads $t$", fontsize='large')
        plt.ylabel(fullnames[self._phase] + " time $T$ (s)", fontsize='large')
        plt.legend(title="Modes") if self._system else plt.legend(title="Systems")

        # Save the figure
        systemidentifier = self._system + "_" if self._system else ""
        filename = "scaling_" + self._phase + "_" + systemidentifier + "times.pdf"
        filepath = os.path.join(self._vispath, filename)
        plt.savefig(filepath, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    # -----------------------------------------------------------------

    ## This function creates a PDF plot showing the speedup in function of thread count.
    #  The speedup is defined as T(1)/T(N). It is a dimensionless quantity, theoretically it is >= 1.
    #  The function takes the following (optional) arguments:
    #
    #  - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
    #  - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
    #  - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
    #
    def plotspeedups(self, figsize=(12,8), xlim=None, ylim=None, plotfit=False):

        # Inform the user of the fact that the speedups are being calculated and plotted
        self._log.info("Calculating and plotting the speedups...")

        # Initialize a figure with the appropriate size
        plt.figure(figsize=figsize)

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Create a dictionary that stores the fitted parameters for each system and mode
        parameters = dict.fromkeys(self._statistics.keys())

        # Loop over the different combinations of systemname and mode (i.e. different curves)
        for (systemname, mode), [nthreads, times, errors] in self._statistics.items():

            # Determine a label to identify this curve
            if "hybrid" in mode: mode = "hybrid " + mode.strip('hybrid') + ":1"
            label = mode if self._system else systemname + " (" + mode + ")"

            # Get the serial runtime and its standard deviation
            serialtime = self._serialruntime[systemname]
            serialerror = self._serialerror[systemname]

            # Calculate the speedups
            speedups = serialtime / times

            # Calculate the errors on the efficiencies
            speedup_errors = [ speedups[i] * math.sqrt( math.pow(serialerror/serialtime, 2) + math.pow(errors[i]/times[i], 2) ) for i in range(len(nthreads)) ]

            # Plot the data points for this curve
            plt.errorbar(nthreads, speedups, speedup_errors, marker='.', label=label)

            # Add the appropriate ticks
            ticks |= set(nthreads)

            # Set the weights of the different speedup points for the fitting procedure
            speedup_weigths = speedup_errors if not np.any(np.isinf(speedup_errors)) else None

            # Fit (standard or modified) Amdahl's law to the speedups
            if len(nthreads) < 10:

                popt, pcov = curve_fit(Amdahl, nthreads, speedups, sigma=speedup_weigths, absolute_sigma=False)
                perr = np.sqrt(np.diag(pcov))
                parameters[(systemname, mode)] = [popt[0], perr[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

            else:

                popt, pcov = curve_fit(modAmdahl, nthreads, speedups, sigma=speedup_weigths, absolute_sigma=False)
                perr = np.sqrt(np.diag(pcov))
                parameters[(systemname, mode)] = [popt[0], perr[0], popt[1], perr[1], popt[2], perr[2], popt[3], perr[3]]

        # Use a logarithmic scale for both axes
        plt.xscale('log')
        plt.yscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1]*2)

        # Open the file that will contain the fitted parameters
        parameterfilepath = os.path.join(self._vispath, "parameters_" + self._phase + ".txt")
        parameterfile = open(parameterfilepath, 'w')
        parameterfile.write("# Fit parameters for the speedups to Amdahl's law:\n")
        parameterfile.write("#     S_n = 1 / ( 1 - p + p/n + a + b*n + c*n^2 ) \n")
        parameterfile.write("# Column 1: System name\n")
        parameterfile.write("# Column 2: Mode (mpi, threads or hybrid)\n")
        parameterfile.write("# Column 3: Parallel fraction p\n")
        parameterfile.write("# Column 4: Error on p\n")
        parameterfile.write("# Column 5: Parameter a\n")
        parameterfile.write("# Column 6: Error on a\n")
        parameterfile.write("# Column 7: Parameter b\n")
        parameterfile.write("# Column 8: Error on b\n")
        parameterfile.write("# Column 9: Parameter c\n")
        parameterfile.write("# Column 10: Error on c\n")

        # Plot the fitted speedup curves and write the parameters to the file
        fit_nthreads = np.logspace(np.log10(ticks[0]), np.log10(ticks[-1]), 50)
        for (systemname, mode), [p, p_err, a, a_err, b, b_err, c, c_err] in parameters.items():

            # If requested, plot the fitted curve
            if plotfit:

                # Calculate the fitted speedups
                fit_speedups = [modAmdahl(n, p, a, b, c) for n in fit_nthreads]

                # Add the plot
                plt.plot(fit_nthreads, fit_speedups)

            # Write the parameters to file
            parameterfile.write(systemname + " " + mode + " " + str(p) + " " + str(p_err) + " " + str(a) + " " + str(a_err)
                                                        + " " + str(b) + " " + str(b_err) + " " + str(c) + " " + str(c_err)
                                                        + "\n")

        # Close the parameter file
        parameterfile.close()

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        plt.ylim(ticks[0], ticks[-1])
        plt.grid(True)

        # Plot a line that denotes linear scaling (speedup = nthreads)
        plt.plot(ticks, ticks, linestyle='--')

        # Set custom axis limits if requested
        if xlim != None: plt.xlim(xlim)
        if ylim != None: plt.ylim(ylim)

        # Add axis labels and a legend
        plt.xlabel("Total number of threads $t$", fontsize='large')
        plt.ylabel("Speedup $S_t$", fontsize='large')
        plt.legend(title="Modes") if self._system else plt.legend(title="Systems")

        # Save the figure
        systemidentifier = self._system + "_" if self._system else ""
        filename = "scaling_" + self._phase + "_" + systemidentifier + "speedups.pdf"
        filepath = os.path.join(self._vispath, filename)
        plt.savefig(filepath, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    ## This function creates a PDF plot showing the efficiency in function of thread count.
    #  Efficiency is defined as T(1)/T(N)/N. It is a dimensionless quantity <= 1.
    #  The function takes the following (optional) arguments:
    #
    #  - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
    #  - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
    #  - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
    #
    def ploteffs(self, figsize=(12,8), xlim=None, ylim=None):

        # Inform the user of the fact that the efficiencies are being calculated and plotted
        self._log.info("Calculating and plotting the efficiencies...")

        # Initialize a figure with the appropriate size
        plt.figure(figsize=figsize)

        # Create a set that stores the tick labels for the plot
        ticks = set()

        # Loop over the different combinations of systemname and mode (i.e. different curves)
        for (systemname, mode), [nthreads, times, errors] in self._statistics.items():

            # Determine a label to identify this curve
            if "hybrid" in mode: mode = "hybrid " + mode.strip('hybrid') + ":1"
            label = mode if self._system else systemname + " (" + mode + ")"

            # Get the serial runtime and its standard deviation
            serialtime = self._serialruntime[systemname]
            serialerror = self._serialerror[systemname]

            # Calculate the efficiencies
            efficiencies = (serialtime / times) / nthreads

            # Calculate the errors on the efficiencies
            eff_errors = [ efficiencies[i] * math.sqrt( math.pow(serialerror/serialtime, 2) + math.pow(errors[i]/times[i], 2) ) for i in range(len(nthreads)) ]

            # Plot the data points for this curve
            plt.errorbar(nthreads, efficiencies, eff_errors, marker='.', label=label)

            # Add the appropriate ticks
            ticks |= set(nthreads)

        # Use a logaritmic scale for the x axis (nthreads)
        plt.xscale('log')

        # Add one more tick for esthetic reasons
        ticks = sorted(ticks)
        ticks.append(ticks[-1]*2)

        # Format the axis ticks and create a grid
        ax = plt.gca()
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlim(ticks[0], ticks[-1])
        plt.ylim(0, 1.1)
        plt.grid(True)

        # Set custom axis limits if requested
        if xlim != None: plt.xlim(xlim)
        if ylim != None: plt.ylim(ylim)

        # Add axis labels and a legend
        plt.xlabel("Total number of threads $t$", fontsize='large')
        plt.ylabel("Efficiency $\epsilon_t$", fontsize='large')
        plt.legend(title="Modes") if self._system else plt.legend(title="Systems")

        # Save the figure
        systemidentifier = self._system + "_" if self._system else ""
        filename = "scaling_" + self._phase + "_" + systemidentifier + "efficiencies.pdf"
        filepath = os.path.join(self._vispath, filename)
        plt.savefig(filepath, bbox_inches='tight', pad_inches=0.25)
        plt.close()

# -----------------------------------------------------------------

## This function defines Amdahl's law for the speedup
def Amdahl(n, p):

    return 1.0/(1 - p + p/n)

# This function defines a modified version of Amdahl's law, which accounts for different kinds of overhead
def modAmdahl(n, p, a, b, c):

    return 1.0/(1 - p + p/n + a + b*n + c*n**2)

# -----------------------------------------------------------------