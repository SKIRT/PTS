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
import matplotlib
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..basics.log import Log

# -----------------------------------------------------------------

class ScalingPlotter(object):

    """
    An instance of the ScalingPlotter is used to create plots of the runtimes, speedups and efficiencies as a function
    of the number of threads, based on the results of (multiple) SKIRT scaling benchmark tests.
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Set the table to None initially
        self.table = None

        # Set the output path to None initially
        self.output_path = None

        # Create a logger
        self.log = Log()

    # -----------------------------------------------------------------

    def run(self, input, output_path):

        """
        This function ...
        :return:
        """

        # If the input is a Table object
        if isinstance(input, Table): self.table = input

        # If the input is a string
        elif isinstance(input, basestring): self.table = Table.read(input, format="ascii")

        # Invalid input
        else: raise ValueError("Input must be either an Astropy Table object or a filename (e.g. memory.dat)")

        # Set the path to the output directory
        self.output_path = output_path

        # Make the plots
        self.plot()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Plot the scaling relation of the runtime
        self.plot_times()

        # Plot the scaling relation of the speedup
        self.plot_speedups()

        # Plot the scaling relation of the efficiency
        self.plot_efficiencies()

        # Plot the scaling relation of the memory usage
        self.plot_memory()

    # -----------------------------------------------------------------

    def plot_times(self, figsize=(12,8), xlim=None, ylim=None):

        """
        This function creates a PDF plot showing the execution time as a function of the number of threads.
        It takes the following arguments:
        :param figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
        :param xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
        :param ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
        :return:
        """

        # Inform the user of the fact that the runtimes are being plotted
        self.log.info("Plotting the runtimes...")

        # Determine the file path for this plot
        file_path = os.path.join(self.output_path, "times.pdf")

        # Create a PDF Pages object
        pp = PdfPages(file_path)

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

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
        #systemidentifier = self._system + "_" if self._system else ""
        #filename = "scaling_" + self._phase + "_" + systemidentifier + "times.pdf"

        # Save the figure
        pp.savefig()
        pp.close()

    # -----------------------------------------------------------------

    def plot_speedups(self, figsize=(12,8), xlim=None, ylim=None, plotfit=False):

        """
        This function creates a PDF plot showing the speedup as a function of the number of threads.
        The speedup is defined as T(1)/T(N). It is a dimensionless quantity, theoretically it is >= 1.
        This function takes the following (optional) arguments:
        :param figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
        :param xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
        :param ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
        :param plotfit:
        :return:
        """

        # Inform the user of the fact that the speedups are being calculated and plotted
        self.log.info("Calculating and plotting the speedups...")

        # Determine the file path for this plot
        file_path = os.path.join(self.output_path, "speedups.pdf")

        # Create a PDF Pages object
        pp = PdfPages(file_path)

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

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
        #systemidentifier = self._system + "_" if self._system else ""
        #filename = "scaling_" + self._phase + "_" + systemidentifier + "speedups.pdf"

        # Save the figure
        pp.savefig()
        pp.close()

    # -----------------------------------------------------------------

    def plot_efficiencies(self, figsize=(12,8), xlim=None, ylim=None):

        """
        This function creates a PDF plot showing the efficiency as a function of the number of threads.
        Efficiency is defined as T(1)/T(N)/N. It is a dimensionless quantity <= 1.
        The function takes the following (optional) arguments:
        :param figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
        :param xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
        :param ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
        :return:
        """

        # Inform the user of the fact that the efficiencies are being calculated and plotted
        self.log.info("Calculating and plotting the efficiencies...")

        # Determine the file path for this plot
        file_path = os.path.join(self.output_path, "efficiencies.pdf")

        # Create a PDF Pages object
        pp = PdfPages(file_path)

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

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
        #systemidentifier = self._system + "_" if self._system else ""
        #filename = "scaling_" + self._phase + "_" + systemidentifier + "efficiencies.pdf"

        # Save the figure
        pp.savefig()
        pp.close()

    # -----------------------------------------------------------------

    def plot_memory(self, figsize=(12,8)):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Plotting the memory scaling...")

        # Determine the file path for this plot
        file_path = os.path.join(self.output_path, "times.pdf")

        # Create a PDF Pages object
        pp = PdfPages(file_path)

        # Initialize figure with the appropriate size
        plt.figure(figsize=figsize)
        plt.clf()

        # ...

        # Save the figure
        pp.savefig()
        pp.close()

# -----------------------------------------------------------------

## This function defines Amdahl's law for the speedup
def Amdahl(n, p): return 1.0/(1 - p + p/n)

# This function defines a modified version of Amdahl's law, which accounts for different kinds of overhead
def modAmdahl(n, p, a, b, c): return 1.0/(1 - p + p/n + a + b*n + c*n**2)

# -----------------------------------------------------------------