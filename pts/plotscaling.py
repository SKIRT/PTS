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

# Use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------

## This function returns a nice human readable label for a benchmark result file name,
# which is assumed be formatted as <system-name>_<#hard-cores>_<#extra-cores>.dat
def label(filename):
    segments = filename.split("_")
    system = segments[0]
    mode = segments[1]
    threads = segments[2].split(".")[0]
    return system + " (" + mode + ")"

# -----------------------------------------------------------------

## This function combines time samples for the same thread count into a single averaged entry with an error bar
def combinesamples(counts,times):
    result = { }        # key: count; value: [ number of entries so far, average time, minimum, maximum ]
    for count,time in zip(counts,times):
        if not count in result: result[count] = [ 0, 0., float('Inf'), 0. ]
        n = result[count][0]
        result[count][0] = n+1
        result[count][1] = (result[count][1]*n + time) / float(n+1)
        result[count][2] = min(result[count][2], time)
        result[count][3] = max(result[count][3], time)

    return (  np.array([count for count in sorted(result.keys())]),
              np.array([result[count][1] for count in sorted(result.keys())]),
              np.array(( [abs(result[count][2]-result[count][1]) for count in sorted(result.keys())],
                         [abs(result[count][3]-result[count][1]) for count in sorted(result.keys())] )) )

# -----------------------------------------------------------------

## This function returns the total number of threads for a benchmark result file name,
# which is assumed be formatted as <system-name>_<mode>_<#processes>_<#threads>.dat
def totalthreads(filename):
    segments = filename.split("_")
    processes = segments[2]
    threads = segments[3].split(".")[0]
    return int(processes)*int(threads)

# -----------------------------------------------------------------

## This function makes the scaling plots
def plotscaling(directory, filename):

    if filename != "":
        filenames = [filename]

    else:
        # Get a list of result files (i.e. all *.dat files in the current directory)
        filenames = sorted(filter(lambda fn: fn.endswith(".dat"), os.listdir(directory)), key=totalthreads)

    # Generate the plots
    plottimes(filenames, "scaling_times.pdf", directory, xlim=(0,40))
    plotspeedups(filenames, "scaling_speedups.pdf", directory, xlim=(0,40))
    ploteffs(filenames, "scaling_effs.pdf", directory, xlim=(0,40))

## This function creates a PDF plot showing the execution time in function of thread count.
# The function takes the following arguments:
# - filenames: a list of names for files containing benchmark results
# - plotfile: the name of the plot file to be generated; this \em must have the \c .pdf filename extension
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
# - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
# - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
#
def plottimes(filenames, plotfile, path, figsize=(12,8), xlim=None, ylim=None):

    # Get the number of logical cores on the system
    logical_cores = multiprocessing.cpu_count()

    assert plotfile.endswith(".pdf")
    figure = plt.figure(figsize=figsize)

    # Loop over files
    for filename in filenames:
        
        # Plot the data points in this file
        filepath = os.path.join(path, filename)
        count,time = np.loadtxt(filepath, usecols=(2,4), unpack=True);
        count,time,error = combinesamples(count,time)
        plt.errorbar(count, time, error, marker='.', label=label(filename))

    # Set axis limits if requested
    plt.grid(True)
    if xlim != None: plt.xlim(xlim)
    if ylim != None: plt.ylim(ylim)

    # Add axis labels and a legend
    plt.xlabel("Total number of threads", fontsize='large')
    plt.ylabel("Stellar emission time (s)", fontsize='large')
    plt.legend(title="Systems")

    # Save the figure
    plotfilepath = os.path.join(path, plotfile)
    plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.25)
    plt.close()

# -----------------------------------------------------------------

## This function creates a PDF plot showing the efficiency in function of thread count.
# Efficiency is defined as T(1)/T(N)/N. It is a dimensionless quantity <= 1.
# The function takes the following arguments:
# - filenames: a list of names for files containing benchmark results
# - plotfile: the name of the plot file to be generated; this \em must have the \c .pdf filename extension
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
# - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
# - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
#
def ploteffs(filenames, plotfile, path, figsize=(12,8), xlim=None, ylim=None):

    # Get the number of logical cores on the system
    logical_cores = multiprocessing.cpu_count()

    assert plotfile.endswith(".pdf")
    figure = plt.figure(figsize=figsize)

    # Loop over files
    for filename in filenames:
        
        # Plot the data points in this file
        filepath = os.path.join(path, filename)
        count,time = np.loadtxt(filepath, usecols=(2,4), unpack=True);
        dummy1,avg,dummy2 = combinesamples(count,time)
        eff = avg[0] / time / count
        count,eff,error = combinesamples(count,eff)
        plt.errorbar(count, eff, error, marker='.', label=label(filename))

    # Set axis limits if requested
    plt.grid(True)
    if xlim != None: plt.xlim(xlim)
    if ylim != None: plt.ylim(ylim)

    # Add axis labels and a legend
    plt.xlabel("Total number of threads", fontsize='large')
    plt.ylabel("Efficiency (dimensionless)", fontsize='large')
    plt.legend(title="Systems")

    # Save the figure
    plotfilepath = os.path.join(path, plotfile)
    plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.25)
    plt.close()

# -----------------------------------------------------------------

def plotspeedups(filenames, plotfile, path, figsize=(12,8), xlim=None, ylim=None):

    # Get the number of logical cores on the system
    logical_cores = multiprocessing.cpu_count()

    assert plotfile.endswith(".pdf")
    figure = plt.figure(figsize=figsize)

    # Loop over files
    for filename in filenames:
        
        # Plot the data points in this file
        filepath = os.path.join(path, filename)
        count,time = np.loadtxt(filepath, usecols=(2,4), unpack=True);
        dummy1,avg,dummy2 = combinesamples(count,time)
        eff = avg[0] / time
        count,eff,error = combinesamples(count,eff)
        plt.errorbar(count, eff, error, marker='.', label=label(filename))

    # Set axis limits if requested
    plt.grid(True)
    if xlim != None: plt.xlim(xlim)
    if ylim != None: plt.ylim(ylim)

    # Add axis labels and a legend
    plt.xlabel("Total number of threads", fontsize='large')
    plt.ylabel("Speedup (dimensionless)", fontsize='large')
    plt.legend(title="Systems")

    # Save the figure
    plotfilepath = os.path.join(path, plotfile)
    plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.25)
    plt.close()

# -----------------------------------------------------------------
