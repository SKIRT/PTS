#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.seds Plot SED information contained in text files
#
# The function in this module creates a PDF plot for one or more SEDs, each provided as a sequence
# of SED data points in a text file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ..tools import archive as arch

# -----------------------------------------------------------------

# This function creates a plot combining the "sed.dat" files in the output of the specified simulation.
# The plots are saved in PDF format and are placed next to the original file(s) with
# a similar name (omitting the instrument name) but a different extension.
# The function takes the following arguments:
# - simulation: the SkirtSimulation object representing the simulation to be handled
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
# - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
# - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
def plotseds(simulation, figsize=(10,6), xlim=None, ylim=None, output_path=None, format="PDF"):
    sedpaths = simulation.seddatpaths()
    if len(sedpaths) > 0:
        labels = [ path.rsplit("_",2)[1] for path in sedpaths ]
        outpath = sedpaths[0].rsplit("_",2)[0] + "_sed." + format.lower()
        if output_path is not None: outpath = os.path.join(output_path, os.path.basename(outpath))
        success = plotseds_impl(sedpaths, outpath, labels, simulation.fluxlabel(), figsize=figsize, xlim=xlim, ylim=ylim)
        if success: print("Created " + format + " SED plot file " + outpath)

# -----------------------------------------------------------------

## This function creates a PDF plot for one or more SEDs, all plotted on the same log/log axes.
# Each SED is provided as a separate text file. Each line in such a file contains two columns,
# respectively specifying the wavelength and the corresponding flux.
# The current implementation assumes that the wavelength is specified in micron.
#
# The function takes the following arguments:
# - sedfiles: a list of filepaths for the input files containing the SED data points
# - plotfile: the filepath of the plot file to be generated; this \em must have the \c .pdf filename extension
# - labels: a list of labels, one for each sed file, in the same order, for use in the plot's legend;
#   if this argument is missing the sed filenames are used as labels
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
# - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
# - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
#
def plotseds_impl(sedfiles, plotfile, labels=None, fluxlabel="Flux", figsize=(10,6), xlim=None, ylim=None):

    # Initialize figure with the appropriate size
    plt.figure(figsize=figsize)
    plt.clf()

    if labels is None: labels = sedfiles

    # setup the figure
    plt.grid(True)

    # loop over sed files and labels
    for sedfile,label in zip(sedfiles,labels):
        data = np.loadtxt(arch.opentext(sedfile))
        if len(data.shape)==2:
            plt.loglog(data[:,0], data[:,1], label=label)
        else:
            return False  # there is just one data point

    # set axis limits if requested
    if xlim != None: plt.xlim(xlim)
    if ylim != None: plt.ylim(ylim)

    # add axis labels and a legend
    plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    plt.xlabel(r"$\lambda\,(\mu \mathrm{m})$", fontsize='large')
    plt.ylabel(fluxlabel, fontsize='large')
    plt.legend()

    # Save the figure
    plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)
    plt.close()

    return True

# -----------------------------------------------------------------
