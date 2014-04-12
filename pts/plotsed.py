#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plotsed Plot SED information contained in text files
#
# The function in this module creates a PDF plot for one or more SEDs, each provided as a sequence
# of SED data points in a text file.

# -----------------------------------------------------------------

import numpy as np
import os.path

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------

## This function creates a PDF plot for one or more SEDs, all plotted on the same log/log axes.
# Each SED is provided as a separate text file. Each line in such a file contains two columns,
# respectively specifying the wavelength \f$\lambda\f$ and the corresponding flux
# \f$\lambda\,F_\lambda\f$.
# The current implementation assumes that the wavelength is specified in micron,
# and that the flux is specified in \f$\text{W}\,\text{m}^{-2}\f$.
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
def plotsed(sedfiles, plotfile, labels=None, figsize=(10,6), xlim=None, ylim=None):
    plotfile = os.path.expanduser(plotfile)
    assert plotfile.endswith(".pdf")
    if labels == None: labels = sedfiles

    # setup the figure
    figure = plt.figure(figsize=figsize)

    # loop over sed files and labels
    for sedfile,label in zip(sedfiles,labels):
        data = np.loadtxt(sedfile)
        if len(data.shape)==2:
            plt.loglog(data[:,0], data[:,1], label=label)
        else:
            return False  # there is just one data point

    # set axis limits if requested
    if xlim != None: plt.xlim(xlim)
    if ylim != None: plt.ylim(ylim)

    # add axis labels and a legend
    plt.xlabel("$\lambda\,(\mu m)$", fontsize='large')
    plt.ylabel("$\lambda\,F_\lambda\,(\mathrm{W}\,\mathrm{m}^{-2})$", fontsize='large')
    plt.legend()

    # save the figure
    plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)

    # close things
    plt.close()
    return True

# -----------------------------------------------------------------
