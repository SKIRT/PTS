#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plotgrid Plot dust grid coordinates dumped in a text file
#
# The function in this module creates a PDF plot from a sequence of grid coordinates
# provided as an input text file.

# -----------------------------------------------------------------

import numpy as np
import os.path

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------

## This function creates a PDF plot from a sequence of grid coordinates provided as an input text file.
# Each line in the file contains two coordinates seperated by whitespace or is empty. Consecutive nonempty
# lines represent a sequence of "lineto" commands; an empty line marks a "moveto" command.
#
# The function takes the following arguments:
# - gridfile: the filepath of the input file containing the grid coordinates
# - plotfile: the filepath of the plot file to be generated; this \em must have the \c .pdf filename extension
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 8 x 8 inch
# - axes: if True the figure contains axes and labels; if False (the default) there are no axis and labels
#
def plotgrid(gridfile, plotfile, figsize=(8,8), axes=False):
    plotfile = os.path.expanduser(plotfile)
    assert plotfile.endswith(".pdf")

    # setup the figure
    figure = plt.figure(figsize=figsize)
    if not axes:
        ax1 = plt.axes(frameon=False)
        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)

    # for each set of consecutive nonempty lines in the file, draw the line segments
    x,y = np.zeros((500)), np.zeros((500))  # this is the maximum number of points in a path
    index = 0
    for line in file(os.path.expanduser(gridfile)):
        coords = line.split()
        if len(coords)==0 and index>0:
            plt.plot(x[0:index], y[0:index], 'k-', linewidth=0.1)
            index = 0
        if len(coords)==2:
            x[index],y[index] = float(coords[0]), float(coords[1])
            index+=1

    # save the figure
    plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)

    # close things
    plt.close()

# -----------------------------------------------------------------
