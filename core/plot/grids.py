#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.grids Plot a dust grid wireframe to PDF
#
# The function in this module creates a PDF plot from a sequence of grid coordinates
# provided as an input text file, directly creating a PDF file (i.e. without using matplotlib).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np
from reportlab.pdfgen import canvas

# Import the relevant PTS classes and modules
from ..tools import archive as arch

# -----------------------------------------------------------------

# This function creates a plot for each "gridxx.dat" file in the output of the specified simulation.
# The plots are saved in PDF format and are placed next to the original file(s) with
# the same name but a different extension.
#
# There are two format variations for 2D and 3D information, respectively. The 2D format describes the
# intersection of a dust grid structure with one of the coordinate planes. The 3D format fully describes
# all or part of the dust cells in the grid. Each line in the file contains two (2D) or three (3D)
# coordinates seperated by whitespace, or is empty. Consecutive nonempty lines represent a sequence of
# "lineto" commands; an empty line marks a "moveto" command.
#
# The function takes the following arguments:
# - simulation: the SkirtSimulation object representing the simulation to be handled
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 8 x 8 inch
#
def plotgrids(simulation, figsize=(8,8), output_path=None, silent=False, prefix=None, linewidth=0.1, maxlevel=None):

    # Set file prefix
    if prefix is None: prefix = ""
    else: prefix = prefix + "_"

    # If max level is defined, get log file
    if maxlevel is not None:
        logfile = simulation.log_file
        cell_tree_distribution = logfile.tree_leaf_distribution
        ncells = cell_tree_distribution.get_ncells_below_level(maxlevel, including=True)
    else: ncells = None

    # Loop over the grid files
    for gridfile in simulation.gridxxdatpaths():

        # Determine output file path
        plotfile = gridfile[:-4] + ".pdf"
        if output_path is not None: plotfile = os.path.join(output_path, prefix + os.path.basename(plotfile))

        # Make the plot
        make_grid_plot(gridfile, plotfile, figsize=figsize, silent=silent, linewidth=linewidth, ncells=ncells)

# -----------------------------------------------------------------

def create_figure(figsize, plotfile, linewidth=0.1):

    # setup the figure with the appropriate size (in points)
    figwidth = 72 * figsize[0]
    figheight = 72 * figsize[1]
    if figwidth == figheight: figheight += 2  # to ensure portrait orientation when printed
    fig = canvas.Canvas(plotfile, pagesize=(figwidth, figheight))
    fig.setAuthor("Python Toolkit for SKIRT")
    fig.setLineWidth(linewidth)

    # Return
    return fig, figwidth, figheight

# -----------------------------------------------------------------

def make_grid_plot2d(gridfile, plotfile, figsize=(8,8), silent=False, linewidth=0.1, ncells=None):

    # Create
    fig, figwidth, figheight = create_figure(figsize, plotfile, linewidth=linewidth)

    # determine the extent of the grid being plotted
    xmin, ymin, xmax, ymax = float('Inf'), float('Inf'), float('-Inf'), float('-Inf')
    for line in arch.opentext(gridfile):
        coords = line.split()
        if len(coords) == 2:
            x, y = float(coords[0]), float(coords[1])
            xmin, ymin, xmax, ymax = min(xmin, x), min(ymin, y), max(xmax, x), max(ymax, y)

    # determine the scales and offsets to fit and center the grid in the figure
    xs = figwidth * 0.95 / (xmax - xmin)
    ys = figheight * 0.95 / (ymax - ymin)
    xo = (figwidth - xs * (xmax - xmin)) / 2. - xmin * xs
    yo = (figheight - ys * (ymax - ymin)) / 2. - ymin * ys

    # for each set of consecutive nonempty lines in the file, draw the line segments
    path = None
    current_ncells = 0
    for line in arch.opentext(gridfile):
        coords = line.split()
        if len(coords) == 0 and path != None:
            fig.drawPath(path)
            path = None
        if len(coords) == 2:
            if ncells is not None and current_ncells > ncells: break
            x, y = xo + xs * float(coords[0]), yo + ys * float(coords[1])
            if path is None:
                path = fig.beginPath()
                path.moveTo(x, y)
            else:
                path.lineTo(x, y)
            current_ncells += 1

    # save the figure
    fig.showPage()
    fig.save()
    if not silent: print("Created PDF grid plot file " + plotfile)

# -----------------------------------------------------------------

def make_grid_plot3d(gridfile, plotfile, figsize=(8,8), silent=False, linewidth=0.1, ncells=None):

    # Create
    fig, figwidth, figheight = create_figure(figsize, plotfile, linewidth=linewidth)

    # determine the extent of the grid being plotted (largest half-width in all directions)
    extent = 0.
    for line in arch.opentext(gridfile):
        coords = line.split()
        if len(coords) == 3:
            extent = max(extent, np.max(np.abs(np.array(map(float, coords)))))

    # determine the scale and offsets to fit and center the grid in the figure
    s = min(figwidth, figheight) * 0.95 / (2 * extent) / np.sqrt(3)
    xo = extent * s + (figwidth - s * (2 * extent)) / 2.
    yo = extent * s + (figheight - s * (2 * extent)) / 2.

    # set the viewing angles
    inclination, azimuth = 65, 40
    costheta = np.cos(inclination * np.pi / 180.)
    sintheta = np.sin(inclination * np.pi / 180.)
    cosphi = np.cos(azimuth * np.pi / 180.)
    sinphi = np.sin(azimuth * np.pi / 180.)

    # for each set of consecutive nonempty lines in the file, draw the line segments
    path = None
    current_ncells = 0
    for line in arch.opentext(gridfile):
        coords = line.split()
        if len(coords) == 0 and path != None:
            fig.drawPath(path)
            path = None
        if len(coords) == 3:
            if ncells is not None and current_ncells > ncells: break
            x, y, z = map(float, coords)
            xp = - sinphi * x + cosphi * y
            yp = - cosphi * costheta * x - sinphi * costheta * y + sintheta * z
            xf = xo + s * xp
            yf = yo + s * yp
            if path is None:
                path = fig.beginPath()
                path.moveTo(xf, yf)
            else: path.lineTo(xf, yf)
            current_ncells += 1

    # save the figure
    fig.showPage()
    fig.save()
    if not silent: print("Created PDF grid plot file " + plotfile)

# -----------------------------------------------------------------

def make_grid_plot(gridfile, plotfile, figsize=(8,8), silent=False, linewidth=0.1, ncells=None):

    # determine the format type from the first nonempty line (3D format has 3 columns, 2D format has 2 columns)
    for line in arch.opentext(gridfile):
        form = len(line.split())
        if form > 0: break

    # 2D
    if form == 2: make_grid_plot2d(gridfile, plotfile, figsize=figsize, silent=silent, linewidth=linewidth, ncells=ncells)

    # 3D
    else: make_grid_plot3d(gridfile, plotfile, figsize=figsize, silent=silent, linewidth=linewidth, ncells=ncells)

# -----------------------------------------------------------------
