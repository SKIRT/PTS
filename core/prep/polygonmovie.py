#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.polygonmovie Creating a movie from a sequence of polygon coordinates.
#
# The function in this module creates a movie from a sequence of polygon coordinates
# provided as an input file.

# -----------------------------------------------------------------

# Import standard modules
import warnings
import numpy as np
import os.path

# Import the relevant PTS classes and modules
from ..basics.log import log

# Use a non-interactive back-end to generate high-quality raster graphics
import matplotlib
with warnings.catch_warnings():
    warnings.filterwarnings('error')
    try:
        if matplotlib.get_backend().lower() != "agg": matplotlib.use("agg")
    except Warning as w: log.warning("An failed attempt of setting the Matplotlib backend has been made because it has already been set")
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ..basics.moviefile import MovieFile

# -----------------------------------------------------------------

## This function creates a movie from a sequence of polygon coordinates provided as an input file.
# The intention is to see a polygon "crawl around the canvas".
#
# Each line in the input file specifies the 2*N coordinates of a planar polygon, as space-separated
# floating point numbers in the order X1 Y1 X2 Y2 ... XN YN. All lines in the file must have the
# same number of coordinates. For each line the output movie will have a frame containing a line
# drawing of the polygon represented by the coordinates on that line.
#
# The function takes the following arguments:
# - coordfile: the filepath of the input file containing the polygon coordinates
# - moviefile: the filepath of the movie to be generated; this \em must have the .mov or .mpg filename extension
# - background: a function with two arguments x,y that returns the intensity of the background, or None (the default)
# - xrange, yrange: the coordinate range in the x and y directions; the default value is (0,1) for both
# - shape: the number of movie frame pixels in the x and y directions; the default value is (500,500)
# - rate: the number of frames per second in the output movie; the default value is 24 fps
#
def createpolygonmovie(coordfile, moviefile, background=None, xrange=(0,1), yrange=(0,1), shape=(500,500), rate=24):
    # read the coordinates
    polygons = np.loadtxt(os.path.expanduser(coordfile))

    # open the movie file
    movie = MovieFile(moviefile, shape=shape, rate=rate)

    # setup the figure and its background
    figure = plt.figure(1, dpi=100, figsize=(shape[0]/100.,shape[1]/100.))
    x0,x1 = xrange
    y0,y1 = yrange
    px,py = shape
    if background != None:
        # transpose and scale i,j indices to x,y values
        backim = np.fromfunction(lambda i,j: background((j+0.5)/px*(x1-x0)+x0,(i+0.5)/py*(y1-y0)+y0), shape[::-1])
        plt.imshow(backim, origin='lower', aspect='auto', interpolation='bicubic', extent=(x0,x1,y0,y1))
    plt.grid(True)
    plt.xlim(xrange[0], xrange[1])
    plt.ylim(yrange[0], yrange[1])

    # for each line in the file, draw the polygon and add a movie frame
    for p in polygons:
        plt.plot(np.concatenate((p[0::2],p[0:1])), np.concatenate((p[1::2],p[1:2])), 'w-')
        plt.draw()                                  # flush the drawing
        movie.add(figure.canvas.buffer_rgba(0,0))   # pass the pixel buffer on as movie frame
        figure.axes[0].lines.pop()                  # remove the plotted line

    # close things
    plt.close()
    movie.close()

# -----------------------------------------------------------------
