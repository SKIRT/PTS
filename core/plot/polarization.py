#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plotpolarization Plot polarization information contained in \c stokes*.fits files
#
# The function in this module creates a PDF polarization map for a set of "prefix_instr_stokes*.fits" files
# produced by a SKIRT simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import pyfits
import matplotlib
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ..tools import archive as arch

# -----------------------------------------------------------------

# This function creates a polarization map for each set of "prefix_instr_stokes*.fits" files in the output of
# the specified simulation. The plots are saved in a PDF file named "prefix_instr_stokes.pdf",
# placed next to the original set of files.
#
# The function takes the following arguments:
# - simulation: the SkirtSimulation object representing the simulation to be handled
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
def plotpolarization(simulation, figsize=(10,6)):

    for fileI, fileQ, fileU, fileV in simulation.stokesfitspaths():

        # setup the figure
        figure = plt.figure(figsize=figsize)

        # read the Stokes frames with shape (ny, nx) or datacubes with shape (nlambda, ny, nx)
        # if datacubes, pick the first frame (for now)
        I = pyfits.getdata(arch.openbinary(fileI))
        Q = pyfits.getdata(arch.openbinary(fileQ))
        U = pyfits.getdata(arch.openbinary(fileU))
        V = pyfits.getdata(arch.openbinary(fileV))
        if len(I.shape)==3:
            index = 0
            I = I[index,:,:]
            Q = Q[index,:,:]
            U = U[index,:,:]
            V = V[index,:,:]

        # compute the (linear or total) polarization degree
        #degree = np.sqrt(Q**2 + U**2)
        degree = np.sqrt(Q**2 + U**2 + V**2)
        degree[degree>0] /= I[degree>0]
        degree[degree<=0] = np.NaN   # so that areas with zero polarization are left blank

        # compute the polarization angle
        angle = 0.5 * np.arctan2(U,Q)
        # add an offset according to Calamai, Landi Degl'Innocenti & Landi Degl'Innocenti (1975)
        #angle[ Q<0 ] += np.pi/2
        #angle[ (Q>0) & (U<0) ] += np.pi
        #angle[ (Q==0) & (U>0) ] = np.pi/4
        #angle[ (Q==0) & (U<0) ] = 3*np.pi/4

        # determine the scale (fixed for now)
        scale = 0.001

        # plot the vector field
        quiverplot = plt.quiver(degree*np.cos(angle), degree*np.sin(angle),
                                pivot='middle', headwidth=0, headlength=1, headaxislength=1,
                                units='xy', scale_units='xy', scale=scale)

        # add a legend
        plt.quiverkey(quiverplot, 0.05, 0.05, scale, "{} %".format(scale*100),
                      coordinates='axes', labelpos='E')

        # save and close the figure
        plotfile = fileQ.replace("Q.fits",".pdf")
        plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        print("Created PDF polarization map" + plotfile)

# -----------------------------------------------------------------
