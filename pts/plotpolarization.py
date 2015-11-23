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

import numpy as np
import pyfits

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

from matplotlib import ticker
import pts.archive as arch

# -----------------------------------------------------------------

# This function creates a polarization map for each set of "prefix_instr_stokes*.fits" files in the output of
# the specified simulation. The plots are saved in a PDF file named "prefix_instr_stokes.pdf",
# placed next to the original set of files.
#
# The function takes the following arguments:
# - simulation: the SkirtSimulation object representing the simulation to be handled
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 6x6 inch
# - binsize: the number of pixels in each bin, in horizontal and vertical directions
# - wavelength: the wavelength for which to create the plot, in micron; if missing the first frame is used
def plotpolarization(simulation, figsize=(6,6), binsize=(10,10), wavelength=None):
    binX = binsize[0]
    binY = binsize[1]
    for (fileI, fileQ, fileU, fileV), name in zip(simulation.stokesfitspaths(), simulation.instrumentnames()):
        print "Creating PDF polarization map for instrument {} of simulation {}".format(name, simulation.prefix())

        # setup the figure
        figure = plt.figure()

        # determine the appropriate frame index
        index = 0 if wavelength==None else simulation.frameindices([wavelength])[0]

        # read the Stokes frames with shape (ny, nx) or datacubes with shape (nlambda, ny, nx)
        # if datacubes, pick the first frame (for now)
        I = pyfits.getdata(arch.openbinary(fileI))
        Q = pyfits.getdata(arch.openbinary(fileQ))
        U = pyfits.getdata(arch.openbinary(fileU))
        V = pyfits.getdata(arch.openbinary(fileV))
        if len(I.shape)==3:
            I = I[index,:,:]
            Q = Q[index,:,:]
            U = U[index,:,:]
            V = V[index,:,:]

        # determine number of pixels dropped by binning and inform user
        orLenX = np.shape(I)[1]
        dropX =  orLenX % binX
        startX = dropX/2
        orLenY = np.shape(I)[0]
        dropY = orLenY % binY
        startY = dropY/2
        dropped = dropX * orLenY + dropY * orLenX - dropX * dropY
        droppedFraction = dropped / orLenX / orLenY
        print "Picture dimensions: %s; bin size: %s by %s; dropping %s pixels (%.2f%%) (l:%s, r:%s, b:%s, t:%s)" \
            % (np.shape(I.T), binX, binY, dropped, 100*droppedFraction, startX, dropX - startX, startY, dropY-startY)

        # actual binning
        posX = np.arange(startX-0.5+binX/2.0, orLenX - dropX + startX - 0.5, binX)
        posY = np.arange(startY-0.5+binY/2.0, orLenY - dropY + startY - 0.5, binY)
        binnedI = np.zeros((len(posY),len(posX)))
        binnedQ = np.zeros((len(posY),len(posX)))
        binnedU = np.zeros((len(posY),len(posX)))
        binnedV = np.zeros((len(posY),len(posX)))
        for x in range(len(posX)):
            for y in range(len(posY)):
                binnedI[y,x] = np.sum(I[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                binnedQ[y,x] = np.sum(Q[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                binnedU[y,x] = np.sum(U[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])
                binnedV[y,x] = np.sum(V[startY+binY*y : startY+binY*(y+1) , startX+binX*x : startX+binX*(x+1)])

        # degree of polarization with low resolution for foreground segments
        # compute the linear polarization degree
        degreeLD = np.sqrt(binnedQ**2 + binnedU**2)
        degreeLD[degreeLD>0] /= binnedI[degreeLD>0]

        # compute the polarization angle
        angle = 0.5 * np.arctan2(binnedU, binnedQ)

        # plot contour of intensity (in HD)
        I[I<=0] = np.NaN
        backGrndPlt = plt.contourf(I, alpha=0.4, locator=ticker.LogLocator())
        cbar = plt.colorbar(backGrndPlt)
        cbar.set_label(simulation.fluxlabel(), labelpad=5)

        # determine a characteristic 'high' degree of polarization in the Frame (For key and scaling)
        charDegree = np.percentile(degreeLD, 99.0) # This has to be done before degreeLD contains 'np.NaN'
        # removing pixels with zero polarization
        degreeLD[degreeLD<=0] = np.NaN   # a runtime error here means that a laser is pointed at your instrument
        if not 0<charDegree<1:
            charDegree = max(np.max(degreeLD), 0.01)

        # calculate the scaling. It is relative to the picture size, thus the high polarization 'charDegree' should
        # be short enough so it does not overlap with the polarization around it.
        scale = (charDegree* 2.2)*max(len(posX)/figsize[0], len(posY)/figsize[1])
        keyLength = 10**(int(np.log10(charDegree))) # a runtime error here means no polarization in the whole picture
        key = "{:.3g}%".format(100 * keyLength)

        # create the polarization vector arrays
        cosPolarization = degreeLD*np.cos(angle)
        sinPolarization = degreeLD*np.sin(angle)

        # plot the vector field (in LD)
        X,Y = np.meshgrid( posX , posY)
        quiverplot = plt.quiver(X,Y, cosPolarization, sinPolarization, pivot='middle', units='inches',
                                angles = 'xy', scale = scale, scale_units = 'inches', headwidth=0, headlength=1,
                                headaxislength=1, minlength = 0.6)

        # add legend, labels and title
        plt.quiverkey(quiverplot, 0.8, 0.02, keyLength, key,
                      coordinates='axes', labelpos='E')
        figure.suptitle(simulation.prefix(), fontsize=14, fontweight='bold')
        plt.title('instrument: ' + name + ', $\lambda=' +
                    str(simulation.wavelengths()[index])+'\ \mu m$', fontsize=12)
        plt.axis('scaled')
        plt.xlabel('x (pixels)')
        plt.ylabel('y (pixels)')

        # save and close the figure
        plotfile = fileQ.replace("Q.fits",".pdf")
        plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        print "Created PDF polarization map " + plotfile

# -----------------------------------------------------------------
