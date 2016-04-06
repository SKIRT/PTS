#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.polarization Plot polarization information contained in \c stokes*.fits files
#
# The function in this module creates a PDF polarization map for a set of "prefix_instr_stokes*.fits" files
# produced by a SKIRT simulation.

# -----------------------------------------------------------------

# Import standard modules
import numpy as np
import pyfits
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import os
from matplotlib import ticker
from matplotlib.colors import LogNorm

# Import the relevant PTS classes and modules
from ..tools import archive as arch

# -----------------------------------------------------------------

# This function creates polarization maps of a simulation.
# It places them in a folder "_polarization" next to the input files.
#
# The function takes the following arguments:
# - simulation: the SkirtSimulation object representing the simulation to be handled
# - instrumentList: NOT YET IMPLEMENTED the instruments of the simulation that are to be plotted
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 6x6 inch
# - binsize: the number of pixels in each bin, in horizontal and vertical directions
# - wavelength: the wavelength for which to create the plot, in micron; if missing the first frame is used
# - polAvY: to plot the over the y-axis averaged polarization degree for all x pixels
# - noMaps: to not plot the polarization map. Only additional options. Speeds up the procedure.
# - export: to write things as .dat files in addition to the .pdf files for all enabled additional options.
def plotpolarization(simulation, instrumentList = 'all', figsize = (6,6), binsize = (7,7), wavelength = "all",
                    polAvY=False, noMaps = False, export = True):
    binX = binsize[0]
    binY = binsize[1]
    instruments = zip(simulation.stokesfitspaths(), simulation.instrumentnames())
    print 'with {} Instruments...'.format(len(instruments))
    for (fileI, fileQ, fileU, fileV), name in instruments:
        print '  Creating PDF polarization map(s) for instrument "{}"'.format(name)

        # read the Stokes frames with shape (ny, nx) or datacubes with shape (nlambda, ny, nx)
        Is = pyfits.getdata(arch.openbinary(fileI))
        Qs = pyfits.getdata(arch.openbinary(fileQ))
        Us = pyfits.getdata(arch.openbinary(fileU))
        Vs = pyfits.getdata(arch.openbinary(fileV))

        # we don't need a full datacube for getting the binning information
        if len(Is.shape)==3:
            I = Is[0,:,:]
        else:
            I = Is

        # determine number of pixels dropped by binning and inform user
        orLenX = np.shape(I)[1]
        dropX =  orLenX % binX
        startX = dropX/2
        orLenY = np.shape(I)[0]
        dropY = orLenY % binY
        startY = dropY/2
        dropped = dropX * orLenY + dropY * orLenX - dropX * dropY
        droppedFraction = dropped / orLenX / orLenY
        print "    Picture dimensions: %s; bin size: %s by %s; dropping %s pixels (%.2f%%) (l:%s, r:%s, b:%s, t:%s)" \
            % (np.shape(I.T), binX, binY, dropped, 100*droppedFraction, startX, dropX - startX, startY, dropY-startY)

        # if all pixels are dropped, skip this instrument and inform user
        if droppedFraction==1:
            print "    Dropped all pixels! Skipping instrument."
            continue

        # create new directory to redirect all plots
        pathparts = fileQ.rsplit('/',1)
        pathparts[0] += '/_polarization'
        try:
            os.makedirs(pathparts[0])
        except OSError:
            if not os.path.isdir(pathparts[0]):
                raise
        path = pathparts[0]+'/'+pathparts[1] # the ending is still "Q.fits", will be changed when saving the figure(s)

        # determine the appropriate frame index or indices
        if wavelength == "all":
          if len(Is.shape)==3:
              indexes = range(len(Is[:,0,0]))
          else:
              indexes = {0}
        else:
            indexes = {simulation.frameindices([float(wavelength)])[0]}

        plotData = np.zeros(0)
        indexnr = 1
        for index in indexes:
            print '    Plotting polarization at {} um, {} of {}'.format(str(simulation.wavelengths()[index]),
                    indexnr, len(indexes))
            indexnr+= 1

            if len(Is.shape)==3:
                I = Is[index,:,:]
                Q = Qs[index,:,:]
                U = Us[index,:,:]
                V = Vs[index,:,:]

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

            # setup the figure
            plt.close('all')
            figure = plt.figure()

            # plot contour of intensity (in HD)
            Inan = I*1.
            Inan[Inan<=0] = np.NaN
            backGrndPlt = plt.imshow(Inan, alpha=0.4, norm=LogNorm(), origin = 'lower')
            #backGrndPlt = plt.contourf(Inan, alpha=0.4, locator=ticker.LogLocator()) # for higher quality
            cbar = plt.colorbar(backGrndPlt)
            cbar.set_label(simulation.fluxlabel(), labelpad=5)

            # determine a characteristic 'high' degree of polarization in the Frame (For key and scaling)
            charDegree = np.percentile(degreeLD, 99.0) # This has to be done before degreeLD contains 'np.NaN'
            # removing pixels with zero polarization
            degreeLD[degreeLD<=0] = np.NaN   # a runtime error here means that a laser is pointed at your instrument
            if not 0<charDegree<1:
                charDegree = np.nanmax((np.max(degreeLD), 0.01)) #  a warning means no polarization at all

            # calculate the scaling. It is relative to the picture size, thus the high polarization 'charDegree' should
            # be short enough so it does not overlap with the polarization around it.
            scale = (charDegree* 2.2)*max(len(posX)/figsize[0], len(posY)/figsize[1])
            keyLength = 10**(int(np.log10(charDegree))) # a runtime error here means no polarization in the whole picture
            key = "{:.3g}%".format(100 * keyLength)

            # create the polarization vector arrays
            xPolarization = degreeLD*np.sin(angle)
            yPolarization = degreeLD*np.cos(angle)

            # plot the vector field (in LD)
            X,Y = np.meshgrid( posX , posY)
            quiverplot = plt.quiver(X,Y, xPolarization, yPolarization, pivot='middle', units='inches',
                                    angles = 'xy', scale = scale, scale_units = 'inches', headwidth=0, headlength=1,
                                    headaxislength=1, minlength = 0.6)

            # add legend, labels and title
            plt.quiverkey(quiverplot, 0.85, 0.02, keyLength, key,
                          coordinates='axes', labelpos='E')
            plt.suptitle(simulation.prefix(), fontsize=14, fontweight='bold')
            plt.title('instrument: ' + name + ', $\lambda=' +
                        str(simulation.wavelengths()[index])+'\ \mu m$', fontsize=12)
            plt.axis('scaled')
            plt.xlabel('x (pixels)')
            plt.ylabel('y (pixels)')
            plt.xlim(xmin=0, xmax=len(I[0,:]))
            plt.ylim(ymin=0, ymax=len(I[:,0]))

            # save and close the figure
            fileEnd = "pol_" + str(simulation.wavelengths()[index]) + "um.pdf"
            plotfile = path.replace("stokesQ.fits",fileEnd)
            if not noMaps:
                plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)
                print "      Created PDF polarization map "# + plotfile
            plt.close()

            #to plot the polarization degree averaged over y-direction for all x-pixels
            if polAvY:
                figure = plt.figure()
                degreeHD = np.sqrt(np.average(Q, axis = 0)**2 + np.average(U, axis=0)**2)
                degreeHD[degreeHD>0] /= np.average(I,axis = 0)[degreeHD>0]
                plt.plot(100. * degreeHD)
                plt.suptitle(simulation.prefix(), fontsize=14, fontweight='bold')
                plt.title('instrument: ' + name + ', $\lambda=' +
                            str(simulation.wavelengths()[index])+'\ \mu m$', fontsize=12)
                plt.xlabel('x (pixels)')
                plt.ylabel('Average polarization (%)')
                polAvYfile = "polAvY".join(plotfile.rsplit("pol", 1))
                plt.savefig(polAvYfile, bbox_inches='tight', pad_inches=0.25)
                if export:
                    np.savetxt(polAvYfile.replace(".pdf", ".dat"), degreeHD,
                    header = 'column 1: y-axis averaged polarization degree for all x pixels')
                plt.close()
                print "      Created integrated polarization plot"

# -----------------------------------------------------------------
