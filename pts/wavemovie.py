#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.wavemovie Creating a movie that runs through all wavelengths in the SKIRT simulation output.
#
# The function in this module creates a movie for the output of a SKIRT simulation. The movie combines
# the SEDs (bottom panel) and the pixel frames (top panel, from left to right) for up to three instruments,
# running through all wavelengths in the simulation.

# -----------------------------------------------------------------

import numpy as np
import matplotlib
if matplotlib.get_backend().lower() != "agg": matplotlib.use("agg")
import matplotlib.pyplot as plt

from pts.moviefile import MovieFile
from pts.rgbimage import RGBImage

# -----------------------------------------------------------------

# This function creates a "wavelength" movie for the output of a SKIRT simulation. The movie combines
# the SEDs (bottom panel) and the pixel frames (top panel, from left to right) for up to three instruments,
# running through all wavelengths in the simulation.
#
# The function takes the following arguments:
#  - sedpaths:    a list of the filepaths to the SED files; the list must have 1, 2 or 3 items
#  - fitspaths:   a list of the filepaths to the FITS files; the list must have same number of items as sedpaths
#  - labels:      a list of human-readable labels; the list must have same number of items as sedpaths
#  - outpath:     the filepath to the output file; must have the ".mov" filename extension
#  - wavelengths: a list of the wavelengths in the simulation, in micron;
#                 the length of the list must match the number of wavelengths in each of the SED and fits files
#  - xlim, ylim:  the lower and upper limits of the x/y axis of the SED plot, specified as a 2-tuple;
#                 if missing the corresponding axis is auto-scaled
#  - from_percentile, to_percentile: the percentile values, in range [0,100], used to clip the luminosity values
#                 loaded from the fits files; the default values are 30 and 100 respectively
#
def makewavemovie(sedpaths, fitspaths, labels, outpath, wavelengths,
                  xlim=None, ylim=None, from_percentile=30, to_percentile=100):

    # verify some aspects of the arguments
    npaths = len(sedpaths)
    nlambda = len(wavelengths)
    if npaths<1 or npaths>3 or len(fitspaths)!=npaths or len(labels)!=npaths:
        raise ValueError("There should be 1, 2 or 3 SEDs, FITS files and labels")
    if len(wavelengths)<4:
        raise ValueError("There should be at least 4 wavelengths")

    # get the image shapes (assume that fits frames in all files have the same shape)
    fitsshape = RGBImage(fitspaths[0]).shape
    sedshape = (npaths*fitsshape[0], fitsshape[1]/2)
    totalshape = (sedshape[0], fitsshape[1]+sedshape[1])

    # determine the appropriate pixel range for ALL images
    print "  preprocessing frames for " + outpath.rsplit("/",1)[1] + "..."
    ranges = []
    for frame in range(nlambda):
        for fitspath in fitspaths:
            im = RGBImage(fitspath, frames=(frame,frame,frame))
            ranges += list(im.percentilepixelrange(from_percentile,to_percentile))
    rmin = min(ranges)
    rmax = max(ranges)

    # open the movie file
    movie = MovieFile(outpath, shape=totalshape, rate=10)

    # for each wavelength, add a movie frame
    for frame in range(nlambda):
        print "  adding frame " + str(frame)+"/"+str(nlambda) + "..."

        # assemble the top panel
        image = None;
        for fitspath in fitspaths:
            im = RGBImage(fitspath, frames=(frame,frame,frame))
            im.setrange(rmin,rmax)
            im.applylog()
            im.applycmap("gnuplot")
            if image==None: image = im
            else: image.addright(im)

        # plot the seds
        dpi = 100
        figure = plt.figure(dpi=dpi, figsize=(sedshape[0]/float(dpi),sedshape[1]/float(dpi)), facecolor='w', edgecolor='w')
        for sedpath,label in zip(sedpaths,labels):
            data = np.loadtxt(sedpath)
            plt.loglog(data[:,0], data[:,1], label=label)
        plt.axvline(wavelengths[frame], color='m')
        plt.figtext(0.49,0.15, r"$\lambda={0:.4g}\,\mu m$".format(wavelengths[frame]), fontsize='large')
        plt.ylabel(getfluxlabel(sedpaths[0]), fontsize='large')
        plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        if xlim != None: plt.xlim(xlim)
        if ylim != None: plt.ylim(ylim)
        plt.legend(loc='best')
        plt.draw()     # flush the drawing
        im = RGBImage(figure)
        plt.close()

        # assemble the final image and add it as a movie frame
        image.addbelow(im)
        image.addto(movie)

    # close the movie file
    movie.close()
    return True

# -----------------------------------------------------------------

## This function returns an appropriate axis label for the flux described in the second line of the specified file.
def getfluxlabel(sedpath):
    # get the second line of the file, which contains the description of the flux column
    sedfile = open(sedpath, 'r')
    sedfile.readline()
    fluxdescription = sedfile.readline()
    sedfile.close()

    # select the appropriate label based on the units given in the description
    if "Jy" in fluxdescription: return r"$F_\nu\,(\mathrm{Jy})$"
    if "W/m2" in fluxdescription: return r"$\lambda\,F_\lambda\,(\mathrm{W}\,\mathrm{m}^{-2})$"
    return "Flux"

# -----------------------------------------------------------------
