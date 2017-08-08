#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.wavemovie Creating a movie that runs through all wavelengths in the SKIRT simulation output.
#
# The function in this module creates a movie for the output of a SKIRT simulation. The movie combines
# the SEDs (bottom panel) and the pixel frames (top panel, from left to right) for up to three instruments,
# running through all wavelengths in the simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import logger
from ..basics.log import log

# Import standard modules
import warnings
import os
import numpy as np
import matplotlib
with warnings.catch_warnings():
    warnings.filterwarnings('error')
    try:
        if matplotlib.get_backend().lower() != "agg": matplotlib.use("agg")
    except Warning as w: log.warning("An failed attempt of setting the Matplotlib backend has been made because it has already been set")
import matplotlib.pyplot as plt

# Import astronomical modules
try: import pyfits
except ImportError: import astropy.io.fits as pyfits

# Import the relevant PTS classes and modules
from ..basics.moviefile import MovieFile
from ..basics.rgbimage import RGBImage
from ..tools import archive as arch

# -----------------------------------------------------------------

# This function creates a movie for the output of the specified simulation. The movie combines
# the SEDs (bottom panel) and the pixel frames (top panel, from left to right) for up to three instruments,
# running through all wavelengths in the simulation. The movie is placed next to the original file(s) with
# a similar name (omitting the instrument name) but a different extension.
# The function takes the following arguments:
#  - simulation:  the SkirtSimulation object representing the simulation to be handled
#  - xlim, ylim:  the lower and upper limits of the x/y axis of the SED plot, specified as a 2-tuple;
#                 if missing the corresponding axis is auto-scaled
#  - from_percentile, to_percentile: the percentile values, in range [0,100], used to clip the luminosity values
#                 loaded from the fits files; the default values are 30 and 100 respectively
def makewavemovie(simulation, xlim=None, ylim=None, from_percentile=30, to_percentile=100, output_path=None):
    sedpaths = simulation.seddatpaths()
    npaths = len(sedpaths)
    if 1 <= npaths <= 3:
        fitspaths = simulation.totalfitspaths()
        wavelengths = simulation.wavelengths()
        nlambda = len(wavelengths)
        if len(fitspaths) == npaths and nlambda > 3:
            labels = [ path.rsplit("_",2)[1] for path in sedpaths ]
            outpath = sedpaths[0].rsplit("_",2)[0] + "_wave.mov"
            if output_path is not None: outpath = os.path.join(output_path, os.path.basename(outpath))
            fluxlabel = simulation.fluxlabel()

            # get the image shapes (assume that fits frames in all files have the same shape)
            fitsshape = simulation.instrumentshape()
            sedshape = (npaths*fitsshape[0], fitsshape[1]/2)
            totalshape = (sedshape[0], fitsshape[1]+sedshape[1])

            # load the data
            print("  loading data for " + outpath.rsplit("/",1)[1] + "...")
            datacubes = [ pyfits.getdata(arch.openbinary(fitspath)).T for fitspath in fitspaths ]
            sedtables = [ np.loadtxt(arch.opentext(sedpath)) for sedpath in sedpaths ]

            # determine the appropriate pixel range for ALL images
            print("  preprocessing frames...")
            ranges = []
            for frame in range(nlambda):
                for data in datacubes:
                    im = RGBImage(np.dstack(( data[:,:,frame],data[:,:,frame],data[:,:,frame] )))
                    ranges += list(im.percentilepixelrange(from_percentile,to_percentile))
            rmin = min(ranges)
            rmax = max(ranges)

            # open the movie file
            movie = MovieFile(outpath, shape=totalshape, rate=10)

            # for each wavelength, add a movie frame
            for frame in range(nlambda):
                print("  adding frame " + str(frame+1)+"/"+str(nlambda) + "...")

                # assemble the top panel
                image = None
                for data in datacubes:
                    im = RGBImage(np.dstack(( data[:,:,frame],data[:,:,frame],data[:,:,frame] )))
                    im.setrange(rmin,rmax)
                    im.applylog()
                    im.applycmap("gnuplot")
                    if image==None: image = im
                    else: image.addright(im)

                # plot the seds
                dpi = 100
                figure = plt.figure(dpi=dpi, figsize=(sedshape[0]/float(dpi),sedshape[1]/float(dpi)),
                                    facecolor='w', edgecolor='w')
                for data,label in zip(sedtables,labels):
                    plt.loglog(data[:,0], data[:,1], label=label)
                plt.axvline(wavelengths[frame], color='m')
                plt.figtext(0.49,0.15, r"$\lambda={0:.4g}\,\mu m$".format(wavelengths[frame]), fontsize='large')
                plt.ylabel(fluxlabel, fontsize='large')
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
            print("Created wave movie file " + outpath)

# -----------------------------------------------------------------
