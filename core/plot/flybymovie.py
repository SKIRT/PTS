#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.flybymovie Creating a flyby movie for a SKIRT model.
#
# The main function in this module creates a flyby movie for the model specified in a particular ski file.
# The other functions each perform a single stage of the action; they can be called directly if needed.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..simulation.skifile import SkiFile
from ..basics.rgbimage import RGBImage
from ..basics.moviefile import MovieFile

# -----------------------------------------------------------------

# This function adds the appropriate instruments to the specified ski file to create a flyby movie
# according to the specifications contained in the provided timeline object.
# For optimal results, the ski file should specify a wavelength grid with three wavelengths corresponding to
# Blue, Green and Red.
#
# Parameters:
#  - skifile: absolute or relative file path of the ski file
#  - timeline: a pts.timeline.Timeline object containing the specifications for the movie
#
# Important notes:
#  - the specified ski file is adjusted by this function (replacing the instruments)
#
def prepareflybymovie(skifile, timeline):
    # insert the appropriate instruments in the ski file
    ski = SkiFile(skifile)
    ski.setperspectiveinstruments(timeline.getframes())
    ski.saveto(skifile)

# -----------------------------------------------------------------

## This function combines the fits files resulting from the specified SKIRT simulation into a single RGB movie file.
# The resulting file is saved next to the fits files with a similar name and the \c .mov filename extension.
#
# Parameters:
#  - simulation: a SkirtSimulation object representing the simuluation to be handled
#  - rate: the number of frames per second (default is 24)
#  - wavelengths: a 3-tuple of wavelengths to be used for the RGB frames;
#                 if missing or None the RGBImage default frames are used
#  - from_percentile and to_percentile: the percentile values, in range [0,100], used to clip the luminosity
#    values loaded from the fits file
#  - contrast: if True, the contrast of each frame is enhanced (nice, but it takes quite a while); default is False
#
def makeflybymovie(simulation, rate=24, wavelengths=None, from_percentile=30, to_percentile=100, contrast=False):
    # verify the simulation status
    if simulation.status() != 'Finished':  raise ValueError("Simulation " + simulation.status())
    fitspaths = simulation.totalfitspaths()
    outpath = simulation.outfilepath("flyby.mov")

    # find the appropriate wavelength indices in the data cube
    indices = simulation.frameindices(wavelengths) if wavelengths is not None else None

    # determine the appropriate pixel range for ALL images
    print("  preprocessing frames for " + outpath.rsplit("/",1)[1] + "...")
    ranges = []
    for fits in fitspaths:
        im = RGBImage(fits, frames=indices)
        ranges += list(im.percentilepixelrange(from_percentile,to_percentile))
    rmin = min(ranges)
    rmax = max(ranges)

    # create the movie file
    movie = MovieFile(outpath, shape=simulation.instrumentshape(), rate=rate)
    nframes = len(fitspaths)
    for frame in range(nframes):
        print("  adding frame " + str(frame+1) + "/" + str(nframes) + "...")
        im = RGBImage(fitspaths[frame], frames=indices)
        im.setrange(rmin,rmax)
        im.applylog()
        #im.applysqrt()
        if contrast: im.applycurve()
        im.addto(movie)
    movie.close()

# -----------------------------------------------------------------
