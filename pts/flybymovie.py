#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.flybymovie Creating a flyby movie for a SKIRT model.
#
# The main function in this module creates a flyby movie for the model specified in a particular ski file.
# The other functions each perform a single stage of the action; they can be called directly if needed.

# -----------------------------------------------------------------

import numpy as np
import os
import os.path

from pts.skirtsimulation import SkirtSimulation
from pts.skirtexec import SkirtExec
from pts.skifile import SkiFile
from pts.rgbimage import RGBImage
from pts.moviefile import MovieFile

# -----------------------------------------------------------------

## This function creates a flyby movie for the simulation model described in the specified ski file,
# according to the specifications contained in the provided timeline object.
#
# The function adds the appropriate instruments to the specified ski file, launches SKIRT to
# perform the simulation, and combines the resulting fits files into a single RGB movie file.
# For optimal results, the ski file should specify a wavelength grid with 3 wavelengths corresponding to
# Blue, Green and Red.
#
# Parameters:
#  - skirtpath: absolute or relative path to the skirt executable
#  - skifile: absolute or relative file path of the ski file
#  - timeline: a pts.timeline.Timeline object containing the specifications for the movie
#  - moviefile: absolute or relative file path of the movie file to be created
#  - from_percentile and to_percentile: the percentile values, in range [0,100], used to clip the luminosity
#    values loaded from the fits file (see the RGBImage.applylog() function)
#  - contrast: if True, the contrast of each frame is enhanced (nice, but it takes quite a while); default is False
#
# Important notes:
#  - the specified ski file is adjusted by this function (replacing the instruments)
#  - any simulation input files should reside in the same folder as the ski file
#  - a simulation output folder called "flyby_out" is created next to the ski file
#
def createflybymovie(skirtpath, skifile, timeline, moviefile, from_percentile=30, to_percentile=100, contrast=False):
    adjustskifile(skifile, timeline)
    simulation = execskirt(skirtpath, skifile)
    createmovie(simulation, moviefile, timeline.rate, from_percentile, to_percentile, contrast)

# -----------------------------------------------------------------

# This function adds the appropriate instruments to the specified ski file to create a flyby movie
# according to the specifications contained in the provided timeline object.
# For optimal results, the ski file should specify a wavelength grid with 3 wavelengths corresponding to
# Blue, Green and Red.
#
# Parameters:
#  - skifile: absolute or relative file path of the ski file
#  - timeline: a pts.timeline.Timeline object containing the specifications for the movie
#
# Important notes:
#  - the specified ski file is adjusted by this function (replacing the instruments)
#
def adjustskifile(skifile, timeline):
    # insert the appropriate instruments in the ski file
    ski = SkiFile(skifile)
    ski.setperspectiveinstruments(timeline.getframes())
    ski.saveto(skifile)

## This function launches SKIRT to perform the simulation for the specified ski file, storing the results
# in a temporary output folder called "flyby_out" next to the ski file. Any simulation input files should
# reside in the same folder as the ski file.
#
# Parameters:
#  - skirtpath: absolute or relative path to the skirt executable
#  - skifile: absolute or relative file path of the ski file
#
# Return value: the SkirtSimulation object returned by the SkirtExec.execute() function
#
def execskirt(skirtpath, skifile):
    skirt = SkirtExec(skirtpath)
    outpath = os.path.join(os.path.dirname(skifile), "flyby_out")
    if not os.path.exists(outpath): os.mkdir(outpath)
    simulation = skirt.execute(skifile, inpath="", outpath=outpath, skirel=True)[0]
    if simulation.status() != 'Finished':  raise ValueError("Simulation " + simulation.status())
    return simulation

## This function combines the fits files resulting from the specified SKIRT simulation
# into a single RGB movie file.
#
# Parameters:
#  - simulation: a SkirtSimulation object, e.g. as returned by the SkirtExec.execute() function,
#    or the absolute or relative file path of the simulation output directory
#  - moviefile: absolute or relative file path of the movie file to be created
#  - rate: the number of frames per second (default is 24)
#  - from_percentile and to_percentile: the percentile values, in range [0,100], used to clip the luminosity
#    values loaded from the fits file (see the RGBImage.applylog() function)
#  - contrast: if True, the contrast of each frame is enhanced (nice, but it takes quite a while); default is False
#
# Important notes:
#  - the specified ski file is adjusted by this function (replacing the instruments)
#  - any simulation input files should reside in the same folder as the ski file
#  - a simulation output folder called "flyby_out" is created next to the ski file
#
def createmovie(simulation, moviefile, rate=24, from_percentile=30, to_percentile=100, contrast=False):
    # verify the simulation status
    if not isinstance(simulation,SkirtSimulation): simulation = SkirtSimulation(outpath=simulation)
    if simulation.status() != 'Finished':  raise ValueError("Simulation " + simulation.status())

    # determine the appropriate pixel range for ALL images
    ranges = []
    for fits in simulation.totalfitspaths():
        im = RGBImage(fits)
        im.applylog(from_percentile, to_percentile)
        ranges += list(im.pixelrange())
    rmin = min(ranges)
    rmax = max(ranges)

    # create the movie file
    movie = MovieFile(moviefile, shape=simulation.instrumentshape(), rate=rate)
    for fits in simulation.totalfitspaths():
        im = RGBImage(fits)
        im.applylog()
        im.setrange(rmin,rmax)
        if contrast: im.applycurve()
        im.addto(movie)
    movie.close()

# -----------------------------------------------------------------
