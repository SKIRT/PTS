#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.polyfitmovie Create a movie of the progress during a polychromatic FitSKIRT fit.
#
# This script creates a movie of the progress during a polychromatic FitSKIRT fit.
# The in/out filenames and other parameters are hardcoded in the script.
# The script is a good example of using the pts.rgbimage and pts.moviefile modules.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.rgbimage import RGBImage
from pts.core.basics.moviefile import MovieFile

# ------------------------------------------------------------------

# wavelengths and names of reference images used during fit; number of output simulations
waves = ["0.345","0.475","0.622","0.763","0.905","1.22","1.63","2.19"]
refs=["U_NORM.fits","G_NORM.fits","R_NORM.fits","I_NORM.fits","Z_NORM.fits","J_extr.fits","H_extr.fits","K_NORM.fits"]
numSimulations = 88

# ------------------------------------------------------------------

# build a list of reference images
refImages = []
for ref in refs:
     im = RGBImage(ref)
     im.setrange(0,8e-5)
     im.applycmap("afmhot")
     refImages.append(im)

# make all reference images the same shape, and remember the final shape of a single image
for i in refImages:
     for j in refImages:
          i.enlargecanvas(j)
imageShape = refImages[0].shape

# stack the reference images on top of each other
fullRef = refImages[0]
for i in range(1,len(refImages)):
     fullRef.addbelow(refImages[i])

# add a column with text labels to the right of the reference stack
fullRef.addright(RGBImage("ugrizjhk.jpg"))

# create a movie, frame by frame
movie = MovieFile("fit.mov", shape=(fullRef.shape[0]+imageShape[0], fullRef.shape[1]), rate=6)
for i in range(numSimulations):
    print("Starting frame", i+1)

    # build a list of simulated images for this frame
    simImages = []
    for j in waves:
        im = RGBImage("Final_result_"+j+"_"+str(i)+".fits")
        im.setrange(0,8e-5)
        im.applycmap("afmhot")
        simImages.append(im)

    # make all simulated images the same shape as the largest reference image
    for j in simImages:
        j.enlargecanvas(imageShape)

    # stack the simulated images on top of each other, and scale the values to match the JPEG image
    fullSim = simImages[0]
    for j in range(1,len(simImages)):
        fullSim.addbelow(simImages[j])

    # glue the simulated stack to the right of a clone of the reference stack
    fullFrame = RGBImage(fullRef.pixelarray())
    fullFrame.addright(fullSim)

    # add the frame twice because movie files don't support lower frame rates
    fullFrame.addto(movie)
    fullFrame.addto(movie)

movie.close()

# ------------------------------------------------------------------
