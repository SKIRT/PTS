#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.makergbimages Creating RGB images for SKIRT simulation output.
#
# The function in this module creates RGB images for the "total.fits" files in the output of a SKIRT simulation.
# The caller can specify the wavelengths corresponding to the R,G,B frames in the image.

# -----------------------------------------------------------------

import numpy as np
from pts.rgbimage import RGBImage

# -----------------------------------------------------------------

# This function creates one or more RGB images for each "total.fits" file in the output of the specified simulation.
# If the \em wavelength_tuples argument is missing, a single image is created
# for each "total.fits" file using the frames loaded by default by the RGBimage constructor. Otherwise,
# the \em wavelength_tuples argument must contain a sequence of 3-tuples with (R,G,B) wavelengths; each of
# these tuples causes an image to be created, loading the frames corresponding to the specified wavelengths.
# The \em from_percentile and \em to_percentile arguments take the percentile values, in range [0,100], used
# to clip the luminosity values loaded from the fits file.
# The images are saved in PNG format and are placed next to the original file(s) with the same name
# but a different extension. If there are multiple images per fits file, a serial number is added.
def makergbimages(simulation, wavelength_tuples=None, from_percentile=30, to_percentile=100):

    # loop over the wavelength tuples
    if wavelength_tuples==None: wavelength_tuples = [ None ]
    for index in range(len(wavelength_tuples)):

        # get the frame indices corresponding to the requested wavelengths
        if wavelength_tuples[index] == None: frames = None
        else: frames = simulation.frameindices(wavelength_tuples[index])

        # get a list of output file names, including extension, one for each instrument
        outnames = simulation.totalfitspaths()
        if len(outnames) > 0:

            # determine the appropriate pixel range for ALL output images for this galaxy
            ranges = []
            for outname in outnames:
                im = RGBImage(outname, frames=frames)
                ranges += list(im.percentilepixelrange(from_percentile,to_percentile))
            rmin = min(ranges)
            rmax = max(ranges)

            # create an RGB file for each output file
            for outname in outnames:
                im = RGBImage(outname, frames=frames)
                im.setrange(rmin,rmax)
                im.applylog()
                im.applycurve()
                savename = outname[:-5] + (str(index+1) if index > 0 else "") + ".png"
                im.saveto(savename)
                print "Created RGB image file " + savename

# -----------------------------------------------------------------
