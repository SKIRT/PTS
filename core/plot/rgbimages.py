#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.rgbimages Creating RGB images for SKIRT simulation output.
#
# The function in this module creates RGB images for the "total.fits" files in the output of a SKIRT simulation.
# The caller can specify the wavelengths corresponding to the R,G,B frames in the image.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import warnings
import numpy as np

# Import astronomical modules
try: import pyfits
except ImportError: import astropy.io.fits as pyfits

# Import the relevant PTS classes and modules
from ..basics.rgbimage import RGBImage
from ..tools import archive as arch

# -----------------------------------------------------------------

## This function creates one or more RGB images for each "total.fits" file in the output of the specified simulation.
# If the \em wavelength_tuples argument is missing, a single image is created
# for each "total.fits" file using the frames loaded by default by the RGBimage constructor. Otherwise,
# the \em wavelength_tuples argument must contain a sequence of 3-tuples with (R,G,B) wavelengths; each of
# these tuples causes an image to be created, loading the frames corresponding to the specified wavelengths.
# The \em from_percentile and \em to_percentile arguments take the percentile values, in range [0,100], used
# to clip the luminosity values loaded from the fits file.
# The images are saved in PNG format and are placed next to the original file(s) with the same name
# but a different extension. If there are multiple images per fits file, a serial number is added.
def makergbimages(simulation, wavelength_tuples=None, from_percentile=30, to_percentile=100, output_path=None):

    # loop over the wavelength tuples
    if wavelength_tuples is None: wavelength_tuples = [None]

    # Loop over the tuples
    for index in range(len(wavelength_tuples)):

        # Tuples in a list
        if isinstance(wavelength_tuples, list):

            # get the frame indices corresponding to the requested wavelengths
            if wavelength_tuples[index] is None: frames = None
            else: frames = simulation.frameindices(wavelength_tuples[index])

            # Set name
            if index > 0: name = "_" + str(index + 1)
            else: name = ""

        # Tuples in a dict
        elif isinstance(wavelength_tuples, dict):

            name = "_" + wavelength_tuples.keys()[index]
            wavelengths = wavelength_tuples.values()[index]
            frames = simulation.frameindices(wavelengths)

        # Invalid
        else: raise ValueError("Invalid value for 'wavelength_tuples'")

        # get a list of output file names, including extension, one for each instrument
        outnames = simulation.totalfitspaths()
        if len(outnames) > 0:

            # determine the appropriate pixel range for ALL output images for this galaxy
            ranges = []
            for outname in outnames:
                try:
                    im = RGBImage(outname, frames=frames)
                    ranges += list(im.percentilepixelrange(from_percentile, to_percentile))
                except ValueError:
                    warnings.warn("Something is wrong with " + outname + ": skipping")
                    outnames.remove(outname)
                    continue

            rmin = min(ranges)
            rmax = max(ranges)

            # create an RGB file for each output file
            for outname in outnames:
                im = RGBImage(outname, frames=frames)
                im.setrange(rmin,rmax)
                im.applylog()
                im.applycurve()
                savename = outname[:-5] + name + ".png"
                if output_path is not None: savename = os.path.join(output_path, os.path.basename(savename))
                im.saveto(savename)
                print("Created RGB image file " + savename)

# -----------------------------------------------------------------

## This function creates an RGB image for each "total.fits" file in the output of the specified simulation,
# integrating the pixel values for each color over multiple output frames according to the specified filters.
# The images are saved in PNG format and by default are placed next to the original file(s) with the same name
# but a different extension, and optionally including an extra postfix.
#
# The function takes the following arguments:
#  - \em simulation: a SkirtSimulation instance representing the simulation for which to make images.
#  - \em filterspecs: a sequence of 4-tuples (filter, r, g, b), each containing a Filter instance and three
#        weights that specify the contribution of the fluxes integrated over this filter to each RGB channel.
#  - \em postfix: a string that will be added to the standard image file name; defaults to the empty string.
#  - \em fmax: the largest flux value that will be shown in the image without clipping;
#        if None or missing the function uses the largest flux value found in any of the channels/images.
#  - \em fmin: the smallest flux value that will be shown in the image without clipping;
#        if None or missing the function determines this value from \em fmax and \em frange
#        using the formula \f$f_\mathrm{range}=\log_{10}(f_\mathrm{max}/f_\mathrm{min})\f$.
#  - \em frange: the dynamic range for the flux values, expressed in order of magnitudes (decades);
#        the default value is 3. If \em fmin is specified, the value of \em frange is ignored.
#  - \em output_path: path to the PNG output directory; default is next to the fits files.
#
def makeintegratedrgbimages(simulation, filterspecs, postfix="", fmax=None, fmin=None, frange=3, output_path=None):

    # get the wavelength grid
    wavelengths = simulation.wavelengths()

    # get a list of output file names, including extension, one for each instrument
    outnames = simulation.totalfitspaths()
    if len(outnames) > 0:

        # construct an image object with integrated color channels for each output file
        # and keep track of the largest flux (in Jansky) in any of the images and channels
        images = []
        fluxmax = 0
        for outname in outnames:

            # get the data cube in per wavelength units
            cube = pyfits.getdata(arch.openbinary(outname)).T
            cube = simulation.convert(cube, to_unit='W/m3', quantity='fluxdensity', wavelength=wavelengths)

            # initialize an RGB frame
            dataRGB = np.zeros( (cube.shape[0], cube.shape[1], 3) )

            # add color for each filter
            for filter,w0,w1,w2 in filterspecs:
                data = filter.convolve(wavelengths, cube)
                data = simulation.convert(data, from_unit='W/m3', to_unit='Jy', wavelength=filter.pivotwavelength())
                dataRGB[:,:,0] += w0*data
                dataRGB[:,:,1] += w1*data
                dataRGB[:,:,2] += w2*data

            # construct the image object
            fluxmax = max(fluxmax, dataRGB.max())
            images.append(RGBImage(dataRGB))

        # determine the appropriate pixel range for all output images
        if fmax==None:
            fmax = fluxmax
        if fmin==None:
            fmin = fmax / 10**frange

        # create an RGB file for each output file
        for outname,im in zip(outnames,images):
            im.setrange(fmin,fmax)
            im.applylog()
            im.applycurve()
            savename = outname[:-5] + postfix + ".png"
            if output_path is not None: savename = os.path.join(output_path, os.path.basename(savename))
            im.saveto(savename)
            print("Created integrated RGB image file " + savename)

        # return the flux range used for these images
        return (fmin,fmax)

# -----------------------------------------------------------------
