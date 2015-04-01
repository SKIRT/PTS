#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.preparedata Prepare data for use with SKIRT or FitSKIRT
#

# -----------------------------------------------------------------

# Import standard modules
import os.path

# Import relevant PTS modules
from pts.image import Image
from pts.log import Log

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Get the full path to the data and preparation directories
datapath = os.path.join(os.getcwd(), "data")
preppath = os.path.join(os.getcwd(), "prep")

## This function
def importimage(filter):

    # Show which image we are importing
    log.info("Importing image for " + filter + " band")

    # Determine the path to this FITS file
    filepath = os.path.join(datapath, filter + ".fits")

    # Create an image object from this FITS file and return it
    return Image(filepath)

## This function masks the nans, edges
def mask(image, edges=True, extra=False, write=False):

    # Create a mask for the pixels in the primary image that have a NaN value
    image.masknans()

    # Select the nans mask for later
    image.masks.nans.select()

    if edges:

        # Assuming the nans are in the outer parts of the image, make an expanded mask that also covers the edges
        image.expandmasks("edges")

        # Deselect the nans mask
        image.masks.nans.deselect()

        # Select the edges mask for later
        image.masks.edges.select()

    if extra:

        # Import 'extra' region
        filepath = "data/extra/" + image.name + ".reg"
        image.importregion(filepath, "extra")

        # Select the extra region
        image.regions.extra.select()

        # Create a mask 'extra' from the 'extra' region
        image.createmask()

        # Deselect the extra region
        image.regions.extra.deselect()

        # Select the extra mask for later
        image.masks.extra.select()

    # Combine the edge and extra mask (the currently active masks)
    image.combinemasks("total")

    # Select the total mask
    image.masks.deselectall()
    image.masks.total.select()

    # Apply the total mask to the primary image
    image.applymasks()

    # Deselect all regions and masks
    image.masks.deselectall()
    image.regions.deselectall()

    # If requested, save the masked primary image
    if write:

        path = os.path.join(preppath, image.name, "masked.fits")
        image.save(path)

## This function interpolates over the stars
def interpolatestars(image, write=False):

    # Create a region object for the stars create a mask from it
    filepath = "data/stars/" + image.name + ".reg"
    image.importregion(filepath, "stars")

    # Select the stars region
    image.regions.stars.select()

    # Create a mask from the stars region
    image.createmask()

    # Select the stars mask
    image.masks.stars.select()

    # Interpolate the image within the stars mask
    image.interpolate()

    # If requested, save the interpolated image
    if write:

        path = os.path.join(preppath, image.name, "interpolated.fits")
        image.save(path)

## This function subtracts the sky of the specified image
def subtractsky(image, write=False):

    pass

## This function scales the image by a certain factor
def scale(image, factor, write=False):

    # Multiply the primary frame by the conversion factor
    image.multiply(factor)

    # If requested, save the scaled image
    if write:

        path = os.path.join(preppath, image.name, "scaled.fits")
        image.save(path)

## This function convolves the image with the PSF of the PACS 160 image
def convolve(image, pixelscale, kernel, write=False, writekernel=False):

    # Convolve to the PACS 160 resolution
    image.convolve(kernel, pixelscale)

    # If requested, save the convolved image
    if write:

        path = os.path.join(preppath, image.name, "convolved.fits")
        image.save(path)

    # If requested, save the new kernel
    if writekernel:

        # Select the kernel frame
        image.frames.deselectall()
        image.frames.kernel.select()

        path = os.path.join(preppath, image.name, "newkernel.fits")
        image.save(path)

## This function rebins the image to the resolution of the Pacs 160 micron image
def rebin(image, write=False):

    # Do the rebinning
    pacs160path = os.path.join(datapath, "PACS160.fits")
    image.rebin(pacs160path)

    # If requested, save the rebinned image
    if write:

        path = os.path.join(preppath, image.name, "rebinned.fits")
        image.save(path)

# -----------------------------------------------------------------

## OLD STARS: IRAC, YOUNG NI STARS: FUV, YOUNG I STARS: Ha + 24micron, DUST: H + PACS70 + PACS160

# -----------------------------------------------------------------






