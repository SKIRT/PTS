#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.image Work with a FITS image
#  bla

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import math
import numpy as np

# Modules for plotting
#if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import pyregion
#import pyds9

# Import astronomic modules
import astropy.io.fits as pyfits
#from photutils import CircularAperture
#from photutils import aperture_photometry

# Import PTS modules
from pts.log import Log

# -----------------------------------------------------------------

## An instance of the Image class represents a FITS image file
#
class Image:

    ## The constructor accepts the following arguments:
    #
    #  - path: the path of the FITS image file
    #
    def __init__(self, path):

        # TODO: allow just a filename
        # TODO: check whether the path points to a FITS file

        # Set the current path of the image
        self._path = path

        # Set the name of the image
        self._name = os.path.splitext(os.path.basename(self._path))[0]

        # Determine the path to the directory of the image
        self._directory = os.path.dirname(self._path)

        # Create a logger
        self._log = Log()

        # Show which image we are importing
        self._log.info("Reading in file: " + self._path)

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(self._path)

        # Get the primary image
        hdu = hdulist[0]

        # Get the image header
        self._header = hdu.header

        # Check for multiple planes with the PLANE keyword
        multiplanes = False
        try:

            plane_id = self._header['PLANE0']
            self._log.warning("Multiple planes detected. Using PLANE0 = " + plane_id)
            multiplanes = True

        # If the PLANE keyword is not found, we assume we only have one plane
        except KeyError: pass

        # Get the data from the HDU
        self._data = hdu.data[0] if multiplanes else hdu.data

        # Close the fits file
        hdulist.close()

        # Print dimension of data
        self._log.info("Dimensions of data cube: " + str(self._data.shape))

        # Print type of data
        self._log.info("Type of data: " + str(self._data.dtype.name))

        # Set the subtracted flag to False
        self._subtracted = False

        # Check whether background has already been subtracted from this image
        # Try to read the 'BACK_SUB' key from the header and set self._subtracted accordingly
        try:

            self._subtracted = self._header['BACK_SUB']

        except KeyError:

            #self._log.warning("Could not find information on subtraction; assuming image is not sky-subtracted yet")
            self._subtracted = False

        # Inform the user if this image has already been subtracted
        if self._subtracted: self._log.success("Background already subtracted for this image")

    ## This function returns the original path of this image
    @property
    def path(self):

        return self._path

    ## This function returns whether the image has been sky-subtracted or not
    @property
    def subtracted(self):

        return self._subtracted

    ## This function returns the number of pixels in the x direction
    @property
    def xsize(self):

        return self._data.shape[1]

    ## This function returns the number of pixels in the y direction
    @property
    def ysize(self):

        return self._data.shape[0]

    ## This function ..
    @property
    def mean(self):

        return np.mean(self._data)

    ## This function ..
    @property
    def median(self):

        return np.median(self._data)

    ## This function ..
    @property
    def min(self):

        return np.min(self._data)

    ## This function ..
    @property
    def max(self):

        return np.max(self._data)

    ## This function
    @property
    def stdev(self, corrected=False):

        # Set the delta degrees of freedom
        ddof = 1 if corrected else 0

        # Return the standard deviation of the data
        return np.std(self._data, ddof=ddof)

    ## This function saves the image at the specified path
    def save(self, path):

        # Create the HDU
        hdu = pyfits.PrimaryHDU(self._data, self._header)

        # Write to file
        hdu.writeto(path, clobber=True)

        # Update the image path
        self._path = path

    ## This function undoes the last operation, if possible
    def undo(self):

        # Check whether the previous state has been saved or not
        if self._prevdata is not None:

            # Replace the data with the previous data
            self._data = self._prevdata

            # Set the previous data to None
            self._prevdata = None

        else:

            self._log.warning("Cannot undo")

    ## This function
    def plot(self, path=None):

        # Plot the data using logaritmic scale
        plt.imshow(self._data, cmap='gray', norm=LogNorm())

        # Add a color bar
        plt.colorbar()

        # Display the result
        plt.show()

    ## This function
    def histogram(self, path=None):

        NBINS = 1000
        plt.hist(self._data.flat, NBINS)

        # Display the result
        plt.show()

    ## This function ..
    #
    #  - xrange: a slice object: slice(lowx, highx)
    #  - yrange: a slice object: slice(lowy, highy)
    def crop(self, xrange, yrange):

        # Save the current state of the image
        self._prevdata = self._data

        # Crop the image
        self._data = self._data[yrange, xrange]

    ## This function masks ...
    def maskregions(self, interpolate=False):

        # Save the current state of the image
        self._prevdata = self._data

        # TODO: automatically detect sources?

        # Determine the path to the regions file
        # TODO: look for a more general regions file
        regionspath = os.path.join(self._directory, self._name + "_mask.reg")

        # From the list of regions, calculate the mask
        regions = pyregion.open(regionspath)
        mask = regions.get_mask(shape=(self.ysize,self.xsize))

        # Mask the image with zeroes
        self._data[mask] = 0

        self.plot()

        # If we need to interpolate
        if interpolate:

            # Set some parameters
            order = 2
            linear = True

            # For each region, we interpolate within a box surrounding the region
            for region in regions:

                # Center and radius of this region
                xcenter = round(region.coord_list[0])
                ycenter = round(region.coord_list[1])
                radius  = round(region.coord_list[2]) + 1

                # Construct a box around this region
                xmin = int(xcenter -1 - 2*radius)
                xmax = int(xcenter -1 + 2*radius)
                ymin = int(ycenter -1 - 2*radius)
                ymax = int(ycenter -1 + 2*radius)
                xrange = slice(xmin,xmax)
                yrange = slice(ymin,ymax)
                box = self._data[yrange, xrange]

                # Create mask and boxminusmask
                mask = np.zeros_like(box)
                boxminusmask = np.ones_like(box)
                mask[box == 0] = 1
                boxminusmask[box == 0] = 0

                # Interpolate
                range = ((xmin,xmax),(ymin,ymax))
                interpolatedbox = self._interpolate(order,linear,range)

                # Fill the data array with the interpolated values
                self._data[yrange, xrange] = self._data[yrange, xrange]*boxminusmask + interpolatedbox*mask

        self.plot()

    ## This function ..
    def _interpolate(self, order, linear, limits=None):

        # Make lists to contain the coordinates and fluxes from which to interpolate from
        xvalues = []
        yvalues = []
        fluxes = []

        # Iterate over the image pixels
        for x in range(self._data.shape[1]):

            # Check whether this x value falls within the range
            if x < limits[0][0] or x > limits[0][1]: continue

            for y in range(self._data.shape[0]):

                # Check whether this y value falls within the range
                if y < limits[1][0] or y > limits[1][1]: continue

                # Get the flux in this pixel if it is not masked
                if self._data[y,x] != 0:

                    xvalues.append(x)
                    yvalues.append(y)
                    fluxes.append(self._data[y,x])

        # Make nummy arrays from x and y coordinates
        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)

        # Fit a 2D polynomial
        self._log.info("Fitting a 2D polynomial to the surrounding pixel values")
        m = polyfit2d(xvalues, yvalues, fluxes, order, linear)

        # Evaluate the polynomial in the specified range
        xx, yy = np.meshgrid(range(limits[0][0],limits[0][1]), range(limits[1][0],limits[1][1]))
        box_fluxes = polyval2d(xx, yy, m)

        # Return the result
        return box_fluxes

    ## This function subtracts the sky from this image ...
    #
    def subtractsky(self, order=3):

        # TODO: automatically determine background?

        # Determine the path to the regions file
        regionspath = os.path.join(self._directory, self._name + "_sky.reg")

        # If a regions file can not be found, look for a more general regions file
        if not os.path.isfile(regionspath):

            regionspath = os.path.join(self._directory, "sky.reg")

            if not os.path.isfile(regionspath):

                self._log.error("A regions file could not be found for this image")
                exit()

        # Create a list for all the sky regions defined for this image
        regions = getregions(regionspath)

        # Loop over all regions, calculate the median flux in each region
        xvalues = []
        yvalues = []
        fluxes = []

        for region in regions:

            # For regions that are circles
            if region[0] == "circle":

                # Get the center and radius for this circular region
                center = region[1]
                radius = region[2]

                # Create an aperture object
                aperture = CircularAperture(center, radius)

                # Calculate the photometry in the aperture
                fluxtable = aperture_photometry(self._data, aperture)

                # Get the sum of the flux within the aperture and the x and y coordinate
                asum = fluxtable['aperture_sum'][0]
                x = fluxtable['xcenter'][0][0]
                y = fluxtable['ycenter'][0][0]

                # Calculate a mean flux within the box
                flux = asum/(math.pi * radius**2)

                # Append the x coordinate, y coordinate and flux to the appropriate lists
                xvalues.append(x)
                yvalues.append(y)
                fluxes.append(flux)

            # For regions that are boxes
            elif region[0] == "box":

                # Get the center, width and height for this box region
                xcenter, ycenter = region[1]
                width = region[2]
                height = region[3]

                # Determine the coordinates of the corner pixels
                xmin = round(xcenter - 0.5*width)
                xmax = round(xcenter + 0.5*width)
                ymin = round(ycenter - 0.5*height)
                ymax = round(ycenter + 0.5*height)

                # Calculate the sum of the flux in the pixels within the box
                asum = sum(self._data[ymin:ymax,xmin:xmax])

                # Calculate a mean flux within the box
                flux = asum / (width*height)

                # Append the x coordinate, y coordinate and flux to the appropriate lists
                xvalues.append(xcenter)
                yvalues.append(ycenter)
                fluxes.append(flux)

        # Compyte mean background for plotting purposes
        meansky = np.mean(fluxes)

        self._log.info("Mean sky value = " + str(meansky))

        # Set some arguments
        plot = True
        linear = True
        xlimit = self._data.shape[1]
        ylimit = self._data.shape[0]

        # Make numpy arrays for x and y coordinates
        xvalues = np.copy(xvalues)
        yvalues = np.copy(yvalues)

        # Fit a 2D polynomial
        self._log.info("Fitting a 2D polynomial to the fluxes")
        m = polyfit2d(xvalues, yvalues, fluxes, order, linear)

        # Evaluate it on a grid...
        nx, ny = 20, 20
        xx, yy = np.meshgrid(np.linspace(0, xlimit, nx), np.linspace(0, ylimit, ny))
        zz = polyval2d(xx, yy, m)

        # Make a sky array and subtract it from the object frame
        sky = np.copy(self._data)
        strip = np.arange(1.0,np.float(xlimit+1.0))
        for i in range(0,ylimit):
            sky[i,:] = polyval2d(strip,i,m)

        if plot:

            # Plot polynomial
            plt.figure(1)
            plt.subplot(121)
            plt.imshow(sky[::-1,:], vmin=0.5*meansky, vmax=10.0*meansky, cmap='hot', extent=(0,xlimit,0,ylimit))
            plt.colorbar(orientation='horizontal')
            plt.scatter(xvalues, yvalues, c='k', marker='+')
            plt.xlim(0,xlimit)
            plt.ylim(0,ylimit)

            # Plot image
            plt.subplot(122)
            plt.imshow(self._data[::-1,:], vmin=0.5*meansky, vmax=10.0*meansky, cmap='hot', extent=(0,xlimit,0,ylimit))
            plt.xlim(0,xlimit)
            plt.ylim(0,ylimit)
            plt.colorbar(orientation='horizontal')
            plt.show()

        #subtractiondir = os.path.join(self._directory, "Subtraction")

        #skypath = os.path.join(subtractiondir, "sky_" + self._name + ".fits")

        # Create sky HDU and write to disk
        #hdu = pyfits.PrimaryHDU(sky)
        #hdu.writeto(skypath, clobber=True)

        # Subtract the sky from the data
        self._prevdata = self._data
        self._data = self._prevdata - sky

        # Path of the new, sky-subtracted image
        #newimagepath = os.path.join(subtractiondir, self._name + ".fits")

        #hdu = pyfits.PrimaryHDU(newdata, self._header)
        #hdu.writeto(newimagepath, clobber=True)

# Two-dimensional polynomial fit. Based upon code provided by Joe Kington, see:
# http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent/7997925#7997925
def polyfit2d(x, y, z, order, linear):

    import itertools

    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
        if linear & (i != 0.) & (j != 0.):
            G[:, k] = 0
    m, _, _, _ = np.linalg.lstsq(G, z)

    return m

# Values to two-dimensional polynomial fit.
# Based uppon code provided by Joe Kington
def polyval2d(x, y, m):

    import itertools

    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

## This function ...
def getregions(regionspath):

    # Create an empty list to contain the regions
    regions = []

    # Open the regions file and read each line
    with open(regionspath) as regionsfile:

        # Split the file in lines
        lines = regionsfile.readlines()

        # For each line
        for line in lines:

            # If the line is not a comment line
            if not line.startswith("#"):

                words = line.split(",")

                shape = words[0].split("(")[0]

                if shape == "circle":

                    center = (float(words[0].split("(")[1]), float(words[1]))
                    radius = float(words[2].split(")")[0])

                    regions.append((shape, center, radius))

                elif shape == "box":

                    center = (float(words[0].split("(")[1]), float(words[1]))
                    width = float(words[2])
                    height = float(words[3])

                    regions.append((shape, center, width, height))

    # Return the list of regions found in the file
    return regions