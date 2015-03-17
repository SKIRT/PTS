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
from scipy import ndimage

# Modules for plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Import astronomic modules
import pyregion
import astropy.io.fits as pyfits
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils.background import Background

# Import PTS modules
from pts.mathematics import fitpolynomial, polynomial, fitgaussian, gaussian
from pts.log import Log
from pts.filter import Filter
from pts.galaxy import GalaxyFinder

# -----------------------------------------------------------------

## An instance of the Image class represents a FITS image file
#
class Image(object):

    ## The constructor accepts the following arguments:
    #
    #  - path: the path of the FITS image file
    #  - log: a Log instance
    #  - cut:
    #
    def __init__(self, path, log=None, cut=True):

        # TODO: allow just a filename (without path)
        # TODO: check whether the path points to a FITS file

        # Set the name of the image
        self._name = os.path.splitext(os.path.basename(path))[0]

        # Create a logger or use a pre-existing logger if possible
        self._log = Log() if log is None else log

        # Layers: primary, sky, fittedsky, errors ...
        self._layers = dict()

        # Masks: nans, edges, galaxy, stars
        self._masks = dict()

        # Regions: sky, stars,
        self._regions = dict()

        # Read in the data and the header. The data is saved as an ImageLayer in self._layers and
        # the header is stored as self.header
        self._import(path, cut)

        # Instrument
        #self._instrument = self._findinheader("INSTRUME")
        # Filter
        #self._filter = self._findinheader("FILTER")

        # TODO: Better way? Intelligent search in name for which filter?
        filtername = self._name

        # Get the filter
        self._filter = Filter(filtername)

        # Units
        self._units = self._findinheader("BUNIT")

        # The entries in this dictionary indicate whether certain operations have been performed on the primary image
        self._history = dict()

        # Check whether background has already been subtracted from this image
        # Try to read the 'BACK_SUB' key from the header and set self._subtracted accordingly
        self._history["sky-subtracted"] = True if self._findinheader("BACK_SUB") else False

        # Inform the user if this image has already been subtracted
        if self._history["sky-subtracted"]: self._log.success("Background already subtracted for this image")

    ## This function returns None if not found, False if found but value F, True if found and value T
    def _findinheader(self, key):

        try:
            # Get the value of this key
            value = self.header[key]

        # If the key doesn't exist
        except KeyError:

            # Set the value to None
            value = None

        # Return the value
        return value

    ## This function ..
    def _import(self, path, cut):

        # Show which image we are importing
        self._log.info("Reading in file: " + path)

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(path)

        # Get the primary image
        hdu = hdulist[0]

        # Get the image header
        self.header = hdu.header

        # Check for multiple planes with the PLANE keyword
        multiplanes = False
        if cut:
            try:

                plane_id = self.header['PLANE0']
                self._log.warning("Multiple planes detected. Using PLANE0 = " + plane_id)
                multiplanes = True

            # If the PLANE keyword is not found, we assume we only have one plane
            except KeyError: pass

        # Get the data from the HDU
        data = hdu.data[0] if multiplanes else hdu.data

        # Close the fits file
        hdulist.close()

        # Save the data in an image layer
        self._layers['primary'] = ImageLayer(data, self._log)

    # -----------------------------------------------------------------

    ## This function ..
    def info(self):

        # Print dimension of data
        self._log.info("Dimensions of data cube: " + str(self.xsize) + " x " + str(self.ysize))

        # Print type of data
        self._log.info("Type of data: " + str(self.datatype))

    @property
    def datatype(self):

        return self._layers['primary'].datatype()

    ## This function returns the sky layer of this image
    @property
    def sky(self):

        return self._layers['sky']

    ## This function
    @property
    def primary(self):

        return self._layers['primary']

    ## This function
    @primary.setter
    def primary(self, newdata):

        self._layers['primary'] = newdata

    ## This function
    @property
    def error(self):

        return self._layers['error']

    ## This function
    @property
    def totalmask(self):

        return self._layers['totalmask']

    # ----------------------------------------------------------------- BASIC PROPERTIES

    ## This function returns the number of pixels in the x direction
    @property
    def xsize(self):

        return self.primary.xsize

    ## This function returns the number of pixels in the y direction
    @property
    def ysize(self):

        return self.primary.ysize

    ## This function ..
    @property
    def mean(self):

        return self.primary.mean

    ## This function ..
    @property
    def median(self):

        return self.primary.median

    ## This function ..
    @property
    def min(self):

        return self.primary.min

    ## This function ..
    @property
    def max(self):

        return self.primary.max

    ## This function
    @property
    def stdev(self):

        return self.primary.stdev

    # ----------------------------------------------------------------- FILE OPERATIONS

    ## This function saves the image at the specified path
    def save(self, layer, path):

        # Create the HDU
        hdu = pyfits.PrimaryHDU(self._layers[layer].data, self.header)

        # Write to file
        hdu.writeto(path, clobber=True)

    # -----------------------------------------------------------------

    ## This function
    def _newaction(self, actionname):

        # Indicate that this action has been performed
        self._history[actionname] = True

    # ----------------------------------------------------------------- VISUALISATION

    ## This function
    def plot(self, path=None):

        # Plot the primary image
        self.primary.plot()

    ## This function
    def histogram(self, path=None):

        self.primary.histogram()

    # -----------------------------------------------------------------

    ## This function masks ...
    def mask(self, regions, interpolate=False):

        # New action
        self._newaction("mask")

        # Create the mask
        mask = regions.get_mask(shape=(self.ysize,self.xsize))

        # Mask the image with zeroes
        self.primary.data[mask] = 0

        # If we need to interpolate
        if interpolate:

            # Set some parameters
            order = 3
            linear = False

            # For each region, we interpolate within a box surrounding the region
            for region in regions:

                # Center and radius of this region
                xcenter = round(region.coord_list[0])
                ycenter = round(region.coord_list[1])
                radius  = round(region.coord_list[2])

                # Construct a box around this region
                xmin = int(xcenter -1 - 2*radius)
                xmax = int(xcenter -1 + 2*radius)
                ymin = int(ycenter -1 - 2*radius)
                ymax = int(ycenter -1 + 2*radius)
                xrange = slice(xmin,xmax)
                yrange = slice(ymin,ymax)
                box = self.primary.data[yrange, xrange]

                # Create mask and boxminusmask
                mask = np.zeros_like(box)
                boxminusmask = np.ones_like(box)
                mask[box == 0] = 1
                boxminusmask[box == 0] = 0

                # Interpolate
                range = ((xmin,xmax),(ymin,ymax))
                interpolatedbox = self._interpolate(order,linear,range)

                # Fill the data array with the interpolated values
                # * = elementwise product!
                self.primary.data[yrange, xrange] = self.primary.data[yrange, xrange]*boxminusmask + interpolatedbox*mask

        return mask

    ## This function ..
    def _interpolate(self, order, linear, limits=None):

        # Make lists to contain the coordinates and fluxes from which to interpolate from
        xvalues = []
        yvalues = []
        fluxes = []

        # Iterate over the image pixels
        for x in range(self.xsize):

            # Check whether this x value falls within the range
            if x < limits[0][0] or x > limits[0][1]: continue

            for y in range(self.ysize):

                # Check whether this y value falls within the range
                if y < limits[1][0] or y > limits[1][1]: continue

                # Get the flux in this pixel if it is not masked
                if self.primary.data[y,x] != 0:

                    xvalues.append(x)
                    yvalues.append(y)
                    fluxes.append(self.primary.data[y,x])

        # Make nummy arrays from x and y coordinates
        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)

        # Fit a 2D polynomial
        self._log.info("Fitting a 2D polynomial to the surrounding pixel values")
        m = fitpolynomial(xvalues, yvalues, fluxes, order, linear)

        # Evaluate the polynomial in the specified range
        xx, yy = np.meshgrid(range(limits[0][0],limits[0][1]), range(limits[1][0],limits[1][1]))
        box_fluxes = polynomial(xx, yy, m)

        # Return the result
        return box_fluxes

    ## This function performs aperture photometry on the region
    def photometry(self, region):

        # If the region is a circle
        if region.name == "circle":

            # Center and radius of this region
            xcenter = round(region.coord_list[0])
            ycenter = round(region.coord_list[1])
            radius  = round(region.coord_list[2])
            center = (xcenter, ycenter)

            # Create an aperture object
            aperture = CircularAperture(center, radius)

            # Calculate the photometry in the aperture
            fluxtable = aperture_photometry(self.primary.data, aperture)

            # Get the sum of the flux within the aperture and return it
            return fluxtable['aperture_sum'][0]

            # x and y coordinate of the center of the aperture
            #x = fluxtable['xcenter'][0][0]
            #y = fluxtable['ycenter'][0][0]

        # If the region is a box
        elif region.name == "box":

            # Get the center, width and height for this box region
            xcenter = region.coord_list[0]
            ycenter = region.coord_list[1]
            width = region.coord_list[2]
            height = region.coord_list[3]

            # Determine the coordinates of the corner pixels
            xmin = int(round(xcenter - 0.5*width))
            xmax = int(round(xcenter + 0.5*width))
            ymin = int(round(ycenter - 0.5*height))
            ymax = int(round(ycenter + 0.5*height))

            # Calculate the sum of the flux in the pixels within the box
            return np.sum(self.primary.data[ymin:ymax,xmin:xmax])

    # ----------------------------------------------------------------- BASIC OPERATIONS

    ## This function crops this image (consisting of all layers, masks and regions) to a specified field.
    #  It takes the following parameters:
    #
    #  - xrange: a slice object: slice(lowx, highx)
    #  - yrange: a slice object: slice(lowy, highy)
    #
    def crop(self, xrange, yrange):

        # Save the current state of the image
        self._newaction("cropped")

        # TODO: crop all layers, masks and edges!

        # Crop the image
        self.primary.data = self.primary.data[yrange, xrange]

    ## This function
    def rotate(self):

        # Nans should be masked !

        ## FROM FITSKIRT_ROTATE

        shift_x = inframe.shape[0]/2 - (left_x + right_x)/2
        shift_y = inframe.shape[1]/2 - (left_y + right_y)/2

        print "Shifted frame to center by " + str(shift_x) + "," + str(shift_y)

        shiftframe = ndimage.interpolation.shift(inframe,(shift_x, shift_y))
        rotangle = math.degrees(math.atan(float(left_y - right_y)/float(left_x - right_x)))
        rotangle += 180.0 if flip else 0.0

        rotframe = ndimage.interpolation.rotate(shiftframe,rotangle)

        self._log.info("Rotated frame over " + str(rotangle) + " degrees")

        ##

    ## This function
    def downsample(self, factor):

        # Use the zoom function to resample
        newdata = ndimage.interpolation.zoom(self.primary.data, zoom=1.0/factor)

        # Add the layer
        self.addlayer(newdata, 'primary_downsampled')

    # ----------------------------------------------------------------- VIEW AND ADD LAYERS, REGIONS AND MASKS

    ## This function
    def layers(self):

        # Return the keys of the layers dictionary
        return self._layers.keys()

    ## This function
    def regions(self):

        # Return the keys of the regions dictionary
        return self._regions.keys()

    ## This function
    def masks(self):

        # Return the keys of the masks dictionary
        return self._masks.keys()

    ## This function ..
    def addlayer(self, data, name):

        # Inform the user
        self._log.info("Adding '" + name + "' to the set of image layers")

        # Add the layer to the layers dictionary
        self._layers[name] = ImageLayer(data, self._log)

    ## This function
    def addregion(self, region, name):

        # Inform the user
        self._log.info("Adding '" + name + "' to the set of regions")

        # Add the region to the regions dictionary
        self._regions[name] = region

    ## This function ..
    def addmask(self, data, name):

        # Inform the user
        self._log.info("Adding '" + name + "' to the set of masks")

        # Add the mask to the masks dictionary
        self._masks[name] = ImageMask(data, self._log)

    # ----------------------------------------------------------------- ADVANCED OPERATIONS

    ## This function ..
    def convolve(self):

        pass

    # ----------------------------------------------------------------- MASKS

    ## This function
    def combinemasks(self, mask1, mask2, name=None):

        # Combine the 2 masks
        mask = self._masks[mask1].data + self._masks[mask2].data

        # If no name is given for this mask, combine the names of the two original masks
        if name is None: name = mask1 + "_and_" + mask2

        # Add the mask
        self.addmask(mask.astype(bool), name)

    ## This function ..
    def maketotalmask(self):

        # For now, we only have the edge mask
        totalmask = self._masks['edges'].data

        # Add the mask
        self.addmask(totalmask, 'total')

    ## This function ..
    def makeskymask(self):

        # Combine the edges and galaxy masks
        self.combinemasks("edges", "galaxy", name="sky")

    ## This function
    def maskedges(self):

        # Structure array
        structure = ndimage.generate_binary_structure(2, 2)

        # Make the new mask, made from 100 iterations with the structure array
        mask = ndimage.binary_dilation(self._masks["nans"].data, structure, iterations=100)

        # Add the mask
        self.addmask(mask, "edges")

    ## This function applies a mask on the primary image
    def makemaskedlayer(self, maskname):

        # Copy the primary image
        maskedprimary = np.copy(self.primary.data)

        # Mask this copy
        maskedprimary[self._masks[maskname].data] = 0

        # Add this masked image to the layers
        self.addlayer(maskedprimary, 'primary_masked_' + maskname)

    ## This function applies the total mask on the primary image
    def applytotalmask(self):

        # Set the pixel value to 0 where the mask is True
        self.primary.data[self._masks['total'].data] = 0

    ## This function masks the NaN values in the primary image
    def masknans(self):

        # Get a map of all the NaNs in the primary image
        mask = np.isnan(self.primary.data)

        # Make a nans mask layer
        self.addmask(mask, "nans")

    # This function creates a new mask from a specified region. It takes the name of this region as the sole argument.
    def createmask(self, regionname):

        # Get the region
        region = self._regions[regionname]

        # Create the mask
        mask = region.get_mask(header=self.header, shape=(self.ysize,self.xsize))

        # Add the mask to the masks list
        self.addmask(mask, regionname)

    # ----------------------------------------------------------------- IMAGE SEGMENTATION

    ## This function
    def findgalaxy(self, plot=True):

        # Find the galaxy
        galaxy = GalaxyFinder(self.primary.data[::-1,:], quiet=True)

        # Plot the ellips onto the image frame
        if plot: galaxy.plot()

        # Get the galaxy parameters
        #positionangle = np.radians(-galaxy.theta)
        positionangle = galaxy.theta
        xcenter = galaxy.xpeak
        ycenter = galaxy.ypeak
        eps = galaxy.eps

        # The length of the major axis of the ellipse
        mjr = 3.0 * galaxy.majoraxis

        # The width and heigth of the ellips
        width = mjr
        height = mjr * (1 - eps)

        # Cretae a string identifying this ellipse
        region_string = "image;ellipse(" + str(xcenter) + "," + str(ycenter) + "," + str(width) + "," + str(height) + "," + str(positionangle) + ")"

        # Create a region consisting of one ellipse
        region = pyregion.parse(region_string)

        # Add this region
        self.addregion(region, "galaxy")

    ## Automatically detect sources
    def detectsources(self):

        from astropy.stats import sigma_clipped_stats

        mean, median, std = sigma_clipped_stats(self.primary.data, sigma=3.0)
        print(mean, median, std)

        from photutils import daofind
        sources = daofind(self.primary.data, fwhm=3.0, threshold=2.*std)

        print sources

        from astropy.visualization import SqrtStretch

        from astropy.visualization.mpl_normalize import ImageNormalize

        positions = (sources['xcentroid'], sources['ycentroid'])

        apertures = CircularAperture(positions, r=4.)

        norm = ImageNormalize(stretch=SqrtStretch())

        plt.imshow(self.primary.data, cmap='Greys', origin='lower', norm=norm)

        apertures.plot(color='blue', lw=1.5, alpha=0.5)

        plt.show()

        #from astropy.visualization import SqrtStretch
        #from astropy.visualization.mpl_normalize import ImageNormalize

        #norm = ImageNormalize(stretch=SqrtStretch())
        #fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        #ax1.imshow(self.primary.data, origin='lower', cmap='Greys_r', norm=norm)
        #ax2.imshow(segm, origin='lower', cmap='jet')

        # Return the regions
        #return regions

    ## This function determines the central peak position of the stars indicated by the region file
    def getstarpositions(self, region, plot=True):

        # Make an empty list of stars
        stars = []

        # Loop over all the shapes in this region and fit the stellar profiles with a 2D Gaussian distribution
        for shape in self._regions[region]:

            # Initially, set the minimum and maximum x and y values to zero
            xmin = xmax = ymin = ymax = 0

            # If the region is a circle
            if shape.name == "circle":

                # Center and radius of this circle
                xcenter = round(shape.coord_list[0])
                ycenter = round(shape.coord_list[1])
                radius  = round(shape.coord_list[2])

                # Determine the coordinates of the circle's bounding box
                xmin = int(round(xcenter - radius))
                xmax = int(round(xcenter + radius))
                ymin = int(round(ycenter - radius))
                ymax = int(round(ycenter + radius))

            # If the region is a box
            elif shape.name == "box":

                # Get the center, width and height for this box
                xcenter = shape.coord_list[0]
                ycenter = shape.coord_list[1]
                width = shape.coord_list[2]
                height = shape.coord_list[3]

                # Determine the coordinates of the corner pixels
                xmin = int(round(xcenter - 0.5*width))
                xmax = int(round(xcenter + 0.5*width))
                ymin = int(round(ycenter - 0.5*height))
                ymax = int(round(ycenter + 0.5*height))

            # Cut out a square of the primary image around the star
            square = self.primary.data[ymin:ymax, xmin:xmax]

            # Fit a 2D Gaussian to the brightness distribution
            params = fitgaussian(square)
            fit = gaussian(*params)

            # Unpack the parameters
            (height, x, y, width_x, width_y) = params

            # Plot the result
            if plot:

                plt.matshow(square, cmap=cm.CMRmap)
                plt.contour(fit(*indices(square.shape)), cmap=cm.Blues)
                ax = plt.gca()

                plt.text(0.95, 0.05, """
                        pxl max: %.1f
                        mod max: %.1f
                        x : %.1f
                        y : %.1f
                        width_x : %.1f
                        width_y : %.1f""" %(square[round(x),round(y)],height, x, y, width_x, width_y),
                            fontsize=16, horizontalalignment='right',color='white',
                            verticalalignment='bottom', transform=ax.transAxes)

                plt.show()

            # Add the fitted parameters to the list of stars.
            # NOTE: for some reason, x and y are interchanged by python.
            x = xmin + y
            y = ymin + x

            # Add the parameters of this star to the stars list
            stars.append((x,y,width_x,width_y))

        # Return the list of stars (their positions)
        return stars

    # ----------------------------------------------------------------- SKY-SUBTRACTION

    ## This function estimes the sky from a region object
    #
    def estimatesky_fitskirt(self, region, order=3, linear=True):

        # Register this new action
        self._newaction("sky-subtraction")

        # For each shape in the region, we calculate the mean flux
        xvalues = []
        yvalues = []
        fluxes = []

        # For each shape in the specified region
        for shape in self._regions[region]:

            # Calculate the mean flux within this shape
            flux = self.photometry(shape) / area(shape)

            # Append the x coordinate, y coordinate and flux to the appropriate lists
            xvalues.append(shape.coord_list[0])
            yvalues.append(shape.coord_list[1])
            fluxes.append(flux)

        # Calculate and log the mean sky value
        meansky = np.mean(fluxes)
        self._log.info("Mean sky value = " + str(meansky))

        # Make numpy arrays
        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)

        # Interpolate the fluxes
        parameters = fitpolynomial(xvalues, yvalues, fluxes, order, linear)

        # Make a sky array and subtract it from the object frame
        sky = np.zeros_like(self.primary.data)

        # Create a row of pixels for the sky
        strip = np.arange(np.float(self.xsize))

        # Calculate the sky values on the grid
        for y in range(0, self.ysize):

            # Get the values of the polynomial
            sky[y,:] = polynomial(strip, y, parameters)

        # Add the sky as a new layer to this image
        self.addlayer(sky, "sky")

    ## This function
    def estimatesky(self):

        # Get the background
        bkg = Background(self.primary.data, (100, 100), filter_shape=(3, 3), method='median', mask=self._masks['sky'].data)
        sky = bkg.background

        # Add a new layer
        self.addlayer(sky, 'sky')

    ## This function subtracts the fitted sky from the primary image
    def subtractsky(self):

        # Determine the negative of the total mask
        negativetotalmask = np.logical_not(self._masks["total"].data)

        # Subtract the sky from the data
        self.primary.data = self.primary.data - self.sky.data*negativetotalmask

    # ----------------------------------------------------------------- PSF DETERMINATION

    ## This function estimates the psf
    def estimatepsf_fitskirt(self, region):

        # Get the stars
        stars = self.getstarpositions(region)

        # Initially, set the average x and y position to
        FWHM_x = FWHM_y = 0.0

        # For each star in the list
        for star in stars:

            # Get the width in the x and y direction
            FWHM_x += star[2]
            FWHM_y += star[3]

        # Average the x and y FWHM over all the reference stars
        FWHM_x = FWHM_x / float(len(stars))
        FWHM_y = FWHM_y / float(len(stars))

        # Circular approximation
        FWHM = (FWHM_x + FWHM_y) / 2.0

        # Inform the user
        self._log.info("FWHM : " + str("{:5.2f}".format(FWHM)) + "  ("+ str("{:5.2f}".format(FWHM_x)) + " in x direction and " + str("{:5.2f}".format(FWHM_y)) + " in y direction )")

        #return aver

    ## This function
    def estimatepsf(self):

        pass

# ----------------------------------------------------------------- USEFUL FUNCTIONS

## This function returns the area for a certain region
def area(region):

    # If this region is a circle
    if region.name == "circle":

        return math.pi * region.coord_list[2] * region.coord_list[2]

    # If this region is a box
    if region.name == "box":

        return region.coord_list[2] * region.coord_list[3]

## This function ..
def plotdata(data, path):

    # Plot the data using logaritmic scale
    plt.imshow(data, cmap='gray', norm=LogNorm(), interpolation='nearest')

    # Add a color bar
    plt.colorbar()

    # Display the result
    plt.show()

# -----------------------------------------------------------------

# CLASS MASK
class ImageMask(object):

    ## This function
    def __init__(self, data, log):

        self._data = data

        self._log = log

    ## This function ..
    @property
    def data(self):

        return self._data

    ## This function ..
    @data.setter
    def data(self, newdata):

        self._data = newdata

    ## This function ..
    def plot(self, path=None):

        plotdata(self._data.astype(int), path)

# -----------------------------------------------------------------

# CLASS IMAGELAYER
class ImageLayer(object):

    ## The constructor
    def __init__(self, data, log):

        # Copy the data
        self._data = data

        # Logger
        self._log = log

    ## This function
    @property
    def data(self):

        return self._data

    ## This function
    @data.setter
    def data(self, newdata):

        self._data = newdata

    ## This function
    def plot(self, path=None):

        plotdata(self._data, path)

    ## This function
    def histogram(self, path=None):

        NBINS = 1000
        plt.hist(self._data.flat, NBINS)

        # Display the result
        plt.show()

    ## This function
    @property
    def xsize(self):

        return self._data.shape[1]

    ## This function
    @property
    def ysize(self):

        return self._data.shape[0]

    ## This function
    def datatype(self):

        return self._data.dtype.name

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
    def stdev(self):

        # Set the delta degrees of freedom
        ddof = 1

        # Return the standard deviation of the data
        return np.std(self._data, ddof=ddof)

    ## This function
    def undo(self):

        # Check whether the previous state has been saved or not
        if self._prevdata is not None:

            # Replace the data with the previous data
            self._data = self._prevdata

            # Set the previous data to None
            self._prevdata = None

        else:

            self._log.warning("Cannot undo")

    ## This function ..
    def backup(self):

        self._prevdata = np.copy(self._data)
