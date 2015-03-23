#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.image This module includes classes used for working with astronomical FITS images.
#  The Image class in this package represents such a FITS image, and can be created by a statement as simple as:
#      im = Image("example.fits")
#  There are numerous operations possible on the image once it has been created with the above statement.

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
from pts.skirtunits import SkirtUnits

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
        try:
            self._filter = Filter(filtername)
        except ValueError:
            self._log.warning("Could not determine the filter used for this image")

        # Units
        self._units = self._findinheader("BUNIT")

        # The FWHM of the PSF
        self._fwhm = None

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
        if multiplanes:

            self.addlayer(hdu.data[0], "primary")
            self.addlayer(hdu.data[1], "errors")

            # Close the fits file
            hdulist.close()

        else:

            self.addlayer(hdu.data, "primary")

            # Close the fits file
            hdulist.close()

            # Look for a seperate errors fits file
            errorspath = os.path.splitext(path)[0] + ".ERRORS.fits"
            if os.path.isfile(errorspath):

                # Open the HDU list for the FITS file
                hdulist = pyfits.open(errorspath)

                # Get the primary data
                hdu = hdulist[0]

                self.addlayer(hdu.data, "errors")

                # Close the fits file
                hdulist.close()

            else:

                self._log.warning("No error data found for this image")

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
    def plot(self, layer):

        # Plot the specified layer
        self._layers[layer].plot()

    ## This function
    def histogram(self, layer):

        # Make a histogram of the specified layer
        self._layers[layer].histogram()

    ## This function
    def contourplot(self, layer):

        # Make a contour plot of the specified layer
        self._layers[layer].contourplot()

    # ----------------------------------------------------------------- PHOTOMETRY

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
    def multiply(self, factor):

        # Multiply the image with the factor
        self.primary.data = self.primary.data * factor

    ## This function
    def rotateandcenter_fitskirt(self, left_x, left_y, right_x, right_y, flip=False):

        # Nans should be masked !

        shift_x = self.ysize/2 - (left_x + right_x)/2
        shift_y = self.xsize/2 - (left_y + right_y)/2

        print "Shifted frame to center by " + str(shift_x) + "," + str(shift_y)

        shiftframe = ndimage.interpolation.shift(self.primary.data,(shift_x, shift_y))
        angle = math.degrees(math.atan(float(left_y - right_y)/float(left_x - right_x)))
        angle += 180.0 if flip else 0.0

        # Create the rotated frame
        rotframe = ndimage.interpolation.rotate(shiftframe, angle)

        # TODO: rotate the other layers, regions and masks!

        # Inform the user of the rotation angle
        self._log.info("Rotated frame over " + str(angle) + " degrees")

        # Add the new, rotated layer
        self.addlayer(rotframe, "primary_rotated")

    ## This function can be used to rotate the frame over a given angle
    def rotate(self, angle):

        # Create the rotated frame
        rotframe = ndimage.interpolation.rotate(self.primary.data, angle)

        # TODO: rotate the other layers, regions and masks!

        # Inform the user of the rotation angle
        self._log.info("Rotated frame over " + str(angle) + " degrees")

        # Add the new, rotated layer
        self.addlayer(rotframe, "primary_rotated")

    ## This function rotates the image so that the galactic plane lies horizontal
    def autorotate(self):

        # Rotate about the position angle of the galaxy
        self.rotate(-self.orientation.theta)

    ## This function makes a new layer where the center of the galaxy is in the center of the plane
    def autocenter(self):

        # Determine the center pixel of the image
        imagecenter_x = self.xsize / 2.0
        imagecenter_y = self.ysize / 2.0

        # Calculate the shift to be made in the x and y directions
        shift_x = imagecenter_x - self.orientation.xpeak
        shift_y = imagecenter_y - self.orientation.ypeak

        # Create a centered frame
        centered = ndimage.interpolation.shift(self.primary.data,(shift_y, shift_x))

        # TODO: center the other layers, regions and masks!

        # Add the new, centered layer
        self.addlayer(centered, "primary_centered")

    ## This function
    def downsample(self, factor):

        # Use the zoom function to resample
        newdata = ndimage.interpolation.zoom(self.primary.data, zoom=1.0/factor)

        # Add the layer
        self.addlayer(newdata, 'primary_downsampled')

    # ----------------------------------------------------------------- UNITS

    ## This function can be used to set the units for this image
    def setunits(self, units):

        # Make a SkirtUnits object
        self._units = SkirtUnits("extragalactic", "frequency")

    ## This function
    def convert(self, units):

        # Calculate the conversion factor
        conversionfactor = self._units.convert(1.0, units)

        # Convert the data
        self.primary.data *= conversionfactor

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
    def convolve(self, name):

        # kernels: from http://www.astro.princeton.edu/~ganiano/Kernels/Ker_2012_May/Kernels_fits_Files/Hi_Resolution/

        self._log.info("Convolving image with the kernel " + os.path.splitext(name)[0])

        # The path to the kernel file
        path = os.path.join(os.getenv("HOME"), "Kernels", name)

        # Inform the user that the kernel was found
        self._log.info("Found kernel file at " + path)

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(path)

        # Get the primary image
        hdu = hdulist[0]

        # Do the convolution
        from astropy.convolution import convolve_fft
        resultdata = convolve_fft(self.primary.data, hdu.data)

        # Close the FITS file
        hdulist.close()

        # Add the new convolved layer
        self.addlayer(resultdata, "primary_convolved")

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

    ## This function interpolates the image within the shapes in a certain region ...
    def interpolate(self, r, order=3, linear=False):

        # Make a new copy of the primary image
        interpolated = np.copy(self.primary.data)

        # For each region, we interpolate within a box surrounding the region
        for shape in self._regions[r]:

            # Center and radius of this region
            xcenter = round(shape.coord_list[0])
            ycenter = round(shape.coord_list[1])
            radius  = round(shape.coord_list[2])

            # Construct a box around this region
            xmin = int(xcenter -1 - 2*radius)
            xmax = int(xcenter -1 + 2*radius)
            ymin = int(ycenter -1 - 2*radius)
            ymax = int(ycenter -1 + 2*radius)
            xrange = slice(xmin,xmax)
            yrange = slice(ymin,ymax)
            box = self.primary.data[yrange, xrange]

            # Create mask and boxminusmask
            mask = np.zeros_like(box, dtype=float)
            boxminusmask = np.ones_like(box, dtype=float)
            mask[box == 0] = 1.0
            boxminusmask[box == 0] = 0

            # Interpolate
            range = ((xmin,xmax),(ymin,ymax))
            interpolatedbox = self._interpolate(order,linear,range)

            # Fill the data array with the interpolated values
            # * = elementwise product!
            interpolated[yrange, xrange] = self.primary.data[yrange, xrange]*boxminusmask + interpolatedbox*mask

        # Add the new layer
        self.addlayer(interpolated, "primary_interpolated")

    # ----------------------------------------------------------------- MASKS

     ## This function masks the NaN values in the primary image
    def masknans(self):

        # Get a map of all the NaNs in the primary image
        mask = np.isnan(self.primary.data)

        # Make a nans mask layer
        self.addmask(mask, "nans")

    ## This function
    def maskedges(self):

        # Structure array
        structure = ndimage.generate_binary_structure(2, 2)

        # Make the new mask, made from 100 iterations with the structure array
        mask = ndimage.binary_dilation(self._masks["nans"].data, structure, iterations=100)

        # Add the mask
        self.addmask(mask, "edges")

    ## This function
    def combinemasks(self, m1, m2, name=None):

        # Combine the 2 masks
        mask = self._masks[m1].data + self._masks[m2].data

        # If no name is given for this mask, combine the names of the two original masks
        if name is None: name = m1 + "_and_" + m2

        # Add the mask
        self.addmask(mask.astype(bool), name)

    ## This function ..
    def maketotalmask(self):

        # Create a total mask
        totalmask = np.zeros_like(self.primary.data, dtype=bool)

        for mask in self.masks():

            totalmask += self._masks[mask].data

        # Add the mask
        self.addmask(totalmask, 'total')

    ## This function applies a mask on the primary image
    def makemaskedlayer(self, m):

        # Copy the primary image
        maskedprimary = np.copy(self.primary.data)

        # Mask this copy
        maskedprimary[self._masks[m].data] = 0

        # Add this masked image to the layers
        self.addlayer(maskedprimary, 'primary_masked_' + m)

    ## This function aplies a certain mask on the primary image
    def applymask(self, m):

        # Set the pixel value to 0 where the mask is True
        self.primary.data[self._masks[m].data] = 0

    ## This function applies the total mask on the primary image
    def applytotalmask(self):

        # Set the pixel value to 0 where the mask is True
        self.primary.data[self._masks['total'].data] = 0

    # This function creates a new mask from a specified region. It takes the name of this region as the sole argument.
    def createmask(self, r):

        # Get the region
        region = self._regions[r]

        # Create the mask
        mask = region.get_mask(header=self.header, shape=(self.ysize,self.xsize))

        # Add the mask to the masks list
        self.addmask(mask, r)

    # ----------------------------------------------------------------- IMAGE SEGMENTATION

    ## This function
    def findgalaxy(self, plot=False):

        # Find the orientation of the galaxy in this iamge
        self.orientation = GalaxyFinder(self.primary.data[::-1,:], quiet=True)

        # Plot the ellips onto the image frame
        if plot: self.orientation.plot()

        # The length of the major axis of the ellipse
        major = 3.0 * self.orientation.majoraxis

        # The width and heigth of the ellips
        width = major
        height = major * (1 - self.orientation.eps)

        # Cretae a string identifying this ellipse
        region_string = "image;ellipse(" + str(self.orientation.xpeak) + "," + str(self.orientation.ypeak) + "," + str(width) + "," + str(height) + "," + str(self.orientation.theta) + ")"

        # Create a region consisting of one ellipse
        region = pyregion.parse(region_string)

        # Add this region
        self.addregion(region, "galaxy")

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
    #def estimatesky(self):

        # Get the background
        #bkg = Background(self.primary.data, (100, 100), filter_shape=(3, 3), method='median', mask=self._masks['sky'].data)
        #sky = bkg.background

        # Add a new layer
        #self.addlayer(sky, 'sky')

    ## This function
    def estimatesky(self):

        # The length of the major axis of the ellipse
        major = 3.0 * self.orientation.majoraxis * 2.0

        # The width and heigth of the ellips
        width = major
        height = major * (1 - self.orientation.eps)

        # Create a string identifying this ellipse
        region_string = "image;ellipse(" + str(self.orientation.xpeak) + "," + str(self.orientation.ypeak) + "," + str(width) + "," + str(height) + "," + str(self.orientation.theta) + ")"

        # Create a region consisting of one ellipse
        region = pyregion.parse(region_string)

        # Create the mask
        mask = np.logical_not(region.get_mask(header=self.header, shape=(self.ysize,self.xsize)))

        # Combine the new mask with the galaxy
        newmask = mask + self._masks["galaxy"].data + self._masks["total"].data

        # Add the mask
        self.addmask(newmask.astype(bool), "sky")

        # Make a masked layer, the background
        self.makemaskedlayer("sky")

        # Create a NumPy masked array
        maskedarray = np.ma.masked_array(self.primary.data, newmask.astype(bool))

        # Calculate the mean, median and standard deviation of the sky around the galaxy
        mean = np.ma.mean(maskedarray)
        median = np.ma.median(maskedarray)
        error = np.ma.std(maskedarray, ddof=1)

        # Return the mean, median and standard deviation
        return mean, median, error

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

        # Initially, set the average x and y fwhm to zero
        fwhm_x = fwhm_y = 0.0

        # For each star in the list
        for star in stars:

            # Get the width in the x and y direction
            fwhm_x += star[2]
            fwhm_y += star[3]

        # Average the x and y FWHM over all the reference stars
        fwhm_x = fwhm_x / float(len(stars))
        fwhm_y = fwhm_y / float(len(stars))

        # Return the fwhm in the x and the y direction
        return fwhm_x, fwhm_y

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
    def contourplot(self, path=None):

        # Make the contours
        plt.contour(xi,yi,zi,15, linewidths=0.5, colors='k')
        plt.contourf(xi,yi,zi,15, cmap=plt.cm.jet)

        # Add a color bar
        plt.colorbar()

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
