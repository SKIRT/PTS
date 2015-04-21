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

# Import astronomical modules
import aplpy
import pyregion
import astropy.io.fits as pyfits
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils import CircularAperture
from photutils import aperture_photometry

# Import relevant PTS modules
from pts.inpaint import replace_nans
from pts.mathematics import fitpolynomial, polynomial, fitgaussian, gaussian
from pts.log import Log
from pts.filter import Filter
from pts.galaxy import GalaxyFinder
from pts.skirtunits import SkirtUnits

# -----------------------------------------------------------------

# Disable astropy logging except for warnings and errors
from astropy import log
log.setLevel("WARNING")

# -----------------------------------------------------------------

# Do not show warnings, to block Canopy's UserWarnings from spoiling the console log
import warnings
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------

## An instance of the Image class represents a FITS image file
#
class Image(object):

    ## The constructor accepts the following arguments:
    #
    #  - path: the path of the FITS image file
    #  - log: a Log instance
    #
    def __init__(self, filename, log=None):

        # Create a logger or use a pre-existing logger if possible
        self._log = Log() if log is None else log

        # Check if the specified file exists, otherwise exit with an error
        if not os.path.isfile(filename):

            self._log.error("No such file: " + filename)
            exit()

        # Set the name of the image
        self.name = os.path.splitext(os.path.basename(filename))[0]

        # Frames: primary, sky, fittedsky, errors ...
        self.frames = layers()

        # Masks: nans, edges, galaxy, stars
        self.masks = layers()

        # Regions: sky, stars,
        self.regions = layers()

        # Read in the image
        self._loadimage(filename)

        # The FWHM of the PSF
        self.fwhm = None

        # The entries in this dictionary indicate whether certain operations have been performed on the primary image
        self._history = dict()

    ## This function ...
    def _loadimage(self, filename):

        # Show which image we are importing
        self._log.info("Reading in file: " + filename)

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(filename)

        # Get the primary HDU
        hdu = hdulist[0]

        # Get the image header
        header = hdu.header

        # Load the frames
        self.pixelscale = getpixelscale(header)

        # Obtain the filter for this image
        self.filter = getfilter(self.name, header)

        # Obtain the units of this image
        self.units = getunits(header)

        # Check whether multiple planes are present in the FITS image
        nframes = getnumberofframes(header)
        if nframes > 1:

            # For each frame
            for i in range(nframes):

                # Get the name of this frame, but the first frame always gets the name 'primary'
                description = getframedescription(header, i) if i else "the primary signal map"
                name = getframename(description) if i else "primary"

                # Add this frame to the frames dictionary
                self._addframe(hdu.data[i], name, description)

        else:

            # Add the primary image frame
            self._addframe(hdu.data, "primary", "the primary signal map")

        # Set the basic header for this image
        self.header = header.copy(strip=True)
        self.header["NAXIS"] = 2
        self.header["NAXIS1"] = self.xsize
        self.header["NAXIS2"] = self.ysize

        # Select the primary image frame
        self.frames.primary.select()

        # Close the FITS file
        hdulist.close()

    # ----------------------------------------------------------------- PRINT INFORMATION

    ## This function ...
    def status(self):

        # Show white line
        self._log.info("")

        # List the different frames
        self._log.info("Frames:")
        self._log.info("------------------------------")
        self.frames.list()

        # List the different regions
        self._log.info("")
        self._log.info("Regions:")
        self._log.info("------------------------------")
        self.regions.list()

        # List the different masks
        self._log.info("")
        self._log.info("Masks:")
        self._log.info("------------------------------")
        self.masks.list()

        # Show white line
        self._log.info("")

    ## This function ...
    def info(self):

        # Print name
        self._log.info("Name: " + self.name)

        # Print dimension of data
        self._log.info("Dimensions of data array: " + str(self.xsize) + " x " + str(self.ysize))

        # Print type of data
        self._log.info("Type of data: " + str(self.dtype))

    # ----------------------------------------------------------------- BASIC PROPERTIES OF THE PRIMARY IMAGE

    ## This function returns the data type for the primary image
    @property
    def dtype(self):

        return self.frames.primary.dtype

    ## This function returns the number of pixels in the x direction
    @property
    def xsize(self):

        return self.frames.primary.xsize

    ## This function returns the number of pixels in the y direction
    @property
    def ysize(self):

        return self.frames.primary.ysize

    ## This function ...
    @property
    def mean(self):

        return self.frames.primary.mean

    ## This function ...
    @property
    def median(self):

        return self.frames.primary.median

    ## This function ...
    @property
    def min(self):

        return self.frames.primary.min

    ## This function ...
    @property
    def max(self):

        return self.frames.primary.max

    ## This function
    @property
    def stdev(self):

        return self.frames.primary.stdev

    # ----------------------------------------------------------------- FILE OPERATIONS

    ## This function saves the currently active frame
    def save(self, path):

        # Find the active layer
        frame = self.frames.getactive()[0]

        # Inform the user
        self._log.info("Saving " + frame + " frame as " + path)

        # Create the HDU
        hdu = pyfits.PrimaryHDU(self.frames[frame].data, self.header)

        # Write to file
        hdu.writeto(path, clobber=True)

    ## This function is used to import the (first) frame of another FITS file into this image
    def importframe(self, path, name):

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(path)

        # Get the primary data
        hdu = hdulist[0]

        # Add the error frame
        self._addframe(hdu.data, name)

        # Close the fits file
        hdulist.close()

    ## This function is used to import a new ImageRegion from a regions file
    def importregion(self, path, name):

        # Create an pyregion object from the regions file
        region = pyregion.open(path)

        # Add the region to the set of regions
        self._addregion(region, name)

    ## This function is used to export the selected region(s) to a region file
    def exportregion(self, path):

        # Find the active region
        region = self.regions.getactive()[0]

        # Inform the user
        self._log.info("Creating " + path + " from " + region)

        # Write the region file
        self.regions[region]._region.write(path)

    # -----------------------------------------------------------------

    ## This function
    def _newaction(self, actionname):

        # Indicate that this action has been performed
        self._history[actionname] = True

    # ----------------------------------------------------------------- VISUALISATION

    ## This function shows a plot of the currently active frame, combined with the active regions and masks
    def plot(self, path = None, color=True, grid=False, blacknan=False, publication=False):

        # Get the currently active frame
        frame = self.frames.getactive()[0]

        # Create a total mask of the currently active masks
        totalmask = self.combinemasks()

        # Mask the frame with nans
        maskedimage = np.ma.array(self.frames[frame].data, mask = totalmask)
        image_with_nans =  maskedimage.filled(np.NaN)

        # Create a HDU from this frame with the image header
        hdu = pyfits.PrimaryHDU(image_with_nans, self.header)

        if path is None:

            # Create a figure canvas
            figure = plt.figure(figsize=(12, 12))

            # Create a figure from this frame
            plot = aplpy.FITSFigure(hdu, figure=figure)

        else:

            # Create a figure from this frame
            plot = aplpy.FITSFigure(hdu)

        if color:

            # Plot in color scale
            plot.show_colorscale()

        else:

            # Plot in gray scale
            plot.show_grayscale()

        # Add a color bar
        plot.add_colorbar()

        if blacknan:

            # Set the nan color to black
            plot.set_nan_color('black')

        if grid:

            # Add a grid
            plot.add_grid()

        if publication:

            # Use the 'publication' theme
            plot.set_theme('publication')

        # Add the regions
        for region in self.regions.getactive():

            # Get the shape list
            shapes = self.regions[region]._region.as_imagecoord(self.header)

            # Add these shapes to the plot
            plot.show_regions(shapes)

        if path is None:

            plt.show()
            plt.close('all')

        else:

            plot.save(path)

    ## This function
    def histogram(self, layer):

        # Make a histogram of the specified layer
        self.frames[layer].histogram()

    ## This function
    def contourplot(self, layer):

        # Make a contour plot of the specified layer
        self.frames[layer].contourplot()

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
            fluxtable = aperture_photometry(self.frames.primary.data, aperture)

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
            return np.sum(self.frames.primary.data[ymin:ymax,xmin:xmax])


    # ----------------------------------------------------------------- ARITHMETIC OPERATIONS

    ## This function multiplies all active layers by a given factor
    def multiply(self, factor):

        # Fore each active frame
        for frame in self.frames.getactive():

            # Inform the user
            self._log.info("Multiplying " + frame + " frame with a factor of " + str(factor))

            # Multiply this frame with the given factor
            self.frames[frame].data = self.frames[frame].data * factor

    # ----------------------------------------------------------------- BASIC IMAGE MANIPULATION

    ## This function crops the currently active frames of this image to a specified field.
    #  It takes the following parameters:
    #
    #  - xrange: a slice object: slice(lowx, highx)
    #  - yrange: a slice object: slice(lowy, highy)
    #
    def crop(self, xrange, yrange):

        # For each active frame
        for frame in self.frames.getactive():

            # Inform the user that this frame is being cropped
            self._log.info("Cropping " + frame + " frame")

            # Crop this frame
            self.frames[frame].data = self.frames[frame].data[yrange, xrange]

    ## This function
    def rotateandcenter_fitskirt(self, left_x, left_y, right_x, right_y, flip=False):

        # Nans should be masked !

        shift_x = self.ysize/2 - (left_x + right_x)/2
        shift_y = self.xsize/2 - (left_y + right_y)/2

        print "Shifted frame to center by " + str(shift_x) + "," + str(shift_y)

        shiftframe = ndimage.interpolation.shift(self.frames.primary.data,(shift_x, shift_y))
        angle = math.degrees(math.atan(float(left_y - right_y)/float(left_x - right_x)))
        angle += 180.0 if flip else 0.0

        # Create the rotated frame
        rotframe = ndimage.interpolation.rotate(shiftframe, angle)

        # TODO: rotate the other layers, regions and masks!

        # Inform the user of the rotation angle
        self._log.info("Rotated frame over " + str(angle) + " degrees")

        # Add the new, rotated layer
        self._addframe(rotframe, "primary_rotated")

    ## This function rotates the currently selected frames over a given angle (in degrees)
    def rotate(self, angle):

        # For each active frame
        for frame in self.frames.getactive():

            # Inform the user that this frame is being rotated
            self._log.info("Rotating frame over " + str(angle) + " degrees")

            # Rotate this frame
            self.frames[frame].data = ndimage.interpolation.rotate(self.frames[frame].data, angle)

    ## This function rotates currently active frames so that the galactic plane lies horizontal
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
        centered = ndimage.interpolation.shift(self.frames.primary.data,(shift_y, shift_x))

        # TODO: center the other layers, regions and masks!

        # Add the new, centered layer
        self._addframe(centered, "primary_centered")

    ## This function downsamples the currently active frames by a specified zooming factor
    def downsample(self, factor):

        # For each active frame
        for frame in self.frames.getactive():

            # Inform the user
            self._log.info("Downsampling " + frame + " frame by a factor of " + str(factor))

            # Use the zoom function to resample
            self.frames[frame].data = ndimage.interpolation.zoom(self.frames[frame].data, zoom=1.0/factor)

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
        self.frames.primary.data *= conversionfactor

    ## This function converts the currently active frame(s) into magnitude scale, using the specified zero-point
    #  magnitude
    def to_magnitudes(self, m0):

        # For each active frame
        for frame in self.frames.getactive():

            # Inform the user
            self._log.info("Converting " + frame + " frame to magnitude scale")

            # Convert to magnitude scale
            self.frames[frame].data = m0 - 2.5 * np.log10(self.frames[frame].data)

    ## This function converts the currently active frame(s) from magnitude to flux scale, using the specified
    #  zero-point flux
    def to_fluxes(self, F0):

        # For each active frame
        for frame in self.frames.getactive():

            # Inform the user
            self._log.info("Converting " + frame + " frame to flux scale")

            # Convert to flux scale
            self.frames[frame].data = F0 * np.power(10.0, - self.frames[frame].data / 2.5)

    ## This function sets the orientation of the galaxy in this image
    def setorientation(self, orientation):

        self.orientation = orientation

    ## This functions sets the FWHM of the PSF for this image
    def setfwhm(self, fwhm):

        # Inform the user
        self._log.info("Setting the FWHM of the PSF for this image to " + str(fwhm) + " pixels")

        # Set the FWHM
        self.fwhm = fwhm

    ## This function sets the pixel scale for this image (in arcseconds)
    def setpixelscale(self, pixelscale):

        # Inform the user
        self._log.info("Setting the pixel scale of this image to " + str(pixelscale) + " arcseconds")

        # Set the pixel scale
        self.pixelscale = pixelscale

    # ----------------------------------------------------------------- ADD AND REMOVE LAYERS, REGIONS AND MASKS

    ## This function removes the currently selected frame(s)
    def deleteframes(self):

        # For each active frame
        for frame in self.frames.getactive():

            # Remove this frame from the frames dictionary
            del self.frames[frame]

    ## This function removes the currently selected region(s)
    def deleteregions(self):

        # For each active region
        for region in self.regions.getactive():

            # Remove this region from the regions dictionary
            del self.regions[region]

    ## This function removes the currently selected mask(s)
    def deletemasks(self):

        # For each active mask
        for mask in self.masks.getactive():

            # Remove this mask from the masks dictionary
            del self.masks[mask]

    ## This function ...
    def _addframe(self, data, name, description=None):

        # Inform the user
        self._log.info("Adding '" + name + "' to the set of image frames")

        # Add the layer to the layers dictionary
        self.frames[name] = ImageFrame(data, description, self._log)

    ## This function ...
    def _addregion(self, region, name):

        # Inform the user
        self._log.info("Adding '" + name + "' to the set of regions")

        # Add the region to the regions dictionary
        self.regions[name] = ImageRegion(region, self._log)

    ## This function ...
    def _addmask(self, data, name):

        # Inform the user
        self._log.info("Adding '" + name + "' to the set of masks")

        # Add the mask to the masks dictionary
        self.masks[name] = ImageMask(data, self._log)

    # ----------------------------------------------------------------- ADVANCED OPERATIONS

    ## This function convolves the currently selected frame(s) with a specified kernel (a FITS file)
    def convolve(self, name):

        # Import the convolution function
        from astropy.convolution import convolve_fft

        # The path to the kernel file
        path = os.path.join(os.getenv("HOME"), "Kernels", name)

        # Inform the user that the kernel was found
        self._log.info("Found kernel file at " + path)

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(path)

        # Get the kernel data and header
        kernel = hdulist[0].data
        header = hdulist[0].header

        # Inform the user
        self._log.info("Rebinning the kernel to the image pixel grid")

        # Get the pixel scale of the kernel
        pixelscale_kernel = header["CD1_1"]*3600

        # Calculate the zooming factor
        factor = self.pixelscale / pixelscale_kernel

        # Rebin the kernel to the same grid of the image
        kernel = ndimage.interpolation.zoom(kernel, zoom=1.0/factor)

        # Add a frame
        self._addframe(kernel, "kernel")

        # For all active frames, do the convolution
        for frame in self.frames.getactive():

            # Inform the user that this frame is being convolved
            self._log.info("Convolving " + frame + " frame with the kernel " + os.path.splitext(name)[0])

            # Do the convolution on this frame
            self.frames[frame].data = convolve_fft(self.frames[frame].data, kernel, normalize_kernel=True)

        # Close the FITS file
        hdulist.close()

    ## This function rebins the currently selected frame(s) based on a certain reference FITS file
    def rebin(self, reference):

        # Import the rebinning function
        from pts.hcongrid import hcongrid

        # Open the HDU list for the reference FITS file
        hdulist = pyfits.open(reference)

        # Get the primary image
        hdu = hdulist[0]

        referenceheader = hdu.header

        referenceheader["NAXIS"] = 2
        referenceheader.pop("NAXIS3", None)

        # For all active frames, do the rebinning
        for frame in self.frames.getactive():

            # Inform the user that this frame is being rebinned
            self._log.info("Rebinning " + frame + " frame to the grid of " + reference)

            # Do the rebinning based on the header of the reference image
            self.frames[frame].data = hcongrid(self.frames[frame].data, self.header, referenceheader)

        # Close the reference FITS file
        hdulist.close()

    ## This function interpolates the image within the combination of the currently active masks
    def interpolate(self):

        # Inform the user
        self._log.info("Interpolating the image in the areas covered by the currently selected masks")

        # Combine the active masks
        totalmask = self.combinemasks()

        # Make a copy of the image where masked pixels are filled with NaNs
        maskedimage = np.ma.array(self.frames.primary.data.astype(float), mask = totalmask)
        image_with_nans =  maskedimage.filled(np.NaN)

        # Interpolate the masked regions
        interpolated = replace_nans(image_with_nans, 5, 0.5, 2, "localmean")

        # Replace the primary image by the interpolated one
        self.frames.primary.data = interpolated

    # ----------------------------------------------------------------- MASKS

     ## This function masks the NaN values in the primary image
    def masknans(self):

        # Get a map of all the NaNs in the primary image
        mask = np.isnan(self.frames.primary.data)

        # Make a nans mask layer
        self._addmask(mask, "nans")

    ## This function
    def expandmasks(self, name, iterations=100):

        # Define the structure for the expansion
        structure = ndimage.generate_binary_structure(2, 2)

        # Get a combination of the active masks
        oldmask = self.combinemasks()

        # Make the new mask, made from 100 iterations with the structure array
        newmask = ndimage.binary_dilation(oldmask, structure, iterations)

        # Add the new, expanded mask
        self._addmask(newmask, name)

    ## This function ...
    def combinemasks(self, name=None):

        # Initialize an boolean array for the total mask
        totalmask = np.zeros_like(self.frames.primary.data, dtype=bool)

        # For each active mask
        for mask in self.masks.getactive():

            # Add this mask to the total
            totalmask += self.masks[mask].data

        # If no name is given for the total mask, return it
        if name is None:

            return totalmask

        # If a name is given, save the total mask under this name
        else:

            self._addmask(totalmask, name)

    ## This function makes a new mask which is the inverse (logical not) of the total currently selected mask
    def invertmask(self, name):

        # Get the total selected mask
        currentmask = self.combinemasks()

        # Calculate the inverse of the this total mask
        newmask = np.logical_not(currentmask)

        # Add the new, inverted mask
        self._addmask(newmask, name)

    ## This function applies the currently active masks to the primary image. Masked pixels are set to zero.
    def applymasks(self):

        # For each active mask
        for name in self.masks.getactive():

            # Set the corresponding image pixels to zero for this mask
            self.frames.primary.data[self.masks[name].data] = 0

    # This function creates a new mask from the currently selected region(s).
    def createmask(self):

        # Initialize an boolean array for the total mask
        totalmask = np.zeros_like(self.frames.primary.data, dtype=bool)

        name = ""

        # For each active region
        for region in self.regions.getactive():

            # Create the mask
            totalmask += self.regions[region]._region.get_mask(header=self.header, shape=(self.ysize,self.xsize))

            name += region + "_"

        # Remove the trailing underscore
        name = name.rstrip("_")

        # Add the mask to the masks list
        self._addmask(totalmask, name)

    # ----------------------------------------------------------------- IMAGE SEGMENTATION

    ## This function
    def findgalaxy(self, plot=False):

        # Find the orientation of the galaxy in this iamge
        self.orientation = GalaxyFinder(self.frames.primary.data[::-1,:], quiet=True)

        # Plot the ellips onto the image frame
        if plot: self.orientation.plot()

        # The length of the major axis of the ellipse
        major = 3.0 * self.orientation.majoraxis * 1.7

        # The width and heigth of the ellips
        width = major
        height = major * (1 - self.orientation.eps)

        # Create a string identifying this ellipse
        region_string = "image;ellipse(" + str(self.orientation.ypeak) + "," + str(self.orientation.xpeak) + "," + str(width) + "," + str(height) + "," + str(self.orientation.theta) + ")"

        # Create a region consisting of one ellipse
        region = pyregion.parse(region_string)

        # Add this region
        self._addregion(region, "galaxy")

    ## This function
    def fetchstars(self):

        from pts.catalog import fetch
        from astropy import wcs

        w = wcs.WCS(self.header)

        # Some pixel coordinates of interest.
        pixels = np.array([[0,0],[self.ysize-1,self.xsize-1]], dtype=float)

        # Convert pixel coordinates to world coordinates (RA and DEC in degrees)
        world = w.wcs_pix2world(pixels, 1)
        coordinate1 = world[0]
        coordinate2 = world[1]
        ra_range = [coordinate2[0], coordinate1[0]]
        dec_range = [coordinate1[1], coordinate2[1]]

        # Determine the center in RA and DEC (in degrees)
        ra_center = 0.5*(ra_range[0] + ra_range[1])
        dec_center = 0.5*(dec_range[0] + dec_range[1])

        # Determine the width in RA and DEC (both in degrees)
        ra_width = ra_range[1] - ra_range[0]
        dec_width = dec_range[1] - dec_range[0]

        # Determine the maximum width (in degrees)
        width = max(ra_width, dec_width)

        # Calculate the radius in arcminutes
        width_minutes = width * 60
        radius = 0.5*width_minutes

        # Fetch the stars, create a region
        region = fetch([ra_center,dec_center], radius, faint=16.0, filename="")

        # Add this region
        self._addregion(region, "stars")

    ## This function determines the central peak position of the stars indicated by the currently selected region file
    def getstarpositions(self, plot=False):

        # Get the currently selected region
        region = self.regions.getactive()[0]

        # Make an empty list of stars
        stars = []

        # Keep track of how many stars could not be fitted
        failed = 0

        # Loop over all the shapes in this region and fit the stellar profiles with a 2D Gaussian distribution
        for shape in self.regions[region]._region.as_imagecoord(self.header):

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

            # Check whether this shape falls within the image frame
            if xmin < 0 or xmax >= self.xsize or ymin < 0 or ymax >= self.ysize: continue

            # Cut out a square of the primary image around the star
            square = self.frames.primary.data[ymin:ymax, xmin:xmax]

            if not np.any(square):
                failed += 1
                continue

            try:
                # Fit a 2D Gaussian to the brightness distribution
                params = fitgaussian(square)
                fit = gaussian(*params)
            except:
                # If the fitting failed, skip this star
                failed += 1
                continue

            # Unpack the parameters
            (height, x, y, width_x, width_y) = params

            # Fix negative sigmas
            if width_x < 0: width_x = -width_x
            if width_y < 0: width_y = -width_y

            # Skip invalid parameter values, the fitting failed
            if np.isnan(width_x) or np.isnan(width_y):
                failed += 1
                continue

            # If the center of the Gaussian falls out of the square, skip this star
            if round(x) > square.shape[0] - 1 or round(y) > square.shape[1] - 1:
                failed += 1
                continue

            # Plot the result
            if plot:

                plt.matshow(square)
                plt.contour(fit(*np.indices(square.shape)))
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
            stars.append((x, y, width_x, width_y))

        # Show for how many stars the fitting failed
        if failed > 1: self._log.warning("Fitting stars failed for " + str(failed) + " out of " + str(self.regions[region].nshapes) + " objects")

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

        ## Define a function that calculates the area (in pixels) of a certain shape
        def area(shape):

            if shape.name == "circle": return math.pi * shape.coord_list[2] * shape.coord_list[2]
            if shape.name == "box": return shape.coord_list[2] * shape.coord_list[3]

        # For each shape in the specified region
        for shape in self.regions[region]:

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
        sky = np.zeros_like(self.frames.primary.data)

        # Create a row of pixels for the sky
        strip = np.arange(np.float(self.xsize))

        # Calculate the sky values on the grid
        for y in range(0, self.ysize):

            # Get the values of the polynomial
            sky[y,:] = polynomial(strip, y, parameters)

        # Add the sky as a new layer to this image
        self._addframe(sky, "sky")

    ## This function
    def findsky(self):

        # The length of the major axis of the ellipse
        major = 3.0 * self.orientation.majoraxis * 2.5

        # The width and heigth of the ellips
        width = major
        height = major * (1 - self.orientation.eps)

        # Create a string identifying this ellipse
        region_string = "image;ellipse(" + str(self.orientation.ypeak) + "," + str(self.orientation.xpeak) + "," + str(width) + "," + str(height) + "," + str(self.orientation.theta) + ")"

        # Create a region for the outer ellipse of the annulus
        region = pyregion.parse(region_string)

        # Add the annulus region
        self._addregion(region, "annulus")

        # Create the annulus mask
        annulusmask = np.logical_not(self.regions["annulus"]._region.get_mask(header=self.header, shape=(self.ysize,self.xsize)))

        # Get a combination of the currently selected masks
        currentmask = self.combinemasks()

        # Combine the currently selected mask, the galaxy mask and the annulus mask
        skymask = currentmask + self.masks.galaxy.data + annulusmask

        # Create a NumPy masked array
        maskedarray = np.ma.masked_array(self.frames.primary.data, skymask.astype(bool))

        # Make a mask of > 3 sigma regions
        newmaskedarray = sigma_clip(maskedarray, sig=3.0, iters=None, copy=True)

        # Add the mask
        self._addmask(newmaskedarray.mask, "sky")

        # Make a masked frame, the (sigma-clipped) sky
        skyframe = np.copy(self.frames.primary.data)

        # Set the sky frame to zero in the pixels masked by the new 'sky' mask
        skyframe[self.masks.sky.data] = 0

        # Add this frame to the set of frames
        self._addframe(skyframe, "sky")

        # Determine the mean, median and error of the sigma-clipped sky
        mean, median, error = sigma_clipped_stats(self.frames.primary.data, mask=skymask.astype(bool), sigma=3.0, iters=None)

        # Return the mean, median and standard deviation
        return mean, median, error

    ## This function fits a polynomial to the sky map
    def fitsky_lines(self):

        x = np.arange(self.xsize)

        bkg = np.zeros_like(self.frames.primary.data)

        for col in np.arange(self.ysize):

            weigths = np.logical_not(self.masks["sky"].data[col, x]).astype(int)

            number_of_ones = np.count_nonzero(weigths)

            if number_of_ones < 100: continue

            pfit = np.polyfit(x, self.frames.primary.data[col, x], 3, w=weigths)
            bkg[col, :] = np.polyval(pfit, x)

        # Add the new layer
        self._addframe(bkg, "fittedsky_lines")

    ## This function fits a polynomial to the sky map and saves the result as a frame called 'fittedsky'
    def fitsky(self, fwhm):

        # Determine the size of each box
        step = int(round(4 * fwhm))

        xvalues = []
        yvalues = []
        fluxes = []

        self._log.info("step = " + str(step))

        xx = np.arange(float(step)/2.0 + 1, float(self.xsize)-float(step)/2.0 - 1, float(step))
        yy = np.arange(float(step)/2.0 + 1, float(self.ysize)-float(step)/2.0 - 1, float(step))

        # Loop over all points in an evenly spaced grid (spacing = step)
        for x in xx:

            for y in yy:

                # Determine x range and y range of the box
                xmin = int(round(x - 0.5*step))
                xmax = int(round(x + 0.5*step))
                ymin = int(round(y - 0.5*step))
                ymax = int(round(y + 0.5*step))
                x_range = slice(xmin,xmax)
                y_range = slice(ymin,ymax)

                # Get the part of the sky mask that lies within this box
                box_mask = self.masks.sky.data[y_range, x_range]

                # Make a masked array from the part of the primary image that lies within this box
                maskedarray = np.ma.masked_array(self.frames.primary.data[y_range, x_range], box_mask)

                # Calculate the number of pixels in this box that are not masked
                covered_pixels = np.sum(np.logical_not(box_mask))

                # If this box does not include any pixels that are not masked, go to the next coordinate
                if covered_pixels == 0: continue

                # Calculate the mean flux in this box
                flux = maskedarray.mean()

                # Add the x coordinate, the y coordinate and the flux to the appropriate lists
                xvalues.append(x)
                yvalues.append(y)
                fluxes.append(flux)

        order = 3
        linear = True

        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)

        # Fit a polynomial
        parameters = fitpolynomial(xvalues, yvalues, fluxes, order, linear)

        # Image grid
        xx, yy = np.meshgrid(range(0, self.xsize), range(0, self.ysize))

        # Evaluate the polynomial on the image grid
        fittedsky = polynomial(xx.astype(float), yy, parameters)

        # Add the new layer
        self._addframe(fittedsky, "fittedsky")

    ## This function subtracts the currently active frame(s) from the primary image, in the pixels not covered by
    #  any of the currently active masks (the currently active masks 'protect' the primary image from this subtraction.
    def subtract(self):

        # Get the currently active mask
        totalmask = self.combinemasks()

        # Determine the negative of the total mask
        negativetotalmask = np.logical_not(totalmask)

        # For each active frame
        for frame in self.frames.getactive():

            # Inform the user
            self._log.info("Subtracting " + frame + " frame from the primary image frame")

            # Subtract the data in this frame from the primary image, in the pixels that the mask does not cover
            self.frames.primary.data -= self.frames[frame].data*negativetotalmask

    # ----------------------------------------------------------------- PSF DETERMINATION

    ## This function estimates the psf
    def estimatepsf_fitskirt(self, plot=False):

        # Get the stars
        stars = self.getstarpositions(plot)

        # Initialize two list to contain the full-width-half-maxima for the different stars
        xfwhm_list = []
        yfwhm_list = []

        # For each star in the list
        for star in stars:

            # Get the width in the x and y direction
            xfwhm_list.append(star[2])
            yfwhm_list.append(star[3])

        # Use sigma-clipped statistics to remove the outliers
        # and determine the average FWHM in the x and y direction for the remaining stars
        xfwhm, _, xfwhm_stdev = sigma_clipped_stats(xfwhm_list, sigma=3.0, iters=None)
        yfwhm, _, yfwhm_stdev = sigma_clipped_stats(yfwhm_list, sigma=3.0, iters=None)

        # Return the fwhm in the x and the y direction
        return xfwhm, yfwhm

# ----------------------------------------------------------------- USEFUL FUNCTIONS

## This function ...
def getpixelscale(header):

    # Initially, set the pixel scale to None
    pixelscale = None

    # Search for 'Pixel scale' keyword
    if 'PIXSCALE' in header:

        pixelscale = header['PIXSCALE']

    elif 'SECPIX' in header:

        pixelscale = header['SECPIX']

    # Search for 'Pixel Field of View' keyword
    elif 'PFOV' in header:

        pixelscale = header['PFOV']

    elif 'CD1_1' in header and 'CD1_2' in header:

        pixelscale = math.sqrt(header['CD1_1']**2 + header['CD1_2']**2 ) * 3600.0

    elif 'CD1_1' in header:

        pixelscale = abs(header['CD1_1']) * 3600.0

    elif 'CDELT1' in header:

        pixelscale = abs(header['CDELT1']) * 3600.0

    else:

        log = Log()
        log.warning("Could not determine the pixel scale from the image header")

    # Return the pixel scale (in arcseconds)
    return pixelscale

## This function
def getfilter(name, header):

    log = Log()

    # Initially, set the filter to None
    filter = None

    # Determine the filter from the information in the header
    if 'INSTRUME' in header and 'FILTER' in header:

        filterid = header['INSTRUME'].lower() + header['FILTER'].lower()

    # If no filter information could be found in the header, try to obtain it from the file name
    else:

        filterid = name.lower()

    # Parse the filterid
    if "fuv" in filterid: filter = Filter("GALEX.FUV")
    elif "pacs" in filterid:

        if '70' in filterid or 'blue' in filterid: filter = Filter("Pacs.blue")
        elif '160' in filterid or 'red' in filterid: filter = Filter("Pacs.red")
        else: log.warning("Could not determine the filter for this image")

    elif "mips" in filterid or "24" in filterid: filter = Filter("MIPS.24")
    elif "2mass" in filterid and "h" in filterid: filter = Filter("2MASS.H")
    elif "irac" in filterid:

        if '3.6' in filterid or 'i1' in filterid: filter = Filter("IRAC.I1")

    else: log.warning("Could not determine the filter for this image")

    # Return the filter
    return filter

## This function
def getunits(header):

    # Initially, set the units to None
    units = None

    if 'BUNIT' in header:

        units = header['BUNIT']

    # Return the units
    return units

## This function
def isskysub(header):

    # Initially, set the boolean to False
    subtracted = False

    if 'BACK_SUB' in header:

        subtracted = header['BACK_SUB']

    # Return the boolean value
    return subtracted

## This function
def getnumberofframes(header):

    # Initially, set the boolean to False
    nframes = 1

    if 'NAXIS' in header:

        # If there are 3 axes, get the size of the third
        if header['NAXIS'] == 3: nframes = header['NAXIS3']

    # Return the boolean value
    return nframes

## This function
def getframedescription(header, i):

    planeX = "PLANE" + str(i)

    # Get the description
    description = header[planeX]

    # Return the description
    return description

## This function
def getframename(description):

    # Convert spaces to underscores and ignore things between parentheses
    name = description.split("(")[0].rstrip(" ").replace(" ", "_")

    # If the frame name contains 'error', use the standard name "errors" for this frame
    if 'error' in name:

        name = "errors"

    # Return the frame name
    return name

# -----------------------------------------------------------------

## This class is a wrapper around the dict class, with the additional benefit of being able to access its values
#  with the 'dot' notation. It is a quite genious way of dealing with a set of layers (image frames, masks or regions
#  in this case), with high user-friendliness and easy programming interface.
#
class layers(dict):

    # Special trick to get elements from this dictionary
    def __getattr__(self, attr):

        return self.get(attr, None)

    # Set an item of the dictionary
    __setattr__= dict.__setitem__

    # Delete an item from the dictionary
    __delattr__= dict.__delitem__

    # This function returns a list of the names of the layers which are currently active
    def getactive(self):

        list = []

        # For each layer
        for name in self.keys():

            # If this layer is currently active, add it to the list
            if self[name].active: list.append(name)

        return list

    ## This function selects all the layers
    def selectall(self):

        # For each layer
        for name in self.keys():

            # Deselect this layer
            self[name].select()

    ## This function deselects all the layers
    def deselectall(self):

        # For each layer
        for name in self.keys():

            # Deselect this layer
            self[name].deselect()

    ## This function
    def list(self):

        log = Log()

        # For each layer
        for name in self.keys():

            # If this layer is active, print the name in green
            if self[name].active:

                log.success("        " + name)

            else:

                log.info("        " + name)

# -----------------------------------------------------------------

## Class ImageRegion
class ImageRegion(object):

    ## This function
    def __init__(self, region, log):

        # Set the internal pyregion object
        self._region = region

        # Logger
        self._log = log

        # Set as unactive initially
        self.active = False

    ## This function
    def select(self):

        self.active = True

    ## This function
    def deselect(self):

        self.active = False

    ## This function returns the number of shapes in this region
    @property
    def nshapes(self):

        return len(self._region)

# -----------------------------------------------------------------

## Class ImageMask
class ImageMask(object):

    ## This function
    def __init__(self, data, log):

        # Set the data array
        self.data = data

        # Logger
        self._log = log

        # Set as unactive initially
        self.active = False

    ## This function
    def select(self):

        self.active = True

    ## This function
    def deselect(self):

        self.active = False

# -----------------------------------------------------------------

## Class ImageFrame
class ImageFrame(object):

    ## The constructor
    def __init__(self, data, description, log):

        # Copy the data
        self.data = data

        # Logger
        self._log = log

        # Set as unactive initially
        self.active = False

        # Set the description
        self.description = description

    ## This function
    def select(self):

        self.active = True

    ## This function
    def deselect(self):

        self.active = False

    ## This function
    def histogram(self, path=None):

        NBINS = 1000
        plt.hist(self.data.flat, NBINS)

        # Display the result
        plt.show()

    ## This function
    @property
    def xsize(self):

        return self.data.shape[1]

    ## This function
    @property
    def ysize(self):

        return self.data.shape[0]

    ## This function
    def dtype(self):

        return self.data.dtype.name

    ## This function ...
    @property
    def mean(self):

        return np.mean(self.data)

    ## This function ...
    @property
    def median(self):

        return np.median(self.data)

    ## This function ...
    @property
    def min(self):

        return np.min(self.data)

    ## This function ...
    @property
    def max(self):

        return np.max(self.data)

    ## This function
    @property
    def stdev(self):

        # Set the delta degrees of freedom
        ddof = 1

        # Return the standard deviation of the data
        return np.std(self.data, ddof=ddof)

# -----------------------------------------------------------------
