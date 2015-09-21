#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import os
import os.path
import math
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

# Import astronomical modules
import aplpy
import pyregion
import astropy.io.fits as pyfits
from astropy import wcs
import astropy.units as u
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astropy import log

from photutils import CircularAperture
from photutils import aperture_photometry

# Import image modules
from tools import general, headers, cropping, interpolation, coordinates
import transformations
import fitting
import plotting
import analysis
import regions
import statistics
import catalogs
from layers import Layers
from frames import Frame
from masks import Mask
from regions import Region

# *****************************************************************

# Do not show warnings, to block Canopy's UserWarnings from spoiling the console log
import warnings
warnings.filterwarnings("ignore")

# *****************************************************************

class Image(object):

    """
    This class ...
    """

    def __init__(self, filename=None):

        """
        The constructor ...
        :param filename:
        :return:
        """

        if filename is None:
            
            # Initialize a set of layers to represent image frames, masks and regions
            self.frames = Layers()
            self.masks = Layers()
            self.regions = Layers()

            # Set default values for other attributes
            self.units = None
            self.fwhm = None
            self.header = None
            
            return

        # Check if the specified file exists, otherwise exit with an error
        if not os.path.isfile(filename): raise IOError("No such file: " + filename)

        # Set the name of the image
        self.name = os.path.splitext(os.path.basename(filename))[0]

        # Initialize a set of layers to represent image frames, masks and regions
        self.frames = Layers()
        self.masks = Layers()
        self.regions = Layers()
        
        # Read in the image
        self._load_image(filename)

        # Set default values for other attributes
        self.units = None
        self.fwhm = None

    # *****************************************************************

    def deselect_all(self):

        """
        This function ...
        :return:
        """

        self.frames.deselect_all()
        self.regions.deselect_all()
        self.masks.deselect_all()

    # *****************************************************************

    def info(self):

        """
        This function ...
        :return:
        """

        log.info("Name: " + self.name)
        log.info("Dimensions of data array: " + str(self.xsize) + " x " + str(self.ysize))
        log.info("Type of data: " + str(self.dtype))

    # *****************************************************************

    @property
    def dtype(self): return self.frames.primary.dtype

    # *****************************************************************

    @property
    def xsize(self): return self.frames.primary.xsize

    # *****************************************************************

    @property
    def ysize(self): return self.frames.primary.ysize

    # *****************************************************************

    @property
    def mean(self): return self.frames.primary.mean

    # *****************************************************************

    @property
    def median(self): return self.frames.primary.median

    # *****************************************************************

    @property
    def min(self): return self.frames.primary.min

    # *****************************************************************

    @property
    def max(self): return self.frames.primary.max

    # *****************************************************************

    @property
    def stddev(self): return self.frames.primary.stddev

    # *****************************************************************

    def import_datacube(self, path, name):

        """
        This function imports the datacube of a FITS file into this image
        :param path:
        :param name:
        :return:
        """

        # TODO: add coordinates !

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(path)

        # Get the primary data
        hdu = hdulist[0]
        
        header = hdu.header

        # Check whether multiple planes are present in the FITS image
        nframes = headers.get_number_of_frames(header)
        if nframes > 1:

            # For each frame
            for i in range(nframes):

                # Get the name of this frame, but the first frame always gets the name 'primary'
                description = headers.get_frame_description(header, i)
                frame_name = name + "_" + headers.get_frame_name(description) if i else name

                # Add this frame to the frames dictionary
                #self._add_frame(hdu.data[i], coordinates, frame_name, description)
                self._add_frame(hdu.data[i], None, frame_name, description)

        else:

            # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
            if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

            # Add the primary image frame
            self._add_frame(hdu.data, None, name)

        # Close the fits file
        hdulist.close()

    # *****************************************************************

    def export_datacube(self, filepath):

        """
        This function exports the currently selected frame(s) as a datacube into FITS file
        :param filepath:
        :return:
        """

        # Create an array to contain the data cube
        datacube = []

        # Get the coordinates of the primary frame
        coordinates = self.frames["primary"].coordinates
        
        # Construct a header for the image based on the coordinates of the primary frame
        header = coordinates.to_header() if coordinates is not None else None
        
        if header is None and self.header: header = self.header
        if header is None: header = pyfits.Header()

        plane_index = 1

        # Export all active frames to the specified file
        for frame_name in self.frames.get_selected():

            # Inform the user that this frame is being rebinned
            log.info("Exporting the " + frame_name + " frame to " + filepath)

            # Add this frame to the data cube, if its coordinates match those of the primary frame
            if coordinates == self.frames[frame_name].coordinates: datacube.append(self.frames[frame_name].data)
            
            # Add the name of the frame to the header
            header["PLANE"+str(plane_index)] = frame_name
            
            plane_index += 1

        # Create the HDU from the data array and the header
        hdu = pyfits.PrimaryHDU(np.array(datacube), header)

        # Write the HDU to a FITS file
        hdu.writeto(filepath, clobber=True)

        # Inform the user that the file has been created
        log.info("File " + filepath + " created")

    # *****************************************************************

    def import_region(self, path, name):

        """
        This function imports a new region from a DS9 region file
        :param path:
        :param name:
        :return:
        """

        # Create an pyregion object from the regions file
        region = pyregion.open(path)

        # Add the region to the set of regions
        self._add_region(region, name)

    # *****************************************************************

    def export_region(self, path):

        """
        This function exports the currently selected region(s) to a DS9 region file
        :param path:
        :return:
        """

        # Find the active region
        region = self.regions.get_selected(require_single=True)

        # Inform the user
        log.info("Creating " + path + " from the " + region + " region")

        # Write the region file
        self.regions[region].region.write(path)

    # *****************************************************************

    def get_state(self):

        """
        This function ...
        :return:
        """

        # Create an empty dictionary to contain the state of the current image
        state = dict()

        # Loop over all frames, regions and masks and record whether they are selected
        for frame_name in self.frames: state["frames/"+frame_name] = self.frames[frame_name].selected
        for region_name in self.regions: state["regions/"+region_name] = self.regions[region_name].selected
        for mask_name in self.masks: state["masks/"+mask_name] = self.masks[mask_name].selected

        # Return the state dictionary
        return state

    # *****************************************************************

    def set_state(self, state):

        """
        This function ...
        :param state:
        :return:
        """

        # Deselect all frames, regions and masks of this image
        self.deselect_all()

        # Loop over the entries in the state dictionary
        for identifier, selected in state.items():

            # Split the layer identifier into the layer type and the actual name of that layer
            layer_type, name = identifier.split("/")

            # Set the appropriate flag
            if layer_type == "frames": self.frames[name].selected = selected
            elif layer_type == "regions": self.regions[name].selected = selected
            elif layer_type == "masks": self.masks[name].selected = selected
            else: raise ValueError("Invalid state dictionary")

    # *****************************************************************

    def plot(self, path=None, color=True, grid=False, blacknan=False, publication=False):

        """
        This function shows a plot of the currently selected frame, combined with the active regions and masks
        :param path:
        :param color:
        :param grid:
        :param blacknan:
        :param publication:
        :return:
        """

        # Get the currently active frame
        frame = self.frames.get_selected()[0]

        # Create a total mask of the currently active masks
        total_mask = self.combine_masks(return_mask=True)

        # Mask the frame with nans
        maskedimage = np.ma.array(self.frames[frame].data, mask = total_mask)
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

        # If requested, use the 'publication' theme
        if publication: plot.set_theme('publication')

        # Add the regions
        for region in self.regions.get_selected():

            # Get the shape list
            shapes = self.regions[region].region.as_imagecoord(self.header)

            # Add these shapes to the plot
            plot.show_regions(shapes)

        if path is None:

            #plt.draw()
            #plt.close('all') # redundant
            #plt.show(block=False)
            plt.show()

        else: plot.save(path)

    # *****************************************************************

    def multiply(self, factor):

        """
        This function multiplies each selected frame with a given factor
        :param factor:
        :return:
        """

        # Fore each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Multiplying " + frame_name + " frame with a factor of " + str(factor))

            # Multiply this frame with the given factor
            self.frames[frame_name].data = self.frames[frame_name].data * factor

    # *****************************************************************

    def crop(self, x_min, x_max, y_min, y_max):

        """
        This function crops the currently selected frame(s) ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user that this frame is being cropped
            log.info("Cropping " + frame_name + " frame")

            # Crop this frame
            self.frames[frame_name].data = cropping.crop_check(self.frames[frame_name].data, x_min, x_max, y_min, y_max)

            # TODO: adjust coordinates!

    # *****************************************************************

    def rotate_and_center(self, left_x, left_y, right_x, right_y, flip=False):

        """
        This function rotates and centers the currently selected frame(s) ...
        :param left_x:
        :param left_y:
        :param right_x:
        :param right_y:
        :param flip:
        :return:
        """

        # Nans should be masked !

        shift_x = self.ysize/2 - (left_x + right_x)/2
        shift_y = self.xsize/2 - (left_y + right_y)/2

        log.info("Shifted frame to center by " + str(shift_x) + "," + str(shift_y))

        shiftframe = ndimage.interpolation.shift(self.frames.primary.data,(shift_x, shift_y))
        angle = math.degrees(math.atan(float(left_y - right_y)/float(left_x - right_x)))
        angle += 180.0 if flip else 0.0

        # Create the rotated frame
        rotframe = ndimage.interpolation.rotate(shiftframe, angle)

        # TODO: rotate the other layers, regions and masks!

        # Inform the user of the rotation angle
        log.info("Rotated frame over " + str(angle) + " degrees")

        # Add the new, rotated layer
        self._add_frame(rotframe, "primary_rotated")

    # *****************************************************************

    def rotate(self, angle):

        """
        This function rotates the currently selected frame(s) over a given angle
        :param angle: in degrees
        :return:
        """

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user that this frame is being rotated
            log.info("Rotating " + frame_name + " frame over " + str(angle) + " degrees")

            # Rotate this frame
            self.frames[frame_name].data = ndimage.interpolation.rotate(self.frames[frame_name].data, angle)

    # *****************************************************************

    def auto_rotate(self):

        """
        This function rotates the currently active frame(s) so that the galactic plane lies horizontal
        :return:
        """

        # Rotate about the position angle of the galaxy
        self.rotate(-self.orientation.theta)

    # *****************************************************************

    def auto_center(self):

        """
        This function makes a new layer where the center of the galaxy is in the center of the plane
        :return:
        """

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
        self._add_frame(centered, "primary_centered")

    # *****************************************************************

    def downsample(self, factor):

        """
        This function downsamples the currently selected frame(s) with a specified factor
        :param factor:
        :return:
        """

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Downsampling " + frame_name + " frame by a factor of " + str(factor))

            # Use the zoom function to resample
            self.frames[frame_name].data = ndimage.interpolation.zoom(self.frames[frame_name].data, zoom=1.0/factor)

    # *****************************************************************

    def set_units(self, units):

        """
        This function sets the units for this image
        :param units:
        :return:
        """

        # TODO: FIX THIS FUNCTION

        # Set the units for this image ...
        self.units = units

    # *****************************************************************

    def convert_to_units(self, units):

        """
        This function ...
        :param units:
        :return:
        """

        # TODO: FIX THIS FUNCTION

        # Calculate the conversion factor
        conversionfactor = self.units.convert(1.0, units)

        # Convert the data
        self.frames.primary.data *= conversionfactor

    # *****************************************************************

    def convert_to_magnitudes(self, m_0):

        """
        This function converts the currently selected frame(s) into magnitude scale, using the specified zero-point
        magnitude
        :param m_0:
        :return:
        """

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Converting " + frame_name + " frame to magnitude scale")

            # Convert to magnitude scale
            self.frames[frame_name].data = m_0 - 2.5 * np.log10(self.frames[frame_name].data)

    # *****************************************************************

    def convert_to_fluxes(self, F_0):

        """
        This function converts the currently selected frame(s) from magnitude to flux scale, using the specified
        zero-point flux
        :param F_0:
        :return:
        """

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Converting " + frame_name + " frame to flux scale")

            # Convert to flux scale
            self.frames[frame_name].data = F_0 * np.power(10.0, - self.frames[frame_name].data / 2.5)

    # *****************************************************************

    def set_orientation(self, orientation):

        """
        This function sets the orientation of the galaxy in this image
        :param orientation:
        :return:
        """

        self.orientation = orientation

    # *****************************************************************

    def set_fwhm(self, fwhm):

        """
        This function sets the FWHM of the PSF for this image
        :param fwhm:
        :return:
        """

        # Inform the user
        log.info("Setting the FWHM of the PSF for this image to " + str(fwhm) + " pixels")

        # Set the FWHM
        self.fwhm = fwhm

    # *****************************************************************

    def set_pixelscale(self, pixelscale):

        """
        This function sets the pixelscale for this image
        :param pixelscale: in arcseconds
        :return:
        """

        # Inform the user
        log.info("Setting the pixel scale of this image to " + str(pixelscale) + " arcseconds")

        # Set the pixel scale
        self.pixelscale = pixelscale

    # *****************************************************************

    def delete_frames(self):

        """
        This function removes the currently selected frame(s)
        :return:
        """

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Deleting the " + frame_name + " frame")

            # Remove this frame from the frames dictionary
            del self.frames[frame_name]

    # *****************************************************************

    def copy_frames(self):

        """
        This function ...
        :return:
        """

        # For each selected frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Copying the " + frame_name + " frame as another frame")

            # Copy the data and add it as a new frame
            data_copy = np.copy(self.frames[frame_name].data)
            coordinates = self.frames[frame_name].coordinates
            self._add_frame(data_copy, coordinates, frame_name+"_copy", description="Copy of the "+frame_name+" frame")

    # *****************************************************************

    def delete_regions(self):

        """
        This function removes the currently selected region(s)
        :return:
        """

        # For each active region
        for region_name in self.regions.get_selected(allow_none=False):

            # Inform the user
            log.info("Deleting the " + region_name + " region")

            # Remove this region from the regions dictionary
            del self.regions[region_name]

    # *****************************************************************

    def delete_masks(self):

        """
        This function removes the currently selected mask(s)
        :return:
        """

        # For each active mask
        for mask_name in self.masks.get_selected(allow_none=False):

            # Inform the user
            log.info("Deleting the " + mask_name + " mask")

            # Remove this mask from the masks dictionary
            del self.masks[mask_name]

    # *****************************************************************

    def uncertainty_from_regions(self):

        region = self.combine_regions(allow_none=False)

        frame_name = self.frames.get_selected(require_single=True)

        means = []
        stddevs = []

        for shape in region:

            x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

            box, x_min, x_max, y_min, y_max = cropping.crop(self.frames[frame_name].data, x_center, y_center, x_radius, y_radius)

            mean = np.mean(box)
            stddev = np.std(box)

            means.append(mean)
            stddevs.append(stddev)

        means = np.array(means)
        stddevs = np.array(stddevs)

        uncertainty = np.sqrt(np.std(means)**2 + np.median(stddevs))

        return uncertainty

    # *****************************************************************

    def convolve_fits(self, name):

        """
        This function convolves the currently selected frame(s) with a given kernel in the form of a FITS file
        :param name: the name of the kernel FITS file
        :return:
        """

        # The path to the kernel file
        path = os.path.join(os.getenv("HOME"), "Kernels", name)

        # Inform the user that the kernel was found
        log.info("Found kernel file at " + path)

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(path)

        # Get the kernel data and header
        kernel = hdulist[0].data
        header = hdulist[0].header

        # Inform the user
        log.info("Rebinning the kernel to the image pixel grid")

        # Get the pixel scale of the kernel
        pixelscale_kernel = header["CD1_1"]*3600

        print self.pixelscale

        self.pixelscale = headers.get_pixelscale(header)

        print self.pixelscale

        # Calculate the zooming factor
        factor = self.pixelscale / pixelscale_kernel

        # Rebin the kernel to the same grid of the image
        kernel = ndimage.interpolation.zoom(kernel, zoom=1.0/factor)

        print kernel

        #print kernel[kernel.shape[0]/2.0,kernel.shape[1]/2.0]

        # For all active frames, do the convolution
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user that this frame is being convolved
            log.info("Convolving " + frame_name + " frame with the kernel " + os.path.splitext(name)[0])

            # Do the convolution on this frame
            self.frames[frame_name].data = convolve_fft(self.frames[frame_name].data, kernel, normalize_kernel=True)

        # Close the FITS file
        hdulist.close()

    # *****************************************************************

    def convolve_model(self, fwhm, model="Gaussian"):

        """
        This function convolves the currently selected frame(s) with an analytical kernel model
        :param fwhm:
        :param model:
        :return:
        """

        # From the FWHM, calculate the standard deviation
        sigma = fwhm / 2.355

        # Construct the Gaussian kernel
        kernel = Gaussian2DKernel(sigma)

        # For all active frames, do the convolution
        for frame_name in self.frames.get_selected():

            # Inform the user that this frame is being convolved
            log.info("Convolving " + frame_name + " frame with a Gaussian kernel with a FWHM of " + str(fwhm))

            # Do the convolution on this frame
            self.frames[frame_name].data = convolve(self.frames[frame_name].data, kernel)

    # *****************************************************************

    def rebin(self, reference):

        """
        This function rebins the currently selected frame(s) based on a certain reference FITS file
        :param reference:
        :return:
        """

        # Get the header of the reference image
        refheader = pyfits.getheader(reference)

        # Obtain the coordinate system from the reference header
        coordinates = wcs.WCS(refheader)

        # For all active frames, do the rebinning
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user that this frame is being rebinned
            log.info("Rebinning " + frame_name + " frame to the grid of " + reference)

            # Do the rebinning based on the header of the reference image
            self.frames[frame_name].data = transformations.align_and_rebin(self.frames[frame_name].data, self.header, refheader)

            # Set the new coordinate system for this frame
            self.frames[frame_name].coordinates = coordinates

    # *****************************************************************

    def interpolate(self):

        """
        This function interpolates the currently selected frame(s) within the currently active mask(s)
        :return:
        """

        # Inform the user
        log.info("Interpolating the image in the areas covered by the currently selected masks")

        # Combine all currently selected masks
        total_mask = self.combine_masks(return_mask=True)

        # Loop over all currently selected frames
        for frame_name in self.frames.get_selected():

            # Perform the interpolation on this frame
            self.frames[frame_name].data = interpolation.in_paint(self.frames[frame_name].data, total_mask)

    # *****************************************************************

    def interpolate_in_regions(self, sigma_clipping=True):

        """
        This function ...
        :param sigma_clipping:
        :return:
        """

        # Combine the active masks
        total_mask = self.combine_masks(allow_none=False, return_mask=True)

        # Combine all the active regions
        region = self.combine_regions(allow_none=False)

        # Make a mask from the region and invert it (to mask outside of the region)
        #region_mask = np.logical_not(regions.create_mask(region, self.header, self.xsize, self.ysize))

        # Loop over all active frames
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Interpolating the " + frame_name + " frame within the areas covered by the currently "
                     "selected masks, if enclosed by any of the currently selected regions")

            # Loop over all shapes
            for shape in region:

                # Get the extents of the box that encloses this shape
                x_min, x_max, y_min, y_max = regions.get_enclosing_box(shape)

                # Cut out the box
                box, x_min, x_max, y_min, y_max = cropping.crop_direct(self.frames[frame_name].data, x_min, x_max, y_min, y_max)

                # Cut out the mask
                box_mask = cropping.crop_check(total_mask, x_min, x_max, y_min, y_max)

                # Cut out the region mask
                #box_region_mask = cropping.crop_check(region_mask, x_min, x_max, y_min, y_max)

                # Make a combined mask for the interpolation
                #interpolation_mask = box_mask | box_region_mask

                interpolation_mask = box_mask

                # If sigma clipping is enabled, mask additional pixels if they are outliers
                if sigma_clipping: interpolation_mask = statistics.sigma_clip_mask(box, mask=interpolation_mask)

                # Calculate the interpolated background
                interpolated_box = interpolation.in_paint(box, interpolation_mask)

                # If the interpolated box contains nans, do not fill in the corresponding pixels of the data with these nans,
                # therefore set the pixels that are nan to False in the box_mask (take the difference between the box_mask
                # and the np.isnan(interpolated_box) mask). Then, set the nans to zero in the interpolated_box because
                # False * nan would otherwise still equal to nan.
                box_mask = masks.subtract(box_mask, np.isnan(interpolated_box))
                interpolated_box[np.isnan(interpolated_box)] = 0.0

                # Insert the interpolated values, for the pixels that are masked by the total currently selected mask
                self.frames[frame_name].data[y_min:y_max,x_min:x_max] = interpolated_box*box_mask + self.frames[frame_name].data[y_min:y_max,x_min:x_max]*np.logical_not(box_mask)

    # *****************************************************************

    def mask_nans(self):

        """
        This function masks the 'NaN' values in the currently selected frame(s)
        :return:
        """

        # Get the name of the currently selected frame
        frame_name = self.frames.get_selected(require_single=True)

        # Get a map of all the NaNs in the primary image
        mask = np.isnan(self.frames[frame_name].data)

        # Make a nans mask layer
        self._add_mask(mask, "nans")

    # *****************************************************************

    def expand_masks(self, name, iterations=100):

        """
        This function ...
        :param name:
        :param iterations:
        :return:
        """

        # Define the structure for the expansion
        structure = ndimage.generate_binary_structure(2, 2)

        # Get a combination of the active masks
        old_mask = self.combine_masks(return_mask=True)

        # Make the new mask, made from 100 iterations with the structure array
        newmask = ndimage.binary_dilation(old_mask, structure, iterations)

        # Add the new, expanded mask
        self._add_mask(newmask, name)

    # *****************************************************************

    def combine_regions(self, name=None, allow_none=True):

        """
        This function ...
        :param name:
        :param allow_none:
        :return:
        """

        # Initialize an empty list of shapes
        total_region = pyregion.ShapeList([])

        # TODO: what to do if one region is in image coordinates and other in physical coordinates?
        # Temporary fix: all in image coordinates

        # Loop over all active regions, adding them together
        for region_name in self.regions.get_selected(allow_none=allow_none):

            # Add all the shapes of this region to the combined region
            for shape in self.regions[region_name].region.as_imagecoord(self.header):

                total_region.append(shape)

        # If no name is provided, return the new region
        if name is None: return total_region

        # Else, add the region to the list of regions, with the appropriate name
        else: self._add_region(total_region, name)

    # *****************************************************************

    def combine_masks(self, name=None, allow_none=True, return_mask=False):

        """
        This function ...
        :param name:
        :param allow_none:
        :return:
        """

        # Initialize an boolean array for the total mask
        total_mask = np.zeros_like(self.frames.primary.data, dtype=bool)

        # For each active mask
        for mask_name in self.masks.get_selected(allow_none=allow_none):

            print self.masks[mask_name].data.shape

            # Add this mask to the total
            total_mask += self.masks[mask_name].data

        # Set the name of the total mask
        name = name if name is not None else "total"

        # Return the mask or add it to this image
        if return_mask: return total_mask
        else: self._add_mask(total_mask, name)

    # *****************************************************************

    def invert_mask(self, name):

        """
        This function makes a new mask which is the inverse (logical NOT) of the total currently selected mask
        :param name:
        :return:
        """

        # Get the total selected mask
        currentmask = self.combine_masks(return_mask=True)

        # Calculate the inverse of the this total mask
        newmask = np.logical_not(currentmask)

        # Add the new, inverted mask
        self._add_mask(newmask, name)

    # *****************************************************************

    def apply_masks(self, fill=0.0):

        """
        This function applies the currently selected mask(s) to the currently selected frame(s)
        :param fill:
        :return:
        """

        # Get the total selected mask
        total_mask = self.combine_masks(allow_none=False, return_mask=True)

        # For each active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # TODO: check if the dimensions of frame and mask match!

            # Inform the user
            log.info("Applying the total selected mask to the " + frame_name + " frame")

            # Set the corresponding image pixels to zero for this mask
            self.frames[frame_name].data[total_mask] = fill

    # *****************************************************************

    def create_mask(self, return_mask=False):

        """
        This function creates a mask from the currently selected region(s)
        :return:
        """

        # TODO: use combine_regions and regions.create_mask() !!

        # Initialize an boolean array for the total mask
        total_mask = np.zeros_like(self.frames.primary.data, dtype=bool)

        name = ""

        # For each active region
        for region_name in self.regions.get_selected():

            # Create the mask
            total_mask += regions.create_mask(self.regions[region_name].region, self.header, self.xsize, self.ysize)
            name += region_name + "_"

        # Remove the trailing underscore
        name = name.rstrip("_")

        # Return the mask or add it to this image
        if return_mask: return total_mask
        else: self._add_mask(total_mask, name)

    # *****************************************************************

    def find_galaxy(self, galaxy_name, name=None, plot=False, return_region=False):

        """
        This function ...
        :param galaxy_name:
        :param name:
        :param plot:
        :param return_region:
        :return:
        """

        # Get the name of the active frame
        frame_name = self.frames.get_selected(require_single=True)

        # Get the region with the galaxy position fetched from an online catalog
        region = self.fetch_galaxy(galaxy_name, return_region=True).as_imagecoord(self.header)

        # Find the galaxy
        galaxy_parameters = analysis.find_galaxy_orientation(self.frames[frame_name].data, region, plot=plot)

        # Create a region with one ellipse corresponding to the extent of the galaxy in this frame
        galaxy_region = regions.one_ellipse(galaxy_parameters)

        # Set the name of the new region
        name = name if name is not None else "galaxy"

        # Return the region or add it to this image
        if return_region: return galaxy_region
        else: self._add_region(galaxy_region, name)

    # *****************************************************************

    def fetch_galaxy(self, galaxy_name, name=None, radius=100, return_region=False):

        """
        This function ...
        :param galaxy_name:
        :param name:
        :param radius:
        :param return_region:
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy position from an online catalog...")

        # Search for the position of the specified galaxy in the image
        region = catalogs.fetch_object_by_name(galaxy_name, radius)

        # Set the name of the new region
        name = name if name is not None else "galaxy"

        # Return the region or add it to this image
        if return_region: return region
        else: self._add_region(region, name)

    # *****************************************************************

    def perform_mge(self):

        pass

    # *****************************************************************

    def get_galactic_extinction(self, galaxy_name):

        """
        This function ...
        """

        return catalogs.fetch_galactic_extinction(galaxy_name, self.filter)

    # *****************************************************************

    def fetch_stars(self, radius, catalog=["UCAC4"], galaxy_name=None, return_region=False, column_filters=None):

        """
        This function fetches the positions of astrophysical objects in the image
        :param name:
        :param color:
        :param return_region:
        :param radius:
        :return:
        """

        # Possible catalogs: "UCAC4", "NOMAD", "PPMXL" ? or combinations

        # Inform the user
        log.info("Fetching star positions from an online catalog...")

        # Get the range of RA and DEC of this image
        ra_center, dec_center, size_ra_deg, size_dec_deg = self._get_coordinate_range()

        # Search for stars
        radius_in_arcsec = radius * self.pixelscale
        box = (ra_center, dec_center, size_ra_deg, size_dec_deg)
        if column_filters is None: region = catalogs.fetch_objects_in_box(box, catalog, ["stars"], radius_in_arcsec)
        else: region = catalogs.fetch_objects_in_box(box, catalog, ["optical", "stars"], radius_in_arcsec, column_filters={"Vmag":"<20"})

        if galaxy_name is not None:

            # Search for the position of the specified galaxy in the image
            galaxy_region = catalogs.fetch_object_by_name(galaxy_name, radius)
            original_len = len(region)
            region = regions.subtract(region, galaxy_region, 3.0, self.header)

            if len(region) == original_len - 1: log.info("Removed star position that matches the galaxy center")
            elif len(region) < original_len - 1: log.warning("Removed multiple star positions")

        # Inform the user
        log.info("Found " + str(len(region)) + " entries")

        # Return the region or add it to this image
        if return_region: return region
        else: self._add_region(region, "stars")

    # *****************************************************************

    def create_annulus(self):

        """
        This function ...
        :return:
        """

        pass

        # Circle --> annulus
        # Ellipse --> ellipse annulus
        # Box --> box annulus

    # *****************************************************************

    def expand_regions(self, factor, combine=False):

        """
        This function expands the currently selected region(s)
        :param factor:
        :param combine:
        :return:
        """

        if combine:

            # Create a combined region
            region = self.combine_regions(allow_none=False)

            # Expand this region
            expanded_region = regions.expand(region, factor)

            # Add this region to the list of regions
            self._add_region(expanded_region, "expanded")

        else:

            # Loop over all active regions
            for region_name in self.regions.get_selected(allow_none=False):

                # Inform the user
                log.info("Expanding the " + region_name + " region by a factor of " + str(factor))

                # Create expanded region
                expanded_region = regions.expand(self.regions[region_name].region, factor)

                # Add the expanded region to the list of regions
                self._add_region(expanded_region, region_name + "_expanded")

    # *****************************************************************

    def rename_region(self, name):

        """
        This function renames a region
        :param name:
        :return:
        """

        # Get the name of the currently selected region
        region_name = self.regions.get_selected(require_single=True)

        # Remove the region of the dictionary of regions and re-add it under a different key
        self.regions[name] = self.regions.pop(region_name)

    # *****************************************************************

    def rename_frame(self, name):

        """
        This function renames a frame
        :param name:
        :return:
        """

        # Get the name of the currently selected frame
        frame_name = self.frames.get_selected(require_single=True)

        # Remove the frame from the dictionary of frames and re-add it under a different key
        self.frames[frame_name] = self.frames.pop(frame_name)

    # *****************************************************************

    def model_stars(self, model_name='Gaussian', background_inner_sigmas=5.0, background_outer_sigmas=10.0, fit_sigmas=5.0,
                    upsample_factor=1.0, interpolate_background=True, sigma_clip_background=True, plot=False):

        """
        This function ...
        :param model:
        :param plot:
        :param background_inner_sigmas:
        :param background_outer_sigmas:
        :param fit_sigmas:
        :param resample_factor:
        :param interpolate_background:
        :param sigma_clip_background:
        :return:
        """

        # TODO: should not only work with regions that contain nothing but ellipses

        # Create a region for stars that were succesfully modeled and a region for objects that could not be fitted to a star model
        modeled = pyregion.ShapeList([])
        unmodeled = pyregion.ShapeList([])

        # Get the name of the active frame
        frame_name = self.frames.get_selected(require_single=True)

        # Inform the user
        log.info("Modeling stars in the " + frame_name + " frame enclosed by any of the currently selected regions")

        # Create a new frame to contain the modeled stars
        stars_frame = np.zeros_like(self.frames[frame_name].data)

        # Combine all the active regions
        total_region = self.combine_regions(allow_none=False)

        # Create the background mask
        annuli_mask = masks.annuli_around(total_region, background_inner_sigmas, background_outer_sigmas, self.header, self.xsize, self.ysize)

        # Create a mask that covers pixels too far away from the center of the star (for fitting the model)
        fit_mask = masks.masked_outside(total_region, self.header, self.xsize, self.ysize, expand_factor=fit_sigmas)

        # For each shape (star)
        for shape in total_region:

            # Try to make an analytical model for the star enclosed by this shape
            success, shape, model, extents = analysis.make_star_model(shape, self.frames[frame_name].data,
                                                                      annuli_mask, fit_mask, background_outer_sigmas,
                                                                      fit_sigmas, model_name,
                                                                      upsample_factor=upsample_factor,
                                                                      interpolate_background=interpolate_background,
                                                                      sigma_clip_background=sigma_clip_background,
                                                                      plot=plot)

            if success:

                # Add the 1-sigma contour of the analytical model to the modeled stars region
                modeled.append(shape)

                # Add the model to the stars frame
                stars_frame[extents[2]:extents[3], extents[0]:extents[1]] += model

            else: unmodeled.append(shape)

        # Add the modelled stars frame to the list of frames
        self._add_frame(stars_frame, self.frames[frame_name].coordinates, "stars")

        # Add the successfully modeled stars and the unmodeled stars to the corresponding regions
        self._add_region(modeled, "modeled_stars")
        self._add_region(unmodeled, "unmodeled_stars")

    # *****************************************************************

    def find_stars(self, galaxy_name=None, plot=False, plot_custom=[False, False, False, False], catalog=["UCAC4"],
                   detection_method="peaks", initial_radius=10.0, in_region=False, split_sigma=None, model="Gaussian",
                   split_brightness=False, brightness_sigma=6.0):

        """
        This function searches for stars in the currently selected frame, by fetching star positions from an
        online catalog and looking for unique sources around these positions
        :param plot: flag that indicates whether the found stars should be plotted
        :param plot_custom: additional flags to indicate whether specific cases should be plotted
        """

        # Get the name of the active frame
        frame_name = self.frames.get_selected(require_single=True)

        # Get the region of objects fetched from the NOMAD stellar catalog and transform it into image coordinates
        if in_region: region = self.combine_regions(allow_none=False)
        else: region = self.fetch_stars(initial_radius, catalog=catalog, galaxy_name=galaxy_name, return_region=True).as_imagecoord(self.header)

        # Look for sources
        sources = analysis.find_sources_in_region(self.frames[frame_name].data, region, model, detection_method, plot=plot, plot_custom=plot_custom)
        log.info("Number of sources = " + str(len(sources)))

        # Remove duplicates
        unique_sources = analysis.remove_duplicate_sources(sources)
        log.info("Number of unique sources = " + str(len(unique_sources)))

        # Split the unique sources into stars and unidentified objects
        if split_sigma is None:

            stars = unique_sources
            ufos = []

        else:

            stars, ufos = statistics.sigma_clip_split(unique_sources, lambda source: general.average_stddev(source), sigma=split_sigma)
            log.info("Number of stars = " + str(len(stars)))
            log.info("Number of unidentified objects = " + str(len(ufos)))

        # Only create a region if any stars were found
        if len(stars) > 0:

            if split_brightness:

                # Find saturated stars
                stars, brightest = statistics.split_percentage(stars, lambda source: source.amplitude, percentage=0.05)

                #stars, brightest = statistics.sigma_clip_split(stars, lambda source: source.amplitude, sigma=brightness_sigma, only_high=True)

                if len(brightest) > 0:

                    # Convert the list of saturated stars to a region and add it to the list of regions
                    brightest_region = regions.ellipses_from_coordinates(brightest)
                    self._add_region(brightest_region, "brightest")

            # Convert the list of stars to a region and add it to the list of regions
            stars_region = regions.ellipses_from_coordinates(stars)
            self._add_region(stars_region, "stars")

        else: log.warning("No stars could be detected")

        # Only create a region if any unidentified objects were found
        if len(ufos) > 0:

            # Convert the list of ufos to a region and add it to the list of regions
            ufos_region = regions.ellipses_from_coordinates(ufos)
            self._add_region(ufos_region, "ufos")

    # *****************************************************************

    def split_region(self, criterium, method, percentage=None, sigma=None):

        """
        This function ...
        """

        #criterium = "flux", "center_brightness", "radius" ...

        # method = "sigma_clip" or "percentage"

        # Get the total selected region
        region = self.combine_regions(allow_none=False)

        frame_name = self.frames.get_selected(require_single=True)

        def center_brightness(shape):

            x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

            print self.frames[frame_name].data.shape

            x = int(round(x_center))
            y = int(round(y_center))

            if x < self.frames[frame_name].data.shape[1] and x >= 0 and y < self.frames[frame_name].data.shape[0] and y >= 0:

                return self.frames[frame_name].data[y, x]

            else: return 0.0

        def flux(shape):

            x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

            aperture = CircularAperture((x_center, y_center), r=x_radius)

            phot_table = aperture_photometry(self.frames[frame_name].data, aperture, mask=np.isnan(self.frames[frame_name].data))

            return phot_table[0]["aperture_sum"]

        if criterium == "flux": key = flux
        elif criterium == "center_brightness": key = center_brightness
        elif criterium == "radius": key = lambda shape: shape.coord_list[2]
        else: raise ValueError("Not a valid criterium")

        if method == "sigma_clip":

            dim_region, bright_region = statistics.sigma_clip_split(region, key, sigma=sigma, only_high=True, nans="high")

        elif method == "percentage":

            dim_region, bright_region = statistics.split_percentage(region, key, percentage=percentage, nans="high")

        else: raise ValueError("Not a valid splitting method")

        self._add_region(dim_region, "dim")

        if len(bright_region) > 0:

            print len(bright_region)

            self._add_region(bright_region, "bright")

    # *****************************************************************

    def create_segmentation_mask(self, kernel_fwhm, kernel_size, threshold_sigmas=2.0, expand=True, plot=False):

        """
        This function ...
        :return:
        """

        # Get the total selected region
        region = self.combine_regions()

        if len(region) == 0: pass # Use whole frame

        # Get the name of the currently selected frame
        frame_name = self.frames.get_selected(require_single=True)

        mask = np.zeros_like(self.frames[frame_name].data, dtype=bool)

        # Loop over all shapes
        for shape in region:

            if not regions.in_box(shape, self.frames[frame_name].data.shape): break

            box_mask, x_min, x_max, y_min, y_max = analysis.find_center_segment_in_shape(self.frames[frame_name].data,
                                                                                         shape, kernel_fwhm, kernel_size,
                                                                                         threshold_sigmas, expand=expand,
                                                                                         expansion_level=1,
                                                                                         max_expansion_level=5,
                                                                                         plot=plot)

            # Adapt the mask
            mask[y_min:y_max, x_min:x_max] = masks.union(mask[y_min:y_max, x_min:x_max], box_mask)

        # Add the new mask
        self._add_mask(mask, "segments")

    # *****************************************************************

    def find_sky(self):

        """
        This function ...
        :return:
        """

        # Get the name of the currently selected frame
        frame_name = self.frames.get_selected(require_single=True)

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
        self._add_region(region, "annulus")

        # Create the annulus mask
        annulusmask = np.logical_not(self.regions["annulus"].region.get_mask(header=self.header, shape=(self.ysize,self.xsize)))

        # Get a combination of the currently selected masks
        current_mask = self.combine_masks(return_mask=True)

        # Combine the currently selected mask, the galaxy mask and the annulus mask
        sky_mask = (current_mask + self.masks.galaxy.data + annulusmask).astype(bool)

        # Make a mask of > 3 sigma regions
        new_mask = statistics.sigma_clip_mask(self.frames[frame_name].data, sigma=3.0, mask=sky_mask)

        # Add the mask
        self._add_mask(new_mask, "sky")

        # Make a masked frame, the (sigma-clipped) sky
        skyframe = np.copy(self.frames.primary.data)

        # Set the sky frame to zero in the pixels masked by the new 'sky' mask
        skyframe[self.masks.sky.data] = 0

        # Add this frame to the set of frames
        self._add_frame(skyframe, None, "sky")

    # *****************************************************************

    def estimate_background(self, downsample_factor, plot=False):

        # Get the name of the currently selected frame
        frame_name = self.frames.get_selected(require_single=True)

        # Combine all currently selected masks
        total_mask = self.combine_masks(return_mask=True)

        # Estimate the background
        background = interpolation.low_res_interpolation(self.frames[frame_name].data, downsample_factor, mask=total_mask)

        # Plot the difference between the data and the model, if requested
        if plot: plotting.plot_difference(self.frames[frame_name].data, background)

        # Add the background frame
        self._add_frame(background, self.frames[frame_name].coordinates, "background")

    # *****************************************************************

    def fit_polynomial(self, plot=False, degree=3, sigma_clipping=True):

        """
        This function fits a polynomial function to each of the currently active frames
        :param plot:
        :param upsample_factor:
        :return:
        """

        # Get the currently active mask
        total_mask = self.combine_masks(return_mask=True)

        # For each currently active frame
        for frame_name in self.frames.get_selected(allow_none=False):

            # Inform the user
            log.info("Fitting a polynomial function to the " + frame_name + " frame")

            if sigma_clipping: new_mask = statistics.sigma_clip_mask(self.frames[frame_name].data, mask=total_mask)
            else: new_mask = total_mask

            # Fit the model
            polynomial = fitting.fit_polynomial(self.frames[frame_name].data, mask=new_mask, degree=degree)

            # Plot the difference between the data and the model, if requested
            if plot: plotting.plot_difference_model(self.frames[frame_name].data, polynomial)

            # Evaluate the model
            evaluated = fitting.evaluate_model(polynomial, 0, self.frames[frame_name].data.shape[1], 0, self.frames[frame_name].data.shape[0])

            # Add the evaluated model as a new frame
            description = "A polynomial fit to the " + frame_name + " primary frame"
            self._add_frame(evaluated, self.frames[frame_name].coordinates, frame_name+"_polynomial", description)

    # *****************************************************************

    def subtract(self):

        """
        This function subtracts the currently active frame(s) from the primary image, in the pixels not covered by
        any of the currently active masks (the currently active masks 'protect' the primary image from this subtraction
        :return:
        """

        # Get the currently active mask
        total_mask = self.combine_masks(return_mask=True)

        # Determine the negative of the total mask
        negativetotalmask = np.logical_not(total_mask)

        # For each active frame
        for frame_name in self.frames.get_selected():

            # Inform the user
            log.info("Subtracting " + frame_name + " frame from the primary image frame")

            # Subtract the data in this frame from the primary image, in the pixels that the mask does not cover
            self.frames.primary.data -= self.frames[frame_name].data*negativetotalmask

    # *****************************************************************

    def mean_radius(self):

        """
        This function calculates the mean radius of the shapes (ellipses, circles) in the currently selected region(s)
        :return:
        """

        # TODO: fix for regions that do not only contain ellipses

        # Initialize an empty list to contain the different sigma values
        sigmas = []

        # Loop over all currently selected regions
        for region_name in self.regions.get_selected(allow_none=False):

            # Loop over all shapes in the region
            for shape in self.regions[region_name].region:

                sigma_x = shape.coord_list[2]
                sigma_y = shape.coord_list[3]

                # Add the sigma, averaged over the x and y directions, to the list of sigmas
                sigmas.append(0.5*(sigma_x + sigma_y))

        # Return the mean value of the sigmas
        return np.mean(sigmas)

    # *****************************************************************

    def _load_image(self, filename):

        """
        This function ...
        :param filename:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file: " + filename)

        # Open the HDU list for the FITS file
        hdulist = pyfits.open(filename)

        # Get the primary HDU
        hdu = hdulist[0]

        # Get the image header
        header = hdu.header

        # Obtain the coordinate system
        coordinates = wcs.WCS(header)

        # Load the frames
        self.pixelscale = headers.get_pixelscale(header)

        # Obtain the filter for this image
        self.filter = headers.get_filter(self.name, header)

        # Obtain the units of this image
        self.units = headers.get_units(header)

        # Check whether multiple planes are present in the FITS image
        nframes = headers.get_number_of_frames(header)
        if nframes > 1:

            # For each frame
            for i in range(nframes):

                # Get the name of this frame, but the first frame always gets the name 'primary'
                description = headers.get_frame_description(header, i) if i else "the primary signal map"    
                name = headers.get_frame_name(description) if i else "primary"

                # Add this frame to the frames dictionary
                self._add_frame(hdu.data[i], coordinates, name, description)

        else:

            # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
            if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

            # Add the primary image frame
            self._add_frame(hdu.data, coordinates, "primary", "the primary signal map")

        # Set the basic header for this image
        self.header = header.copy(strip=True)
        self.header["NAXIS"] = 2
        self.header["NAXIS1"] = self.xsize
        self.header["NAXIS2"] = self.ysize

        # Select the primary image frame
        self.frames.primary.select()

        # Close the FITS file
        hdulist.close()

    # *****************************************************************

    def _add_frame(self, data, coordinates, name, description=None):

        """
        This function ...
        :param data:
        :param coordinates:
        :param name:
        :param description:
        :return:
        """

        # Inform the user
        log.info("Adding '" + name + "' to the set of image frames")

        # Add the layer to the layers dictionary
        self.frames[name] = Frame(data, coordinates, description)

    # *****************************************************************

    def _add_region(self, region, name):

        """
        This function ...
        :param region:
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adding '" + name + "' to the set of regions")

        # Add the region to the regions dictionary
        self.regions[name] = Region(region)

    # *****************************************************************

    def _add_mask(self, data, name):

        """
        This function ...
        :param data:
        :param name:
        :return:
        """

        # Inform the user
        log.info("Adding '" + name + "' to the set of masks")

        # Add the mask to the masks dictionary
        self.masks[name] = Mask(data)

    # *****************************************************************

    def _get_coordinate_range(self):

        """
        This function ...
        :return:
        """

        w = wcs.WCS(self.header)

        # Some pixel coordinates of interest.
        pixels = np.array([[0,0],[self.xsize-1,self.ysize-1]], dtype=float)

        world = w.wcs_pix2world(pixels, 0)  # Convert pixel coordinates to world coordinates (RA and DEC in degrees)
        coordinate1 = world[0]
        coordinate2 = world[1]
        ra_range = [coordinate2[0], coordinate1[0]]
        dec_range = [coordinate1[1], coordinate2[1]]

        # Determine the center in RA and DEC (in degrees)
        ra_center = 0.5*(ra_range[0] + ra_range[1])
        dec_center = 0.5*(dec_range[0] + dec_range[1])

        # Determine the width in RA and DEC (both in degrees)
        dec_width = dec_range[1] - dec_range[0]
        ra_width = ra_range[1] - ra_range[0]   # WRONG!

        # Calculate the start and end RA coordinates (in degrees)
        ra_begin = ra_center - 0.5*ra_width
        ra_end = ra_center + 0.5*ra_width

        # Calculate the start and end DEC coordinates (in degrees)
        dec_begin = dec_center - 0.5*dec_width
        dec_end = dec_center + 0.5*dec_width

        # Calculate the
        ra_distance = coordinates.ra_distance(dec_center, ra_begin, ra_end)
        dec_distance = dec_end - dec_begin

        # Calculate the pixel scale of this image in degrees
        pixelscale = self.pixelscale * u.arcsec
        pixelscale_deg = pixelscale.to("deg").value

        # Get the center pixel
        ref_pix = w.wcs.crpix
        ref_world = w.wcs.crval

        # Get the number of pixels
        size_dec_deg = self.ysize * pixelscale_deg
        size_ra_deg = self.xsize * pixelscale_deg

        # Check whether the two different ways of calculating the RA width result in approximately the same value
        assert np.isclose(ra_distance, size_ra_deg, rtol=0.02), "The coordinate system and pixel scale do not match"
        assert np.isclose(dec_distance, size_dec_deg, rtol=0.02), "The coordinate system and pixel scale do not match"

        # Return ...
        return (ra_center, dec_center, size_ra_deg, size_dec_deg)

    # *****************************************************************