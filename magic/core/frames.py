#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

# Import Astromagic modules
from ..tools import coordinates
from ..tools import cropping
from ..tools import transformations
from .vector import Position, Extent

# Import astronomical modules
import aplpy
import astropy.io.fits as pyfits
import astropy.units as u
import astropy.coordinates as coord
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel

# *****************************************************************

class Frame(np.ndarray):

    """
    This class ...
    """

    # *****************************************************************

    def __new__(cls, data, wcs=None, pixelscale=None, description=None, selected=False, unit=None):

        """
        This function ...
        :param cls:
        :param input_array:
        :param info:
        :return:
        """

        obj = np.asarray(data).view(cls)
        obj.wcs = wcs
        obj.pixelscale = pixelscale
        obj.description = description
        obj.selected = selected
        obj.unit = unit

        return obj

    # *****************************************************************

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        if obj is None: return
        self.wcs = getattr(obj, 'wcs', None)
        self.pixelscale = getattr(obj, 'wcs', None)
        self.description = getattr(obj, 'description', None)
        self.selected = getattr(obj, 'selected', False)
        self.unit = getattr(obj, 'unit', None)

    # *****************************************************************

    def select(self):

        """
        This function ...
        :return:
        """

        self.selected = True

    # *****************************************************************

    def deselect(self):

        """
        This function ...
        :return:
        """

        self.selected = False

    # *****************************************************************

    @property
    def xsize(self): return self.shape[1]

    # *****************************************************************

    @property
    def ysize(self): return self.shape[0]

    # *****************************************************************

    def set_unit(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        self.unit = unit

    # *****************************************************************

    def set_fwhm(self, fwhm):

        """
        This function ...
        :param fwhm:
        :return:
        """

        self.fwhm = fwhm

    # *****************************************************************

    def convert_to(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        # Convert the data
        conversion_factor = self.unit / unit
        self *= conversion_factor

        # Set the new unit
        self.unit = unit

    # *****************************************************************

    def to_magnitude(self, m_0):

        """
        This function ...
        :return:
        """

        # Do the conversion
        self = m_0 - 2.5 * np.log10(self)

        # TODO: adapt the unit

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

    # *****************************************************************

    def to_flux(self, f_0):

        """
        This function ...
        :return:
        """

        # Do the conversion
        self = f_0 * np.power(10.0, - self / 2.5)

        # TODO: adapt the unit

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

     # *****************************************************************

    def convolve(self, kernel):

        """
        This function ...
        :param frame:
        :return:
        """

        # Get the pixel scale of the kernel
        pixelscale_kernel = kernel.pixelscale

        # Calculate the zooming factor
        factor = self.pixelscale / pixelscale_kernel

        # Rebin the kernel to the same grid of the image
        kernel = ndimage.interpolation.zoom(kernel, zoom=1.0/factor)

        # Do the convolution on this frame
        self = convolve_fft(self, kernel, normalize_kernel=True)

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

    # *****************************************************************

    def rebin(self, ref_frame):

        """
        This function ...
        :param ref_frame:
        :return:
        """

        # Create headers
        header = self.wcs.to_header()
        ref_header = ref_frame.wcs.to_header()

        # Do the rebinning
        self = transformations.align_and_rebin(self, header, ref_header)

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

    # *****************************************************************

    def crop(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        # Crop the frame
        self = cropping.crop_check(self, x_min, x_max, y_min, y_max)

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

        # TODO: Adapt the wcs!

    # *****************************************************************

    def downsample(self, factor):

        """
        This function ...
        :return:
        """

        self = ndimage.interpolation.zoom(self, zoom=1.0/factor)

        # TODO: Adapt the wcs!

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

    # *****************************************************************

    def rotate(self, angle):

        """
        This function ...
        :param angle:
        :return:
        """

        # Do the rotation
        ndimage.interpolation.rotate(self, angle)

        # TODO: Adapt the wcs!

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

    # *****************************************************************

    def shift(self, extent):

        """
        This function ...
        :return:
        """

        # Transform the data
        self = ndimage.interpolation.shift(self, (extent.y, extent.x))

        # TODO: Adapt the wcs!

        # Check whether this object remains of type Frame
        assert isinstance(self, Frame)

    # *****************************************************************

    def center_around(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        center = Position(x=0.5*self.xsize, y=0.5*self.ysize)
        shift = position - center

        self.shift(shift)

    # *****************************************************************

    def coordinate_range(self):

        """
        This function ...
        :return:
        """

        # Some pixel coordinates of interest.
        pixels = np.array([[0,0],[self.xsize-1, self.ysize-1]], dtype=float)

        world = self.wcs.all_pix2world(pixels, 0)  # Convert pixel coordinates to world coordinates (RA and DEC in degrees)

        #print(world)

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
        pixelscale = self.pixelscale
        pixelscale_deg = pixelscale.to("deg").value

        # Get the center pixel
        ref_pix = self.wcs.wcs.crpix
        ref_world = self.wcs.wcs.crval

        # Get the number of pixels
        size_dec_deg = self.ysize * pixelscale_deg
        size_ra_deg = self.xsize * pixelscale_deg

        # Check whether the two different ways of calculating the RA width result in approximately the same value
        assert np.isclose(ra_distance, size_ra_deg, rtol=0.02), "The coordinate system and pixel scale do not match"
        assert np.isclose(dec_distance, size_dec_deg, rtol=0.02), "The coordinate system and pixel scale do not match"

        center = coord.SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg), frame='fk5')
        ra_span = size_ra_deg * u.deg
        dec_span = size_dec_deg * u.deg

        # Return the center coordinate and the RA and DEC span
        return center, ra_span, dec_span

    # *****************************************************************

    def contains(self, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        pixel_coordinate = self.wcs.all_world2pix([[coordinate.ra.value, coordinate.dec.value]], 0)

        x = pixel_coordinate[0][0]
        y = pixel_coordinate[0][1]

        return 0.0 <= x < self.xsize and 0.0 <= y < self.ysize

    # *****************************************************************

    def replace_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set all NaN pixels to the specified value
        self[np.isnan(self)] = value

    # *****************************************************************

    def plot(self, mask=None, color=True, nan_color='black', grid=False):

        """
        This function ...
        :param mask:
        :return:
        """

        # Mask the frame with nans
        maskedimage = np.ma.array(self, mask=mask)
        frame_with_nans = maskedimage.filled(np.NaN)

        # Create a HDU from this frame with the image header
        hdu = pyfits.PrimaryHDU(frame_with_nans, self.wcs.to_header())

        # Create a figure canvas
        figure = plt.figure(figsize=(12, 12))

        # Create a figure from this frame
        plot = aplpy.FITSFigure(hdu, figure=figure)

        # Set color scale
        if color: plot.show_colorscale()
        else: plot.show_grayscale()

        # Add a color bar
        plot.add_colorbar()

        # Set color for NaN values
        plot.set_nan_color(nan_color)

        # Set grid
        if grid: plot.add_grid()

        # Show the plot on screen
        plt.show()

# *****************************************************************