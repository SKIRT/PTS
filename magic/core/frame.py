#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.frame Contains the Frame class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os.path
import copy
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

# Import astronomical modules
import aplpy
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.convolution import convolve_fft

# Import the relevant AstroMagic classes and modules
from . import Box
from ..basics import Position, Extent, Rectangle, CoordinateSystem
from ..tools import coordinates, cropping, transformations, interpolation, headers, fitting
from ...core.tools.logging import log

# -----------------------------------------------------------------

class Frame(np.ndarray):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __new__(cls, data, wcs=None, description=None, selected=False, unit=None, name=None, filter=None, sky_subtracted=False, zero_point=None):

        """
        This function ...
        :param cls:
        :param data:
        :param wcs:
        :return:
        """

        obj = np.asarray(data).view(cls)
        obj.wcs = wcs
        obj.description = description
        obj.selected = selected
        obj.unit = unit
        obj.name = name
        obj.filter = filter
        obj.sky_subtracted = sky_subtracted
        obj.zero_point = zero_point
        obj.fwhm = None

        return obj

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, index=0, name=None, description=None, plane=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Open the HDU list for the FITS file
        hdulist = fits.open(path)

        # Get the primary HDU
        hdu = hdulist[0]

        # Get the image header
        header = hdu.header

        # Check whether multiple planes are present in the FITS image
        nframes = headers.get_number_of_frames(header)

        # Remove references to a potential third axis
        flat_header = headers.flattened(header)

        # Obtain the world coordinate system from the 'flattened' header
        wcs = CoordinateSystem(flat_header)

        # Load the frames
        header_pixelscale = headers.get_pixelscale(header) # NOTE: SOMETIMES PLAIN WRONG IN THE HEADER !!
        pixelscale = wcs.pixelscale

        # Check whether pixelscale defined in the header is correct
        if header_pixelscale is not None:

            x_isclose = np.isclose(header_pixelscale.x.to("arcsec/pix").value, pixelscale.x.to("arcsec/pix").value)
            y_isclose = np.isclose(header_pixelscale.y.to("arcsec/pix").value, pixelscale.y.to("arcsec/pix").value)

            if not (x_isclose or y_isclose):

                print("WARNING: the pixel scale defined in the header is WRONG:")
                print("           - header pixelscale: (", header_pixelscale.x.to("arcsec/pix"), header_pixelscale.y.to("arcsec/pix"), ")")
                print("           - actual pixelscale: (", pixelscale.x.to("arcsec/pix"), pixelscale.y.to("arcsec/pix"), ")")

        # Obtain the filter for this image
        filter = headers.get_filter(os.path.basename(path[:-5]), header)

        # Obtain the units of this image
        unit = headers.get_unit(header)

        # Get the magnitude zero-point
        zero_point = headers.get_zero_point(header)

        # Check whether the image is sky-subtracted
        sky_subtracted = headers.is_sky_subtracted(header)

        if nframes > 1:

            # Get the description of this frame index
            description = headers.get_frame_description(header, index)

            if plane is not None:

                description = plane
                index = headers.get_frame_index(header, plane)

            # Get the name from the file path
            if name is None: name = os.path.basename(path[:-5])

            # Return the frame
            return cls(hdu.data[index], wcs, description, False, unit, name, filter, sky_subtracted, zero_point)

        else:

            # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
            if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

            # Get the name from the file path
            if name is None: name = os.path.basename(path[:-5])

            # Return the frame
            return cls(hdu.data, wcs, description, False, unit, name, filter, sky_subtracted, zero_point)

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.wcs.pixelscale

    # -----------------------------------------------------------------

    @property
    def xy_average_pixelscale(self):

        """
        This function ...
        :return:
        """

        x_pixelscale = abs(self.pixelscale.x.to("arcsec/pix"))
        y_pixelscale = abs(self.pixelscale.y.to("arcsec/pix"))

        if not np.isclose(x_pixelscale.value, y_pixelscale.value, rtol=0.0005):

            print("WARNING: averaging the pixelscale over the x and y direction may not be a good approximation:")
            print("          x pixelscale (absolute value) =", x_pixelscale)
            print("          y pixelscale (absolute value) =", y_pixelscale)

        # Return a single value for the pixelscale in arcseconds
        return 0.5 * (x_pixelscale + y_pixelscale)

    # -----------------------------------------------------------------

    @classmethod
    def zeros_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Return a zero-filled copy of the frame
        return np.zeros_like(frame)

    # -----------------------------------------------------------------

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        if obj is None: return
        self.wcs = getattr(obj, 'wcs', None)
        self.description = getattr(obj, 'description', None)
        self.selected = getattr(obj, 'selected', False)
        self.unit = getattr(obj, 'unit', None)
        self.name = getattr(obj, 'name', None)
        self.filter = getattr(obj, 'filter', None)
        self.sky_subtracted = getattr(obj, 'sky_subtracted', False)
        self.zero_point = getattr(obj, 'zero_point', None)
        self.fwhm = getattr(obj, 'fwhm', None)

    # -----------------------------------------------------------------

    def select(self):

        """
        This function ...
        :return:
        """

        self.selected = True

    # -----------------------------------------------------------------

    def deselect(self):

        """
        This function ...
        :return:
        """

        self.selected = False

    # -----------------------------------------------------------------

    @property
    def xsize(self): return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self): return self.shape[0]

    # -----------------------------------------------------------------

    @property
    def header(self):

        """
        This function ...
        """

        # If the WCS for this frame is defined, use it to create a header
        if self.wcs is not None: header = self.wcs.to_header()

        # Else, create a new empty header
        else: header = fits.Header()

        # Add properties to the header
        header['NAXIS'] = 2
        header['NAXIS1'] = self.xsize
        header['NAXIS2'] = self.ysize

        # Return the header
        return header

    # -----------------------------------------------------------------

    @property
    def wavelength(self):

        """
        This function ...
        :return:
        """

        # Return the pivot wavelength of the frame's filter, if defined
        if self.filter is None: return None
        else: return self.filter.effectivewavelength() * u.Unit("micron")

    # -----------------------------------------------------------------

    def set_unit(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        self.unit = unit

    # -----------------------------------------------------------------

    def set_fwhm(self, fwhm):

        """
        This function ...
        :param fwhm:
        :return:
        """

        self.fwhm = fwhm

    # -----------------------------------------------------------------

    def convert_to(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        # Convert the data
        conversion_factor = self.unit / unit
        frame = self * conversion_factor

        # Set the new unit
        frame.unit = unit

        # Return the converted frame
        return frame

    # -----------------------------------------------------------------

    def to_magnitude(self, m_0):

        """
        This function ...
        :return:
        """

        # TODO: adapt the unit

        # Do the conversion
        return m_0 - 2.5 * np.log10(self)

    # -----------------------------------------------------------------

    def to_flux(self, f_0):

        """
        This function ...
        :return:
        """

        # TODO: adapt the unit

        # Do the conversion
        return f_0 * np.power(10.0, - self / 2.5)

    # -----------------------------------------------------------------

    def convolved(self, kernel):

        """
        This function ...
        :param kernel:
        :return:
        """

        # Calculate the zooming factor
        factor = (self.xy_average_pixelscale.to("arcsec/pix").value / kernel.xy_average_pixelscale.to("arcsec/pix").value)

        # Rebin the kernel to the same grid of the image
        kernel = ndimage.interpolation.zoom(kernel, zoom=1.0/factor)

        # Do the convolution on this frame
        data = convolve_fft(self, kernel, normalize_kernel=True)

        # Return the convolved frame
        return Frame(data, self.wcs, self.description, self.selected, self.unit, self.name, self.filter, self.sky_subtracted)

    # -----------------------------------------------------------------

    def rebinned(self, reference_wcs):

        """
        This function ...
        :param reference_wcs:
        :return:
        """

        # Do the rebinning
        data = transformations.new_align_and_rebin(self, self.wcs, reference_wcs)

        # Return the rebinned frame
        return Frame(data, reference_wcs, self.description, self.selected, self.unit, self.name, self.filter, self.sky_subtracted)

    # -----------------------------------------------------------------

    def crop(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        # TODO: change the WCS

        # Crop the frame
        data = cropping.crop_check(self, x_min, x_max, y_min, y_max)

        # Return the cropped frame
        return Frame(data, None, self.description, self.selected, self.unit, self.name, self.filter, self.sky_subtracted)

    # -----------------------------------------------------------------

    def downsample(self, factor):

        """
        This function ...
        :return:
        """

        # TODO: change the WCS

        # Calculate the downsampled array
        data = ndimage.interpolation.zoom(self, zoom=1.0/factor)

        # Return the downsampled frame
        return Frame(data, None, self.description, self.selected, self.unit, self.name, self.filter, self.sky_subtracted)

    # -----------------------------------------------------------------

    def rotate(self, angle):

        """
        This function ...
        :param angle:
        :return:
        """

        # Calculate the rotated array
        #frame[np.isnan(frame)] = 0.0
        data = ndimage.interpolation.rotate(self, angle, reshape=False, order=1, mode='constant', cval=float('nan'))
        #new_frame = misc.imrotate(frame, angle, interp="bilinear")

        # Convert the wcs to header
        header = self.wcs.to_header()

        # Rotate the header (Sebastien's script)
        from ..tools import rotation
        rotated_header = rotation.rotate_header(header, angle)

        # Create the new WCS
        rotated_wcs = CoordinateSystem(rotated_header)

        # Return the rotated frame
        return Frame(data, rotated_wcs, self.description, self.selected, self.unit)

    # -----------------------------------------------------------------

    def shift(self, extent):

        """
        This function ...
        :return:
        """

        # TODO: change the WCS

        # Transform the data
        data = ndimage.interpolation.shift(self, (extent.y, extent.x))

        # Return the shifted frame
        return Frame(data, None, self.description, self.selected, self.unit)

    # -----------------------------------------------------------------

    def center_around(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        center = Position(x=0.5*self.xsize, y=0.5*self.ysize)
        shift = position - center

        # Return the shifted frame
        return self.shift(shift)

    # -----------------------------------------------------------------

    def bounding_box(self, unit="deg"):

        """
        This function ...
        :param unit:
        :return:
        """

        # Get coordinate range
        center, ra_span, dec_span = self.coordinate_range()

        ra = center.ra.to(unit).value
        dec = center.dec.to(unit).value

        ra_span = ra_span.to(unit).value
        dec_span = dec_span.to(unit).value

        # Create rectangle
        center = Position(ra, dec)
        radius = Extent(0.5 * ra_span, 0.5 * dec_span)
        box = Rectangle(center, radius)

        # Return the box
        return box

    # -----------------------------------------------------------------

    def coordinate_range(self, silent=False):

        """
        This function ...
        :return:
        """

        coor1 = self.wcs.wcs_pix2world(0.0, 0.0, 0)
        coor2 = self.wcs.wcs_pix2world(self.xsize - 1.0, self.ysize - 1.0, 0)

        #print(float(coor1[0]), float(coor1[1]))
        #print(float(coor2[0]), float(coor2[1]))

        # Some pixel coordinates of interest.
        #pixels = np.array([[0.0, 0.0], [self.xsize - 1.0, self.ysize - 1.0]])
        #world = self.wcs.all_pix2world(pixels, 0)  # Convert pixel coordinates to world coordinates (RA and DEC in degrees)
        #print(world)

        co1 = coord.SkyCoord(ra=float(coor1[0]), dec=float(coor1[1]), unit=(u.Unit("deg"), u.Unit("deg")), frame='fk5')
        co2 = coord.SkyCoord(ra=float(coor2[0]), dec=float(coor2[1]), unit=(u.Unit("deg"), u.Unit("deg")), frame='fk5')

        #print("co1=", co1.to_string('hmsdms'))
        #print("co2=", co2.to_string('hmsdms'))

        #coordinate1 = world[0]
        #coordinate2 = world[1]
        #ra_range = [coordinate2[0], coordinate1[0]]
        #dec_range = [coordinate2[1], coordinate1[1]]

        ra_range = [co1.ra.value, co2.ra.value]
        dec_range = [co1.dec.value, co2.dec.value]

        # Determine the center in RA and DEC (in degrees)
        ra_center = 0.5 * (ra_range[0] + ra_range[1])
        dec_center = 0.5 * (dec_range[0] + dec_range[1])

        # New
        dec_begin = dec_range[0]
        dec_end = dec_range[1]
        ra_begin = ra_range[0]
        ra_end = ra_range[1]

        # Calculate the actual RA and DEC distance in degrees
        ra_distance = abs(coordinates.ra_distance(dec_center, ra_begin, ra_end))
        dec_distance = abs(dec_end - dec_begin)

        # Calculate the pixel scale of this image in degrees
        x_pixelscale_deg = self.pixelscale.x.to("deg/pix").value
        y_pixelscale_deg = self.pixelscale.y.to("deg/pix").value

        # Get the center pixel
        ref_pix = self.wcs.wcs.crpix
        ref_world = self.wcs.wcs.crval

        # Get the number of pixels
        size_dec_deg = self.ysize * x_pixelscale_deg
        size_ra_deg = self.xsize * y_pixelscale_deg

        if not silent:

            # Check whether the two different ways of calculating the RA width result in approximately the same value
            #assert np.isclose(ra_distance, size_ra_deg, rtol=0.05), "The coordinate system and pixel scale do not match: ra_distance=" + str(ra_distance) + ",size_ra_deg=" + str(size_ra_deg)
            #assert np.isclose(dec_distance, size_dec_deg, rtol=0.05), "The coordinate system and pixel scale do not match: dec_distance=" + str(dec_distance) + ",size_dec_deg=" + str(size_dec_deg)

            if not np.isclose(ra_distance, size_ra_deg, rtol=0.05):
                print("ERROR: the coordinate system and pixel scale do not match: ra_distance = " + str(ra_distance) + ", size_ra_deg = " + str(size_ra_deg))
            if not np.isclose(dec_distance, size_dec_deg, rtol=0.05):
                print("ERROR: the coordinate system and pixel scale do not match: dec_distance = " + str(dec_distance) + ",size_dec_deg = " + str(size_dec_deg))

        center = coord.SkyCoord(ra=ra_center, dec=dec_center, unit=(u.Unit("deg"), u.Unit("deg")), frame='fk5')

        # Create RA and DEC span as quantities
        ra_span = ra_distance * u.Unit("deg")
        dec_span = dec_distance * u.Unit("deg")

        # Return the center coordinate and the RA and DEC span
        return center, ra_span, dec_span

    # -----------------------------------------------------------------

    def contains(self, coordinate, transformation="wcs"):

        """
        This function ...
        :param coordinate:
        :param transformation:
        :return:
        """

        x, y = coordinate.to_pixel(self.wcs, origin=0, mode=transformation)
        return 0.0 <= x < self.xsize and 0.0 <= y < self.ysize

    # -----------------------------------------------------------------

    def replace_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set all NaN pixels to the specified value
        self[np.isnan(self)] = value

    # -----------------------------------------------------------------

    def fit_polynomial(self, order, mask=None):

        """
        This function ...
        :return:
        """

        # Do the fitting
        polynomial = fitting.fit_polynomial(self, order, mask=mask)

        # Evaluate the polynomial
        data = fitting.evaluate_model(polynomial, 0, self.xsize, 0, self.ysize)

        # Return a new Frame
        return Frame(data, self.wcs, self.description, self.selected, self.unit, self.name, self.filter, self.sky_subtracted)

    # -----------------------------------------------------------------

    def interpolated(self, mask, method):

        """
        This function ...
        :param mask:
        :return:
        """

        # Fit a polynomial to the data
        if method == "polynomial":
            try:
                return self.fit_polynomial(3, mask=mask)
            except TypeError:
                mask = mask.eroded(2, 1)
                return self.fit_polynomial(3, mask=mask)

        # Interpolate using the local mean method
        elif method == "local_mean":

            # Calculate the interpolated data
            data = interpolation.in_paint(self, mask)

            # Return a new Frame
            return Frame(data, self.wcs, self.description, self.selected, self.unit, self.name, self.filter, self.sky_subtracted)

        # Interpolate using inverse distance weighing
        elif method == "idw":

            # Calculate the interpolated data
            data = interpolation.in_paint(self, mask, method="idw")

            # Return a new Frame
            return Frame(data, self.wcs, self.description, self.selected, self.unit, self.name, self.filter, self.sky_subtracted)

        # Calculate the mean value of the data
        elif method == "mean":

            mean = np.ma.mean(np.ma.masked_array(self, mask=mask))
            return self.full(mean)

        # Calculate the median value of the data
        elif method == "median":

            median = np.ma.median(np.ma.masked_array(self, mask=mask))
            return self.full(median)

        # Invalid option
        else: raise ValueError("Unknown interpolation method")

    # -----------------------------------------------------------------

    def box_like(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        data = self[box.y_min:box.y_max, box.x_min:box.x_max]

        # Create the new box and return it
        return Box(data, box.x_min, box.x_max, box.y_min, box.y_max)

    # -----------------------------------------------------------------

    def save(self, path, header=None, origin=None):

        """
        This function ...
        """

        if header is None: header = self.header

        # Add origin description
        if origin is not None: header["ORIGIN"] = origin
        else: header["ORIGIN"] = "Frame class of PTS package"

        # Create the HDU
        hdu = fits.PrimaryHDU(self, header)

        # Write the HDU to a FITS file
        hdu.writeto(path, clobber=True)

    # -----------------------------------------------------------------

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
        hdu = fits.PrimaryHDU(frame_with_nans, self.wcs.to_header())

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

# -----------------------------------------------------------------
