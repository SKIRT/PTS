#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.frame Contains the Frame class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import numpy as np
from scipy import ndimage

# Import astronomical modules
from reproject import reproject_exact, reproject_interp
from astropy.io import fits
from astropy.units import Unit
from astropy.convolution import convolve_fft
from astropy.nddata import NDDataArray

# Import the relevant PTS classes and modules
from .box import Box
from ..basics.vector import Position, Extent
from ..basics.geometry import Rectangle
from ..basics.skygeometry import SkyCoordinate
from ..basics.coordinatesystem import CoordinateSystem
from ..tools import coordinates, cropping, transformations, headers
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ..basics.mask import Mask

# -----------------------------------------------------------------

def get_frame_names(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Load the header
    header = fits.getheader(path)

    # Get the number of planes
    nplanes = headers.get_number_of_frames(header)

    # Initialize a dictionary to contain the frame names and corresponding descriptions
    frames = dict()

    # Look at the properties of each plane
    for i in range(nplanes):

        # Get name and description of plane
        name, description, plane_type = headers.get_frame_name_and_description(header, i, always_call_first_primary=False)
        if plane_type == "frame": frames[name] = description

    # Return the frames with their name and description
    return frames

# -----------------------------------------------------------------

class Frame(NDDataArray):

    """
    This class ...
    """

    def __init__(self, data, *args, **kwargs):

        """
        This function ...
        :param data:
        :param kwargs:
        """

        wcs = kwargs.pop("wcs", None)
        self.name = kwargs.pop("name", None)
        self.description = kwargs.pop("description", None)
        self.unit = kwargs.pop("unit", None)
        self.zero_point = kwargs.pop("zero_point", None)
        self.filter = kwargs.pop("filter", None)
        self.sky_subtracted = kwargs.pop("sky_subtracted", False)
        self.fwhm = kwargs.pop("fwhm", None)
        self._pixelscale = kwargs.pop("pixelscale", None)
        self._wavelength = kwargs.pop("wavelength", None)

        # Call the constructor of the base class
        super(Frame, self).__init__(data, *args, **kwargs)

        # Set the WCS
        self._wcs = wcs

    # -----------------------------------------------------------------

    def __getitem__(self, item):

        """
        This function ...
        :param item:
        :return:
        """

        return self._data[item]

    # -----------------------------------------------------------------

    def __setitem__(self, item, value):

        """
        This function ...
        :param item:
        :return:
        """

        self._data[item] = value

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, index=None, name=None, description=None, plane=None, hdulist_index=None, no_filter=False):

        """
        This function ...
        :param path:
        :param index:
        :param name:
        :param description:
        :param plane:
        :param hdulist_index: if None, is automatically decided based on where the imageHDU is.
        :param no_filter:
        :return:
        """

        # Show which image we are importing
        log.info("Reading in file " + path + " ...")

        # Open the HDU list for the FITS file
        hdulist = fits.open(path)

        # Look for the first HDU with data
        if hdulist_index is None:

            index = 0
            while True:

                if hdulist[index].data is not None:
                    hdulist_index = index
                    break
                index += 1

            if hdulist_index is None: raise ValueError("The FITS file does not contain any data")

        # Get the primary HDU
        hdu = hdulist[hdulist_index]

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

        if no_filter: fltr = None
        else:

            # Obtain the filter for this image
            fltr = headers.get_filter(fs.name(path[:-5]), header)

        # Obtain the units of this image
        unit = headers.get_unit(header)

        # Obtain the FWHM of this image
        fwhm = headers.get_fwhm(header)

        # Get the magnitude zero-point
        zero_point = headers.get_zero_point(header)

        # Check whether the image is sky-subtracted
        sky_subtracted = headers.is_sky_subtracted(header)

        if nframes > 1:

            if plane is not None:

                for i in range(nframes):

                    # Get name and description of frame
                    name, description, plane_type = headers.get_frame_name_and_description(header, i, always_call_first_primary=False)

                    if plane == name and plane_type == "frame":
                        index = i
                        break

                # If a break is not encountered, a matching plane name is not found
                else: raise ValueError("Plane with name '" + plane + "' not found")

            elif index is not None:

                name, description, plane_type = headers.get_frame_name_and_description(header, index, always_call_first_primary=False)

            else: # index and plane is None

                for i in range(nframes):
                    # Get name and description of frame
                    name, description, plane_type = headers.get_frame_name_and_description(header, i, always_call_first_primary=False)
                    if name == "primary": index = i
                    break

                if index is None: index = 0 # if index is still None, set it to zero (take the first plane)

            # Get the name from the file path
            if name is None: name = fs.name(path[:-5])

            # Return the frame
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            return cls(hdu.data[index],
                       wcs=wcs,
                       name=name,
                       description=description,
                       unit=unit,
                       zero_point=zero_point,
                       filter=fltr,
                       sky_subtracted=sky_subtracted,
                       fwhm=fwhm)

        else:

            # Sometimes, the 2D frame is embedded in a 3D array with shape (1, xsize, ysize)
            if len(hdu.data.shape) == 3: hdu.data = hdu.data[0]

            # Get the name from the file path
            if name is None: name = fs.name(path[:-5])

            # Return the frame
            # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
            return cls(hdu.data,
                       wcs=wcs,
                       name=name,
                       description=description,
                       unit=unit,
                       zero_point=zero_point,
                       filter=fltr,
                       sky_subtracted=sky_subtracted,
                       fwhm=fwhm)

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def apply_mask(self, mask, fill=0.0):

        """
        This function ...
        :param mask:
        :param fill:
        """

        self[mask] = fill

    # -----------------------------------------------------------------

    @property
    def all_zero(self):

        """
        This function ...
        :return:
        """

        return np.all(self._data == 0)

    # -----------------------------------------------------------------

    @property
    def all_nonzero(self):

        """
        This function ...
        :return:
        """

        return not np.any(self._data == 0)

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        # Return the pixelscale of the WCS is WCS is defined
        if self.wcs is not None: return self.wcs.pixelscale
        else: return self._pixelscale  # return the pixelscale

    # -----------------------------------------------------------------

    @property
    def xy_average_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.wcs.xy_average_pixelscale

    # -----------------------------------------------------------------

    @classmethod
    def zeros_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Return a zero-filled copy of the frame
        new_frame = frame.copy()
        new_frame._data = np.zeros_like(frame._data)
        return new_frame

    # -----------------------------------------------------------------

    @classmethod
    def nans_like(cls, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Return a NaN-filled copy of the frame
        nans = cls.zeros_like(frame)
        nans.fill(np.nan)
        return nans

    # -----------------------------------------------------------------

    def is_constant(self):

        """
        This function ...
        :return:
        """

        np.nanmax(self._data) == np.nanmin(self._data)

    # -----------------------------------------------------------------

    @property
    def xsize(self): return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self): return self.shape[0]

    # -----------------------------------------------------------------

    @property
    def fwhm_pix(self):

        """
        This function ...
        :return:
        """

        return (self.fwhm / self.xy_average_pixelscale).to("pix").value if self.fwhm is not None else None

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

        # ISSUE: see bug #4592 on Astropy GitHub (WCS.to_header issue)
        # temporary fix !!
        # I don't know whether this is a good fix.. but it seems to fix it for a particular situation
        #if "PC1_1" in header:

            #if "NAXIS1" in header: header.remove("NAXIS1")
            #if "NAXIS2" in header: header.remove("NAXIS2")
            #if "CDELT1" in header: header.remove("CDELT1")
            #if "CDELT2" in header: header.remove("CDELT2")
            #header.rename_keyword("PC1_1", "CD1_1")
            #header.rename_keyword("PC2_2", "CD2_2")

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
        if self.filter is not None: return self.filter.effectivewavelength() * Unit("micron")
        else: return self._wavelength # return the wavelength (if defined, is None otherwise)

    # -----------------------------------------------------------------

    @wavelength.setter
    def wavelength(self, value):

        """
        This function ...
        :return:
        """

        self._wavelength = value

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

    def normalize(self, to=1.0):

        """
        This function ...
        :param to:
        :return:
        """

        # Calculate the sum of all the pixels
        sum = np.nansum(self)

        # Calculate the conversion factor
        factor = to / sum

        # Multiply the frame with the conversion factor
        self.__imul__(factor)

    # -----------------------------------------------------------------

    def convolved(self, kernel, allow_huge=False):

        """
        This function ...
        :param kernel:
        :param allow_huge:
        :return:
        """

        return self.copy().convolve(kernel, allow_huge)

    # -----------------------------------------------------------------

    def convolve(self, kernel, allow_huge=False):

        """
        This function ...
        :param kernel:
        :param allow_huge:
        :return:
        """

        kernel_fwhm = kernel.fwhm

        # Skip the calculation for a constant frame
        if self.is_constant():
            copy = self.copy()
            copy.fwhm = kernel_fwhm
            return copy

        nans_mask = np.isnan(self)

        # Do the convolution on this frame
        new_data = convolve_fft(self._data, kernel, normalize_kernel=True, interpolate_nan=True, allow_huge=allow_huge)

        # Put back NaNs
        new_data[nans_mask] = float("nan")

        # Replace the data and FWHM
        self._data = new_data
        self.fwhm = kernel_fwhm

    # -----------------------------------------------------------------

    def rebinned(self, reference_wcs):

        """
        This function ...
        :param reference_wcs:
        :return:
        """

        return self.copy().rebin(reference_wcs)

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=True):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :return:
        """

        # TODO: use the 'Reproject' package here: http://reproject.readthedocs.org/en/stable/

        # Rebin the data
        #new_data = transformations.new_align_and_rebin(self._data, self.wcs, reference_wcs)

        new_data = reproject_exact(self._data, self.wcs)

        # Replace the data and WCS
        self._data = new_data
        self._wcs = reference_wcs

    # -----------------------------------------------------------------

    def cropped(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        return self.copy().crop(x_min, x_max, y_min, y_max)

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

        # Crop the frame
        new_data = cropping.crop_check(self._data, x_min, x_max, y_min, y_max)

        # Change the WCS
        new_wcs = self.wcs.copy()

        # Change the center pixel position
        new_wcs.wcs.crpix[0] -= x_min
        new_wcs.wcs.crpix[1] -= y_min

        # Change the number of pixels
        new_wcs.naxis1 = x_max - x_min
        new_wcs.naxis2 = y_max - y_min

        new_wcs._naxis1 = new_wcs.naxis1
        new_wcs._naxis2 = new_wcs.naxis2

        # Check shape of data
        assert new_data.shape[1] == (x_max - x_min) and new_data.shape[0] == (y_max - y_min)

        # Replace the data and WCS
        self._data = new_data
        self._wcs = new_wcs

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return Position(self.wcs.wcs.crpix[0], self.wcs.wcs.crpix[1])

    # -----------------------------------------------------------------

    def downsampled(self, factor, order=3):

        """
        This function ...
        :return:
        """

        # Calculate the downsampled array
        data = ndimage.interpolation.zoom(self, zoom=1.0/factor, order=order)

        new_xsize = data.shape[1]
        new_ysize = data.shape[0]

        relative_center = Position(self.center.x / self.xsize, self.center.y / self.ysize)

        new_center = Position(relative_center.x * new_xsize, relative_center.y * new_ysize)

        # Make a copy of the current WCS
        #new_wcs = self.wcs.copy() ## THIS IS NOT A REAL DEEP COPY!! THIS HAS CAUSED ME BUGS THAT GAVE ME HEADACHES
        new_wcs = copy.deepcopy(self.wcs)

        # Change the center pixel position
        new_wcs.wcs.crpix[0] = new_center.x
        new_wcs.wcs.crpix[1] = new_center.y

        # Change the number of pixels
        new_wcs.naxis1 = new_xsize
        new_wcs.naxis2 = new_ysize

        new_wcs._naxis1 = new_wcs.naxis1
        new_wcs._naxis2 = new_wcs.naxis2

        # Change the pixel scale
        new_wcs.wcs.cdelt[0] *= float(self.xsize) / float(new_xsize)
        new_wcs.wcs.cdelt[1] *= float(self.ysize) / float(new_ysize)

        # Return the downsampled frame
        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        return Frame(data,
                     wcs=new_wcs,
                     name=self.name,
                     description=self.description,
                     unit=self.unit,
                     zero_point=self.zero_point,
                     filter=self.filter,
                     sky_subtracted=self.sky_subtracted,
                     fwhm=self.fwhm)

    # -----------------------------------------------------------------

    def upsampled(self, factor, integers=False):

        """
        This function ...
        :param factor:
        :param integers:
        :return:
        """

        if integers:

            # Check whether the upsampling factor is an integer or not
            if int(factor) == factor:

                data = ndimage.zoom(self, factor, order=0)

                new_xsize = data.shape[1]
                new_ysize = data.shape[0]

                relative_center = Position(self.center.x / self.xsize, self.center.y / self.ysize)

                new_center = Position(relative_center.x * new_xsize, relative_center.y * new_ysize)

                new_wcs = self.wcs.copy()
                # Change the center pixel position
                new_wcs.wcs.crpix[0] = new_center.x
                new_wcs.wcs.crpix[1] = new_center.y

                # Change the number of pixels
                new_wcs.naxis1 = new_xsize
                new_wcs.naxis2 = new_ysize

                new_wcs._naxis1 = new_wcs.naxis1
                new_wcs._naxis2 = new_wcs.naxis2

                # Change the pixel scale
                new_wcs.wcs.cdelt[0] *= float(self.xsize) / float(new_xsize)
                new_wcs.wcs.cdelt[1] *= float(self.ysize) / float(new_ysize)

                return Frame(data, wcs=new_wcs, name=self.name, description=self.description, unit=self.unit, zero_point=self.zero_point, filter=self.filter, sky_subtracted=self.sky_subtracted, fwhm=self.fwhm)

            # Upsampling factor is not an integer
            else:

                # Create a new Frame
                #frame = Frame.zeros_like(self)

                #print("self.shape", self.shape)

                frame = self.downsampled(1./factor)

                #print("frame.shape", frame.shape)

                #print("Checking indices ...")
                indices = np.unique(self)

                #print("indices:", indices)

                # Loop over the indices
                for index in list(indices):

                    #print(index)

                    index = int(index)

                    where = Mask(self == index)

                    # Rebin this mask
                    #data = transformations.new_align_and_rebin(self.masks[mask_name].astype(float), original_wcs, reference_wcs)

                    # Calculate the downsampled array
                    data = ndimage.interpolation.zoom(where.astype(float), zoom=factor)

                    upsampled_where = data > 0.5

                    frame[upsampled_where] = index

                return frame

        else: return self.downsampled(1./factor)

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
        #header = self.wcs.to_header()

        # Rotate the header (Sebastien's script)
        #from ..tools import rotation
        #rotated_header = rotation.rotate_header(header, angle)

        # Create the new WCS
        #rotated_wcs = CoordinateSystem(rotated_header)

        rotated_wcs = self.wcs.deepcopy()
        rotated_wcs.rotateCD(angle) # STILL UNTESTED

        # Return the rotated frame
        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        return Frame(data,
                     wcs=rotated_wcs,
                     name=self.name,
                     description=self.description,
                     unit=self.unit,
                     zero_point=self.zero_point,
                     filter=self.filter,
                     sky_subtracted=self.sky_subtracted,
                     fwhm=self.fwhm)

    # -----------------------------------------------------------------

    def shift(self, extent):

        """
        This function ...
        :param extent:
        :return:
        """

        # TODO: change the WCS !!!

        # Transform the data
        data = ndimage.interpolation.shift(self, (extent.y, extent.x))

        # Return the shifted frame
        # data, wcs=None, name=None, description=None, unit=None, zero_point=None, filter=None, sky_subtracted=False, fwhm=None
        return Frame(data,
                     wcs=None,
                     name=self.name,
                     description=self.description,
                     unit=self.unit,
                     zero_point=self.zero_point,
                     filter=self.filter,
                     sky_subtracted=self.sky_subtracted,
                     fwhm=self.fwhm)

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
        center, ra_span, dec_span = self.coordinate_range

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

    @property
    def corners(self):

        """
        This function ...
        :return:
        """

        # Get coordinate values
        coordinate_values = self.wcs.calc_footprint(undistort=True)

        # Initialize a list to contain the coordinates of the corners
        corners = []

        for ra_deg, dec_deg in coordinate_values:

            # Create sky coordinate
            coordinate = SkyCoordinate(ra=ra_deg, dec=dec_deg, unit="deg", frame="fk5")

            # Add the coordinate of this corner to the list
            corners.append(coordinate)

        # Return the list of coordinates
        return corners

    # -----------------------------------------------------------------

    @property
    def coordinate_range(self):

        """
        This property ...
        :return:
        """

        return self.wcs.coordinate_range

    # -----------------------------------------------------------------

    def contains(self, coordinate):

        """
        This function ...
        :param coordinate:
        :return:
        """

        pixel = coordinate.to_pixel(self.wcs)
        return 0.0 <= pixel.x < self.xsize and 0.0 <= pixel.y < self.ysize

    # -----------------------------------------------------------------

    def replace_nans(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set all NaN pixels to the specified value
        self._data[np.isnan(self._data)] = value

    # -----------------------------------------------------------------

    def replace_infs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set all inf pixels to the specified value
        self._data[np.isinf(self._data)] = value

    # -----------------------------------------------------------------

    def box_like(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        data = self._data[box.y_min:box.y_max, box.x_min:box.x_max]

        # Create the new box and return it
        return Box(data, box.x_min, box.x_max, box.y_min, box.y_max)

    # -----------------------------------------------------------------

    def save(self, path, header=None, origin=None):

        """
        This function ...
        :param path:
        :param header:
        :param origin:
        """

        if header is None: header = self.header

        # Set unit, FWHM and filter description
        if self.unit is not None: header.set("SIGUNIT", str(self.unit), "Unit of the map")
        if self.fwhm is not None: header.set("FWHM", self.fwhm.to("arcsec").value, "[arcsec] FWHM of the PSF")
        if self.filter is not None: header.set("FILTER", str(self.filter), "Filter used for this observation")

        # Add origin description
        if origin is not None: header["ORIGIN"] = origin
        else: header["ORIGIN"] = "Frame class of PTS package"

        # Create the HDU
        hdu = fits.PrimaryHDU(self._data, header)

        # Write the HDU to a FITS file
        hdu.writeto(path, clobber=True)

# -----------------------------------------------------------------
