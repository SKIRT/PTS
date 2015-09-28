#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np

# Import Astromagic modules
from ..tools import coordinates

# Import astronomical modules
import astropy.units as u
import astropy.coordinates as coord

# *****************************************************************

class Frame(np.ndarray):

    """
    This class ...
    """

    # *****************************************************************

    def __new__(cls, data, wcs=None, pixelscale=None, description=None, selected=False):

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