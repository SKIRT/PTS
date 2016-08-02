#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.mask Contains the Mask class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.io import fits

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..basics.mask import MaskBase

# -----------------------------------------------------------------

class Mask(MaskBase):
    
    """
    This class ...
    """

    def __init__(self, data, **kwargs):

        """
        The constructor ...
        :param data:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Mask, self).__init__(data, **kwargs)

        # Set the WCS
        self.wcs = kwargs.pop("wcs", None)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, index=None, plane=None, hdulist_index=None):

        """
        This function ...
        :param path:
        :param index:
        :param plane:
        :param hdulist_index:
        :return:
        """

        name = None
        description = None
        no_filter = True
        fwhm = None
        add_meta = False

        from . import io  # Import here because io imports SegmentationMap

        # PASS CLS TO ENSURE THIS CLASSMETHOD WORKS FOR ENHERITED CLASSES!!
        return io.load_frame(cls, path, index, name, description, plane, hdulist_index, no_filter, fwhm, add_meta=add_meta)

    # -----------------------------------------------------------------

    @property
    def header(self):

        """
        This function ...
        :return:
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

    def save(self, path, header=None):

        """
        This function ...
        :param path:
        :param header:
        :return:
        """

        # If a header is not specified, created it from the WCS
        if header is None: header = self.header

        from . import io  # Import here because io imports SegmentationMap

        # Write to a FITS file
        io.write_frame(self._data.astype(int), header, path)

# -----------------------------------------------------------------
