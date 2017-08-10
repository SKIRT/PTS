#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.info Contains functions for getting info of images.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.stringify import tostr
from . import headers
from ..basics.coordinatesystem import CoordinateSystem
from ...core.tools import types
from ...core.basics.log import log

# -----------------------------------------------------------------

def get_image_info(image_name, frame, **kwargs):

    """
    This function ...
    :param name:
    :param frame:
    :param kwargs:
    :return:
    """

    fltr = frame.filter
    wavelength = fltr.wavelength if fltr is not None else None
    pixelscale = frame.average_pixelscale
    fwhm = frame.fwhm

    # Get filesize
    filesize = fs.file_size(frame.path).to("MB")

    # Set the info
    info = OrderedDict()
    if kwargs.pop("name", True): info["Name"] = image_name

    if kwargs.pop("path", True): info["Path"] = frame.path

    if kwargs.pop("filter", True): info["Filter"] = fltr
    if kwargs.pop("wavelength", True): info["Wavelength"] = wavelength
    if kwargs.pop("unit", True): info["Unit"] = frame.unit
    if kwargs.pop("pixelscale", True): info["Pixelscale"] = pixelscale
    if kwargs.pop("psf_filter", True): info["PSF filter"] = frame.psf_filter_name
    if kwargs.pop("fwhm", True): info["FWHM"] = fwhm
    #if kwargs.pop("shape", True): info["Dimensions"] = (frame.xsize, frame.ysize)
    if kwargs.pop("xsize", True): info["xsize"] = frame.xsize
    if kwargs.pop("ysize", True): info["ysize"] = frame.ysize
    if kwargs.pop("filesize", True): info["File size"] = filesize

    # Return the info
    return info

# -----------------------------------------------------------------

def get_image_info_strings(image_name, frame, **kwargs):

    """
    This function ...
    :param name:
    :param frame:
    :param kwargs:
    :return:
    """

    strings = []
    info = get_image_info(image_name, frame, **kwargs)
    for name in info:
        string = name + ": " + tostr(info[name], round=True)
        strings.append(string)
    return strings

# -----------------------------------------------------------------

def get_image_info_from_header(image_name, header, **kwargs):

    """
    This function ...
    :param image_name:
    :param header:
    :param path:
    :param kwargs:
    :return:
    """

    # Get the filter and wavelength
    fltr = headers.get_filter(image_name, header)
    wavelength = fltr.wavelength if fltr is not None else headers.get_wavelength(header)

    unit = headers.get_unit(header)
    pixelscale = headers.get_pixelscale(header)
    if pixelscale is None:
        wcs = CoordinateSystem(header)
        pixelscale = wcs.average_pixelscale
    else:
        pixelscale = pixelscale.average
    fwhm = headers.get_fwhm(header)
    nxpixels = header["NAXIS1"]
    nypixels = header["NAXIS2"]

    path = kwargs.pop("image_path", False)

    # Set the info
    info = OrderedDict()
    if kwargs.pop("name", True): info["Name"] = image_name

    if path is not None and kwargs.pop("path", True): info["Path"] = path

    if kwargs.pop("filter", True): info["Filter"] = fltr
    if kwargs.pop("wavelength", True): info["Wavelength"] = wavelength
    if kwargs.pop("unit", True): info["Unit"] = unit
    if kwargs.pop("pixelscale", True): info["Pixelscale"] = pixelscale
    if kwargs.pop("fwhm", True): info["FWHM"] = fwhm
    #if kwargs.pop("shape", True): info["Dimensions"] = (nxpixels, nypixels)
    if kwargs.pop("xsize", True): info["xsize"] = nxpixels
    if kwargs.pop("ysize", True): info["ysize"] = nypixels

    if path is not None and kwargs.pop("filesize", True):
        filesize = fs.file_size(path).to("MB")
        info["File size"] = filesize

    # Return the info
    return info

# -----------------------------------------------------------------

def get_image_info_strings_from_header(image_name, header, **kwargs):

    """
    This function ...
    :param name:
    :param header:
    :return:
    """

    strings = []
    info = get_image_info_from_header(image_name, header, **kwargs)
    for name in info:
        string = name + ": " + tostr(info[name], round=True)
        strings.append(string)
    return strings

# -----------------------------------------------------------------
