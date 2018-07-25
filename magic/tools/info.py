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

# Import astronomical modules
from astropy.io.fits import getheader

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.stringify import tostr
from . import headers
from ..basics.coordinatesystem import CoordinateSystem
from ...core.filter.filter import parse_filter
from ...core.units.unit import parse_unit as u
from ...core.units.unit import parse_quantity
from ...core.tools import strings

# -----------------------------------------------------------------

def get_image_info_file(image_name, frame_path, **kwargs):

    """
    This function ...
    :param image_name:
    :param frame_path:
    :param kwargs:
    :return:
    """

    from ..core.frame import Frame
    frame = Frame.from_file(frame_path)
    return get_image_info(image_name, frame, **kwargs)

# -----------------------------------------------------------------

def get_image_info(image_name, frame, **kwargs):

    """
    This function ...
    :param image_name:
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

def get_image_info_from_remote_header_file(image_name, frame_path, session_or_remote, **kwargs):

    """
    This function ...
    :param image_name:
    :param frame_path:
    :param session:
    :param kwargs:
    :return:
    """

    from ...core.remote.remote import Remote
    from ...core.remote.python import RemotePythonSession

    # Make session
    if isinstance(session_or_remote, Remote): session = session_or_remote.start_python_session(output_path=session_or_remote.pts_temp_path, attached=True)
    elif isinstance(session_or_remote, RemotePythonSession): session = session_or_remote
    else: raise ValueError("Invalid value for 'session_or_remote': " + str(session_or_remote))

    # Import
    session.import_package("getheader", from_name="astropy.io.fits")
    session.import_package("headers", from_name="pts.magic.tools")

    # Load header
    session.send_line("header = getheader('" + frame_path + "')")

    # Get filter
    session.send_line("fltr = headers.get_filter('" + image_name + "', header)")
    filter_name = session.get_simple_variable("str(fltr)")

    # Get wavelength
    if strings.unquote(filter_name) == "None": fltr = wavelength = None
    else:
        fltr = parse_filter(filter_name)
        wavelength = fltr.wavelength

    # Get unit
    session.send_line("unit = headers.get_unit(header)")
    unit_string = session.get_simple_variable("str(unit)")
    if strings.unquote(unit_string) == "None": unit = None
    else: unit = u(unit_string)

    # Get npixels
    nxpixels = session.get_simple_variable('header["NAXIS1"]')
    nypixels = session.get_simple_variable('header["NAXIS2"]')

    # Get the pixelscale
    session.send_line('pixelscale = headers.get_pixelscale(header)')
    pixelscale_string = session.get_simple_variable('str(pixelscale)')
    if strings.unquote(pixelscale_string) == "None": pixelscale = None
    else: pixelscale = parse_quantity(pixelscale_string)

    if pixelscale is None:
        session.import_package("CoordinateSystem", from_name="pts.magic.basics.coordinatesystem")
        session.send_line("wcs = CoordinateSystem(header)")
        pixelscale_string = session.get_simple_variable('str(wcs.average_pixelscale)')
        pixelscale = parse_quantity(pixelscale_string)

    # Get the FWHM
    session.send_line('fwhm = headers.get_fwhm(header)')
    fwhm_string = session.get_simple_variable('str(fwhm)')
    if strings.unquote(fwhm_string) == "None": fwhm = None
    else: fwhm = parse_quantity(fwhm_string)

    session.send_line("psf_filter = headers.get_psf_filter(header)")
    psf_filter_string = session.get_simple_variable('str(psf_filter)')
    if strings.unquote(psf_filter_string) == "None": psf_filter_name = None
    else: psf_filter_name = psf_filter_string

    # Set the info
    info = OrderedDict()

    if kwargs.pop("name", True): info["Name"] = image_name
    if kwargs.pop("path", True): info["Path"] = frame_path
    if kwargs.pop("filter", True): info["Filter"] = fltr
    if kwargs.pop("wavelength", True): info["Wavelength"] = wavelength
    if kwargs.pop("unit", True): info["Unit"] = unit
    if kwargs.pop("pixelscale", True): info["Pixelscale"] = pixelscale
    if kwargs.pop("fwhm", True): info["FWHM"] = fwhm
    if kwargs.pop("psf_filter", True): info["PSF filter"] = psf_filter_name
    if kwargs.pop("xsize", True): info["xsize"] = nxpixels
    if kwargs.pop("ysize", True): info["ysize"] = nypixels
    if kwargs.pop("filesize", True): info["File size"] = session.file_size(frame_path)

    # Return the info
    return info

# -----------------------------------------------------------------

def get_image_info_from_header_file(image_name, frame_path, **kwargs):

    """
    This function ...
    :param image_name:
    :param frame_path:
    :param kwargs:
    :return:
    """

    header = getheader(frame_path)
    return get_image_info_from_header(image_name, header, **kwargs)

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
        # Get flattened form of the header
        flattened_header = headers.flattened(header)
        wcs = CoordinateSystem(flattened_header)
        #print(wcs)
        pixelscale = wcs.average_pixelscale
    else:
        pixelscale = pixelscale.average
    fwhm = headers.get_fwhm(header)
    nxpixels = header["NAXIS1"]
    nypixels = header["NAXIS2"]

    psf_filter = headers.get_psf_filter(header)
    psf_filter_name = str(psf_filter)

    path = kwargs.pop("image_path", None)

    # Set the info
    info = OrderedDict()
    if kwargs.pop("name", True): info["Name"] = image_name

    if path is not None and kwargs.pop("path", True): info["Path"] = path

    if kwargs.pop("filter", True): info["Filter"] = fltr
    if kwargs.pop("wavelength", True): info["Wavelength"] = wavelength
    if kwargs.pop("unit", True): info["Unit"] = unit
    if kwargs.pop("pixelscale", True): info["Pixelscale"] = pixelscale
    if kwargs.pop("fwhm", True): info["FWHM"] = fwhm
    if kwargs.pop("psf_filter", True): info["PSF filter"] = psf_filter_name
    if kwargs.pop("xsize", True): info["xsize"] = nxpixels
    if kwargs.pop("ysize", True): info["ysize"] = nypixels

    if path is not None and kwargs.pop("filesize", True):
        filesize = fs.file_size(path).to("MB")
        info["File size"] = filesize

    # Return the info
    return info

# -----------------------------------------------------------------

def get_image_info_strings_from_header_file(image_name, frame_path, **kwargs):

    """
    This function ...
    :param image_name:
    :param frame_path:
    :param kwargs:
    :return:
    """

    header = getheader(frame_path)
    return get_image_info_strings_from_header(image_name, header, **kwargs)

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
