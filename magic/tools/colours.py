#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.colours Provides functions related to colours.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ...core.filter.filter import parse_filter
from ...core.tools import strings
from ..core.list import check_uniformity
from ...core.filter.filter import represent_filter

# -----------------------------------------------------------------

def calculate_colour(flux_a, flux_b):

    """
    This function ...
    :param flux_a:
    :param flux_b:
    :return:
    """

    # Check units
    flux_a = flux_a.to("Jy").value
    flux_b = flux_b.to("Jy").value

    # Calculate and return
    return -2.5 * np.log10(flux_a / flux_b)

# -----------------------------------------------------------------

def get_filters_for_colour(colour, delimiter="auto"):

    """
    This function ...
    :param colour:
    :param delimiter:
    :return: 
    """

    # Find the delimiter
    if delimiter == "auto": delimiter = find_delimiter(colour)

    str_a, str_b = colour.split(delimiter)
    fltr_a, fltr_b = parse_filter(str_a), parse_filter(str_b)
    return fltr_a, fltr_b

# -----------------------------------------------------------------

def same_colour(colour_a, colour_b):

    """
    This function ...
    :param colour_a:
    :param colour_b:
    :return:
    """

    # Get the filters
    fltr_aa, fltr_ab = get_filters_for_colour(colour_a)
    fltr_ba, fltr_bb = get_filters_for_colour(colour_b)

    #print(fltr_aa, fltr_ba)
    #print(fltr_ab, fltr_bb)

    # Return
    return fltr_aa == fltr_ba and fltr_ab == fltr_bb

# -----------------------------------------------------------------

def get_wavelengths_for_colour(colour, delimiter="auto"):

    """
    This function ...
    :param colour:
    :param delimiter:
    :return:
    """

    # Get filters
    fltr_a, fltr_b = get_filters_for_colour(colour, delimiter=delimiter)

    # Return wavelengths
    return fltr_a.wavelength, fltr_b.wavelength

# -----------------------------------------------------------------

def is_fir_colour(name):

    """
    This function ...
    :param name:
    :return:
    """

    from .wavelengths import is_fir
    wav_a, wav_b = get_wavelengths_for_colour(name)
    return is_fir(wav_a) and is_fir(wav_b)

# -----------------------------------------------------------------

def is_fir_or_submm_colour(name):

    """
    This function ...
    :param name:
    :return:
    """

    from .wavelengths import is_fir, is_submm
    wav_a, wav_b = get_wavelengths_for_colour(name)

    if not (is_fir(wav_a) or is_submm(wav_a)): return False
    if not (is_fir(wav_b) or is_submm(wav_b)): return False
    return True

# -----------------------------------------------------------------

def is_optical_colour(name):

    """
    This function ...
    :param name:
    :return:
    """

    from .wavelengths import is_optical
    wav_a, wav_b = get_wavelengths_for_colour(name)
    return is_optical(wav_a) and is_optical(wav_b)

# -----------------------------------------------------------------

def is_uv_optical_colour(name):

    """
    This function ...
    :param name:
    :return:
    """

    from .wavelengths import is_uv, is_optical
    wav_a, wav_b = get_wavelengths_for_colour(name)
    return is_uv(wav_a) and is_optical(wav_b)

# -----------------------------------------------------------------

def is_uv_nir_colour(name):

    """
    This function ...
    :param name:
    :return:
    """

    from .wavelengths import is_uv, is_nir
    wav_a, wav_b = get_wavelengths_for_colour(name)
    return is_uv(wav_a) and is_nir(wav_b)

# -----------------------------------------------------------------

def get_colour_name_for_filters(filter_a, filter_b, delimiter="-", filter_delimiter=" ", short=False):

    """
    This function ...
    :param filter_a:
    :param filter_b:
    :param delimiter:
    :param filter_delimiter:
    :param short:
    :return:
    """

    name_a = represent_filter(filter_a, delimiter=filter_delimiter, short=short)
    name_b = represent_filter(filter_b, delimiter=filter_delimiter, short=short)
    return name_a + delimiter + name_b

# -----------------------------------------------------------------

def get_colour_name_for_colour(colour, delimiter="-", filter_delimiter=" ", short=False):

    """
    Thisf unction ...
    :param colour:
    :param delimiter:
    :param filter_delimiter:
    :param short:
    :return:
    """

    fltr_a, fltr_b = get_filters_for_colour(colour)
    return get_colour_name_for_filters(fltr_a, fltr_b, delimiter=delimiter, filter_delimiter=filter_delimiter, short=short)

# -----------------------------------------------------------------

def make_colour_map(frame_a, frame_b):

    """
    This function ...
    :param frame_a:
    :param frame_b:
    :return:
    """

    from ...core.basics.log import log
    from ...core.tools.stringify import tostr

    # Check uniformity, THE UNIT HAS TO BE THE SAME OBVIOUSLY
    unit, wcs, pixelscale, psf_filter, fwhm, distance = check_uniformity(frame_a, frame_b)

    # Debugging
    log.debug("Both frames have a unit of " + tostr(unit, add_physical_type=True))

    # Make the colour map and return it
    # UNIT SHOULD BE NONE
    colour = Frame(-2.5 * np.log10(frame_a / frame_b), unit=None, wcs=wcs, pixelscale=pixelscale, psf_filter=psf_filter, fwhm=fwhm, distance=distance)

    # Set infinity values to Nan
    colour.replace_infs(float("nan"))

    # Return the colour map
    return colour

# -----------------------------------------------------------------

def find_delimiter(colour):

    """
    This function ...
    :param colour:
    :return:
    """

    # Find the delimiter for the colour string, so the string has to contain exact one of this delimiter!
    return strings.find_delimiter(colour, 1)

# -----------------------------------------------------------------
