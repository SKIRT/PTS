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
    flux_a = flux_a.to("Jy")
    flux_b = flux_b.to("Jy")

    return -2.5 * np.log10(flux_a / flux_b)

# -----------------------------------------------------------------

#def make_colour_map(frame_a, frame_b):

    #"""
    #This function ...
    #:param frame_a:
    #:param frame_b:
    #:return:
    #"""

    # Check units
    #frame_a = frame_a.copy()
    #frame_a.convert_to("Jy")
    #frame_b = frame_b.copy()
    #frame_b.convert_to("Jy")

    #return Frame(-2.5 * np.log10(frame_a / frame_b), wcs=frame_a.wcs)

# -----------------------------------------------------------------

def get_filters_for_colour(colour, delimiter="auto"):

    """
    This function ...
    :param colour:
    :param delimiter:
    :return: 
    """

    if delimiter == "auto": delimiter = find_delimiter(colour)

    str_a, str_b = colour.split(delimiter)
    fltr_a, fltr_b = parse_filter(str_a), parse_filter(str_b)
    return fltr_a, fltr_b

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
    :param frame_a: HAS TO BE IN JANSKY
    :param frame_b: HAS TO BE IN JANSKY
    :return:
    """

    # Check uniformity
    unit, wcs, pixelscale, psf_filter, fwhm, distance = check_uniformity(frame_a, frame_b)

    # Make the colour map and return it
    colour = Frame(-2.5 * np.log10(frame_a / frame_b), unit=unit, wcs=wcs, pixelscale=pixelscale, psf_filter=psf_filter, fwhm=fwhm, distance=distance)

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
