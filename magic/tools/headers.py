#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import other relevant PTS modules
from pts.filter import Filter

# *****************************************************************

def get_pixelscale(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the pixel scale to None
    pixelscale = None

    # Search for 'PIXSCALE' keyword
    if 'PIXSCALE' in header: pixelscale = header['PIXSCALE']

    # Search for the 'SECPIX' keyword
    elif 'SECPIX' in header: pixelscale = header['SECPIX']

    # Search for 'Pixel Field of View' keyword
    elif 'PFOV' in header: pixelscale = header['PFOV']

    # Search for the CD matrix elements
    elif 'CD1_1' in header and 'CD1_2' in header: pixelscale = math.sqrt(header['CD1_1']**2 + header['CD1_2']**2 ) * 3600.0

    # Search for the diagonal CD matrix elements
    elif 'CD1_1' in header: pixelscale = abs(header['CD1_1']) * 3600.0

    # Search for the 'CDELT1' keyword
    elif 'CDELT1' in header: pixelscale = abs(header['CDELT1']) * 3600.0

    # Return the pixel scale (in arcseconds)
    return pixelscale

# *****************************************************************

def get_filter(name, header):

    """
    This function ...
    :param name:
    :param header:
    :return:
    """

    # Initially, set the filter to None
    filter = None

    # Determine the filter from the information in the header
    if 'INSTRUME' in header and 'FILTER' in header: filterid = header['INSTRUME'].lower() + header['FILTER'].lower()

    # If no filter information could be found in the header, try to obtain it from the file name
    else: filterid = name.lower()

    # Create a filter object from the filterid
    if "fuv" in filterid: filter = Filter("GALEX.FUV")
    elif "pacs" in filterid:

        if '70' in filterid or 'blue' in filterid: filter = Filter("Pacs.blue")
        elif '100' in filterid or 'green' in filterid: filter = Filter("Pacs.green")
        elif '160' in filterid or 'red' in filterid: filter = Filter("Pacs.red")

    elif "mips" in filterid or "24" in filterid: filter = Filter("MIPS.24")
    elif "2mass" in filterid and "h" in filterid: filter = Filter("2MASS.H")
    elif "irac" in filterid:

        if '3.6' in filterid or 'i1' in filterid: filter = Filter("IRAC.I1")

    # Return the filter
    return filter

# *****************************************************************

def get_units(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the units to None
    units = None

    if 'BUNIT' in header:

        units = header['BUNIT']

    # Return the units
    return units

# *****************************************************************

def is_sky_subtracted(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the boolean to False
    subtracted = False

    if 'BACK_SUB' in header:

        subtracted = header['BACK_SUB']

    # Return the boolean value
    return subtracted

# *****************************************************************

def get_number_of_frames(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the boolean to False
    nframes = 1

    if 'NAXIS' in header:

        # If there are 3 axes, get the size of the third
        if header['NAXIS'] == 3: nframes = header['NAXIS3']

    # Return the boolean value
    return nframes

# *****************************************************************

def get_frame_description(header, i):

    """
    This function ...
    :param header:
    :param i:
    :return:
    """

    planeX = "PLANE" + str(i)

    # Get the description
    description = header[planeX]

    # Return the description
    return description

# *****************************************************************

def get_frame_name(description):

    """
    This function ...
    :param description:
    :return:
    """

    # Convert spaces to underscores and ignore things between parentheses
    name = description.split("(")[0].rstrip(" ").replace(" ", "_")

    # If the frame name contains 'error', use the standard name "errors" for this frame
    if 'error' in name:

        name = "errors"

    # Return the frame name
    return name

# *****************************************************************