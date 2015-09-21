#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import numpy as np
import math

# Import other relevant PTS modules
from pts.filter import Filter

# Import astronomical modules
import astropy.wcs as pywcs
from astropy import coordinates
from astropy import units as u

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

    planeX = "PLANE" + str(i+1)

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

def load_wcs_from_header(header):

    """
    This function ...
    :param header:
    :return:
    """

    if issubclass(pywcs.WCS, header.__class__):
        wcs = header
    else:
        try:
            wcs = pywcs.WCS(header)
        except:
            raise TypeError("header must either be a pyfits.Header or pywcs.WCS instance")

        if not hasattr(wcs,'naxis1'):
            wcs.naxis1 = header['NAXIS1']
        if not hasattr(wcs,'naxis2'):
            wcs.naxis2 = header['NAXIS2']

    return wcs

# *****************************************************************

def check_header_matches_image(image, header):

    """
    This function ...
    :param image:
    :param header:
    :return:
    """

    wcs = load_wcs_from_header(header)

    # wcs.naxis attributes are deprecated, so we perform this check conditionally
    if ((hasattr(wcs,'naxis1') and hasattr(wcs,'naxis2')) and not
            (wcs.naxis1 == image.shape[1] and wcs.naxis2 == image.shape[0])):
        raise Exception("Image shape must match header shape.")

# *****************************************************************

def get_pixel_mapping(header1, header2):

    """
    This function determines the mapping from pixel coordinates in header1 to pixel coordinates in header2 (the
    reference header). It takes the following arguments:
    :param header1:
    :param header2:
    :return: a NumPy array describing a grid of y,x pixel locations in the input header's pixel units but the output
    header's world units. It raises a TypeError if neither header is not a Header or WCS instance, and a
    NotImplementedError if the CTYPE in the header is not recognized.
    """

    # Get the WCS from the two headers
    wcs1 = load_wcs_from_header(header1)
    wcs2 = load_wcs_from_header(header2)

    # Convert the coordinates
    if not all([w1==w2 for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
        allowed_coords = ('GLON','GLAT','RA','DEC')
        if all([(any(word in w1 for word in allowed_coords) and
                 any(word in w2 for word in allowed_coords))
                for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
            csys1 = ctype_to_csys(wcs1.wcs)
            csys2 = ctype_to_csys(wcs2.wcs)
            convert_coordinates = True
        else:
            # do unit conversions
            raise NotImplementedError("Unit conversions between {0} and {1} have not yet been implemented.".format(wcs1.wcs.ctype,wcs2.wcs.ctype))
    else:
        convert_coordinates = False

    # sigh... why does numpy use matrix convention?  Makes everything so
    # much harder...
    # WCS has naxis attributes because it is loaded with
    # _load_wcs_from_header
    outshape = [wcs2.naxis2,wcs2.naxis1]

    yy2,xx2 = np.indices(outshape)

    # get the world coordinates of the output image
    lon2,lat2 = wcs2.wcs_pix2world(xx2, yy2, 0)

    # Alternative
    #x = np.arange(wcs2.naxis1)
    #y = np.arange(wcs2.naxis2)
    #X, Y = np.meshgrid(x, y)
    #lon2, lat2 = wcs2.wcs_pix2world(X, Y, 0)

    if convert_coordinates:

        # Transform the world coordinates from the output image into the coordinate
        # system of the input image
        C2 = coordinates.SkyCoord(lon2,lat2,unit=(u.deg,u.deg),frame=csys2)
        C1 = C2.transform_to(csys1)
        lon2,lat2 = C1.spherical.lon.deg,C1.spherical.lat.deg

    xx1,yy1 = wcs1.wcs_world2pix(lon2, lat2, 0)
    grid = np.array([yy1.reshape(outshape),xx1.reshape(outshape)])

    # Return the grid
    return grid

# *****************************************************************

def ctype_to_csys(wcs):

    """
    This function ...
    :param wcs:
    :return:
    """

    ctype = wcs.ctype[0]
    if 'RA' in ctype or 'DEC' in ctype:
        if wcs.equinox == 2000:
            return 'fk5'
        elif wcs.equinox == 1950:
            return 'fk4'
        else:
            raise NotImplementedError("Non-fk4/fk5 equinoxes are not allowed")
    elif 'GLON' in ctype or 'GLAT' in ctype:
        return 'galactic'

# *****************************************************************