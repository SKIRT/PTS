#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import math

# Import PTS modules
from pts.filter import Filter

# Import astronomical modules
import astropy.wcs as pywcs
from astropy import coordinates
from astropy import units as u
from astropy import log

# -----------------------------------------------------------------

def get_pixelscale(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Search for 'PIXSCALE' keyword
    if 'PIXSCALE' in header: return header['PIXSCALE'] * u.arcsec

    # Search for the 'SECPIX' keyword
    elif 'SECPIX' in header: return header['SECPIX'] * u.arcsec

    # Search for 'Pixel Field of View' keyword
    elif 'PFOV' in header: return header['PFOV'] * u.arcsec

    # Search for the CD matrix elements
    elif 'CD1_1' in header and 'CD1_2' in header: return math.sqrt(header['CD1_1']**2 + header['CD1_2']**2 ) * 3600.0 * u.arcsec

    # Search for the diagonal CD matrix elements
    elif 'CD1_1' in header: return abs(header['CD1_1']) * 3600.0 * u.arcsec

    # Search for the 'CDELT1' keyword
    elif 'CDELT1' in header: return abs(header['CDELT1']) * 3600.0 * u.arcsec

    # If none of the above keywords were found, return None
    else: return None

# -----------------------------------------------------------------

def get_filter(name, header):

    """
    This function ...
    :param name:
    :param header:
    :return:
    """

    filterid = name.lower()

    # Determine the filter from the information in the header
    if 'INSTRUME' in header: filterid += header['INSTRUME'].lower()
    if 'FILTER' in header: filterid += header['FILTER'].lower()
    if 'FLTRNM' in header: filterid += header['FLTRNM'].lower()

    # Create a filter object from the filterid
    if "fuv" in filterid: return Filter("GALEX.FUV")
    elif "pacs" in filterid:

        if '70' in filterid or 'blue' in filterid: return Filter("Pacs.blue")
        elif '100' in filterid or 'green' in filterid: return Filter("Pacs.green")
        elif '160' in filterid or 'red' in filterid: return Filter("Pacs.red")
        else:
            log.warning("Could not determine which PACS filter was used for this image")
            return None

    elif "mips" in filterid:

        if "24" in filterid: return Filter("MIPS.24")
        elif "70" in filterid: return Filter("MIPS.70")
        elif "160" in filterid: return Filter("MIPS.160")
        else:
            log.warning("Could not determine which MIPS filter was used for this image")
            return None

    elif "2mass" in filterid and "h" in filterid: return Filter("2MASS.H")
    elif "2mass" in filterid and "j" in filterid: return Filter("2MASS.J")
    elif "2mass" in filterid and "k" in filterid: return Filter("2MASS.Ks")
    elif "irac" in filterid:

        if '3.6' in filterid or 'i1' in filterid: return Filter("IRAC.I1")
        elif '4.5' in filterid or 'i2' in filterid: return Filter("IRAC.I2")
        elif '5.8' in filterid or 'i3' in filterid: return Filter("IRAC.I3")
        elif '8.0' in filterid or 'i4' in filterid: return Filter("IRAC.I4")
        else:
            log.warning("Could not determine which IRAC filter was used for this image")
            return None

    elif "alpha" in filterid or "6561" in filterid: return Filter("656_1")
    elif "r" in filterid and "kpno" in filterid: return Filter("KPNO.Mosaic.R")

    # The filter could not be determined from the specified header
    else: return None

# -----------------------------------------------------------------

def get_unit(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Look for the 'BUNIT' keyword
    if "BUNIT" in header:

        value = header["BUNIT"].split("   / ")[0].rstrip()
        return u.Unit(value)

    # Look for the 'SIGUNIT' keyword
    elif "SIGUNIT" in header:

        value = header["SIGUNIT"].split("   / ")[0].rstrip()
        return u.Unit(value)

    # Look for the 'ZUNITS' keyword
    elif "ZUNITS" in header:

        value = header["ZUNITS"].split("   / ")[0].rstrip()
        return u.Unit(value)

    # No unit information was found
    else: return None

# -----------------------------------------------------------------

def is_sky_subtracted(header):

    """
    This function ...
    :param header:
    :return:
    """

    # Initially, set the boolean to False
    subtracted = False

    if 'BACK_SUB' in header: subtracted = header['BACK_SUB']

    # Return the boolean value
    return subtracted

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

def get_frame_index(header, name):

    """
    This function ...
    """

    for key in header:

        # Skip keys not ...
        if not "PLANE" in key: continue

        if header[key] == name: return int(key.split("PLANE")[1])

    return None

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------
