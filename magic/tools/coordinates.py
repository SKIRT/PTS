#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.coordinates Contains functions for dealing with sky coordinates.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.coordinates import SkyCoord

# -----------------------------------------------------------------

def pixel_mapping(wcs1, wcs2):

    """
    This function determines the mapping from pixel coordinates in header1 to pixel coordinates in header2 (the
    reference header). It takes the following arguments:
    :param wcs1:
    :param wcs2:
    :return: a NumPy array describing a grid of y,x pixel locations in the input header's pixel units but the output
    header's world units. It raises a TypeError if neither header is not a Header or WCS instance, and a
    NotImplementedError if the CTYPE in the header is not recognized.
    """

    # Convert the coordinates
    if not all([w1==w2 for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
        allowed_coords = ('GLON','GLAT','RA','DEC')
        if all([(any(word in w1 for word in allowed_coords) and
                 any(word in w2 for word in allowed_coords))
                for w1,w2 in zip(wcs1.wcs.ctype,wcs2.wcs.ctype)]):
            csys1 = wcs1.csys
            csys2 = wcs2.csys
            convert_coordinates = True
        else:
            # do unit conversions
            raise NotImplementedError("Unit conversions between {0} and {1} have not yet been implemented.".format(wcs1.wcs.ctype,wcs2.wcs.ctype))

    else: convert_coordinates = False

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
        C2 = SkyCoord(lon2, lat2, unit="deg", frame=csys2)
        C1 = C2.transform_to(csys1)
        lon2,lat2 = C1.spherical.lon.deg,C1.spherical.lat.deg

    xx1,yy1 = wcs1.wcs_world2pix(lon2, lat2, 0)
    grid = np.array([yy1.reshape(outshape),xx1.reshape(outshape)])

    # Return the grid
    return grid

# -----------------------------------------------------------------

def relative_coordinate(x, y, x_shift, y_shift):

    """
    This function ...
    :param x:
    :param y:
    :param x_delta:
    :param y_delta:
    :return:
    """

    rel_x = x - x_shift
    rel_y = y - y_shift

    return (rel_x, rel_y)

# -----------------------------------------------------------------

def absolute_coordinate(x, y, x_shift, y_shift):

    """
    This function ...
    :param x:
    :param y:
    :param x_delta:
    :param y_delta:
    :return:
    """

    abs_x = x + x_shift
    abs_y = y + y_shift

    return (abs_x, abs_y)

# -----------------------------------------------------------------

def distance_points(x_pos1, y_pos1, x_pos2, y_pos2):

    """
    This function ...
    :param x_pos1:
    :param y_pos1:
    :param x_pos2:
    :param y_pos2:
    :return:
    """

    diff_x = x_pos1 - x_pos2
    diff_y = y_pos1 - y_pos2

    return np.sqrt(diff_x**2 + diff_y**2)

# -----------------------------------------------------------------

def distance_model_point(model, x_pos, y_pos):

    """
    This function ...
    :param model:
    :param x_pos:
    :param y_pos:
    :return:
    """

    return distance_points(model.x_mean.value, model.y_mean.value, x_pos, y_pos)

# -----------------------------------------------------------------

def distance_models(model_a, model_b):

    """
    This function ...
    :param model_a:
    :param model_b:
    :return:
    """

    if type(model_a).__name__ == "Gaussian2D":

        x_mean_a = model_a.x_mean.value
        y_mean_a = model_a.y_mean.value

    elif type(model_a).__name__ == "AiryDisk2D":

        x_mean_a = model_a.x_0.value
        y_mean_a = model_a.y_0.value

    else: raise ValueError("Models other than Gaussian2D or AiryDisk2D are not yet supported")

    if type(model_b).__name__ == "Gaussian2D":

        x_mean_b = model_b.x_mean.value
        y_mean_b = model_b.y_mean.value

    elif type(model_b).__name__ == "AiryDisk2D":

        x_mean_b = model_b.x_0.value
        y_mean_b = model_b.y_0.value

    return distance_points(x_mean_a, y_mean_a, x_mean_b, y_mean_b)

# -----------------------------------------------------------------

def ra_distance(declination, ra_a, ra_b):

    """
    This function ...
    :param declination: dec in degrees
    :param ra_a: ra in degrees
    :param ra_b: ra in degrees
    :return:
    """

    cos_ra_distance = np.sin(np.radians(declination))**2 + np.cos(np.radians(declination))**2 * np.cos(np.radians(ra_b-ra_a))

    if cos_ra_distance > 1.0 and np.isclose(cos_ra_distance, 1.0): cos_ra_distance = 1.0 # Avoid crashes of np.arcos

    # Return ...
    return np.degrees(np.arccos(cos_ra_distance))

# -----------------------------------------------------------------

def ra_around(ra_center, ra_distance, declination):

    """
    This function ...
    :param ra_center: 
    :param ra_distance: 
    :param declination: 
    :return: 
    """

    difference = abs(ra_difference(ra_distance, declination))
    return ra_center - difference, ra_center + difference

# -----------------------------------------------------------------

def ra_difference(ra_distance, declination):

    """
    This function ...
    :param ra_distance:
    :param declination
    :return: 
    """

    cos_ra_difference = ( np.cos(np.radians(ra_distance)) - np.sin(np.radians(declination))**2 ) / np.cos(np.radians(declination))**2

    if cos_ra_difference > 1.0 and np.isclose(cos_ra_difference, 1.0): cos_ra_difference = 1.0

    return np.degrees(np.arccos(cos_ra_difference))

# -----------------------------------------------------------------

def degrees_to_hms(ra='', dec='', round=False, separator=":"):

    """
    This function ...
    :param ra:
    :param dec:
    :param round:
    :param separator:
    :return:
    """

    RA, DEC, rs, ds = '', '', '', ''

    if dec:

        if str(dec)[0] == '-':
            ds, dec = '-', abs(dec)
        deg = int(dec)
        decM = abs(int((dec-deg)*60))
        if round:
            decS = int((abs((dec-deg)*60)-decM)*60)
        else:
            decS = (abs((dec-deg)*60)-decM)*60
        DEC = ('{0}{1}' + separator + '{2}' + separator + '{3}').format(ds, deg, decM, decS)

    if ra:

        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra/15)
        raM = int(((ra/15)-raH)*60)
        if round:
            raS = int(((((ra/15)-raH)*60)-raM)*60)
        else:
            raS = ((((ra/15)-raH)*60)-raM)*60
        RA = ('{0}{1}' + separator + '{2}' + separator + '{3}').format(rs, raH, raM, raS)

    # Return ...
    if ra and dec: return (RA, DEC)
    else: return RA or DEC

# -----------------------------------------------------------------

def hms_to_degrees(ra='', dec=''):

    """
    This function ...
    :param ra:
    :param dec:
    :return:
    """

    RA, DEC, rs, ds = '', '', 1, 1

    if dec:

        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)

    if ra:

        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)

    if ra and dec: return (RA, DEC)
    else: return RA or DEC

# -----------------------------------------------------------------
