#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import numpy as np

# *****************************************************************

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

# *****************************************************************

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

# *****************************************************************

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

# *****************************************************************

def distance_model_point(model, x_pos, y_pos):

    """
    This function ...
    :param model:
    :param x_pos:
    :param y_pos:
    :return:
    """

    return distance_points(model.x_mean.value, model.y_mean.value, x_pos, y_pos)

# *****************************************************************

def distance_models(model_a, model_b):

    """
    This function ...
    :param model_a:
    :param model_b:
    :return:
    """

    return distance_points(model_a.x_mean.value, model_a.y_mean.value, model_b.x_mean.value, model_b.y_mean.value)

# *****************************************************************

def ra_distance(declination, ra_a, ra_b):

    """
    This function ...
    :param declination:
    :param ra_a:
    :param ra_b:
    :return:
    """

    cos_ra_distance = np.sin(np.radians(declination))**2 + np.cos(np.radians(declination))**2 * np.cos(np.radians(ra_b-ra_a))

    # Return ...
    return np.degrees(np.arccos(cos_ra_distance))

# *****************************************************************

def degrees_to_hms(ra='', dec='', round=False):

    """
    This function ...
    :param ra:
    :param dec:
    :param round:
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
        DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)

    if ra:

        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra/15)
        raM = int(((ra/15)-raH)*60)
        if round:
            raS = int(((((ra/15)-raH)*60)-raM)*60)
        else:
            raS = ((((ra/15)-raH)*60)-raM)*60
        RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)

    # Return ...
    if ra and dec: return (RA, DEC)
    else: return RA or DEC

# *****************************************************************

def hms_to_degrees(ra='', dec=''):

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

# *****************************************************************