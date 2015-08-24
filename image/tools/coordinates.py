#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import numpy as np

# *****************************************************************

def relative_coordinate(x, y, x_delta, y_delta):

    rel_x = x - x_delta
    rel_y = y - y_delta

    return (rel_x, rel_y)

# *****************************************************************

## This function ...
def absolute_coordinate(x, y, x_delta, y_delta):

    abs_x = x + x_delta
    abs_y = y + y_delta

    return (abs_x, abs_y)

# *****************************************************************

## This function ...
def distance_points(x_pos1, y_pos1, x_pos2, y_pos2):

    diff_x = x_pos1 - x_pos2
    diff_y = y_pos1 - y_pos2

    return np.sqrt(diff_x**2 + diff_y**2)

# *****************************************************************

## This function ...
def distance_model_point(model, x_pos, y_pos):

    return distance_points(model.x_mean.value, model.y_mean.value, x_pos, y_pos)

# *****************************************************************

## This function ...
def distance_models(model_a, model_b):

    return distance_points(model_a.x_mean.value, model_a.y_mean.value, model_b.x_mean.value, model_b.y_mean.value)

# *****************************************************************

## This function ...
def ra_distance(declination, ra_a, ra_b):

    cos_ra_distance = np.sin(np.radians(declination))**2 + np.cos(np.radians(declination))**2 * np.cos(np.radians(ra_b-ra_a))

    # Return ...
    return np.degrees(np.arccos(cos_ra_distance))

# *****************************************************************
