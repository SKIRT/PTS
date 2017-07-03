#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.units.helper Contains helper functions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .parsing import parse_quantity

# -----------------------------------------------------------------

def is_mass_density(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "mass density"

# -----------------------------------------------------------------

def is_length(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "length"

# -----------------------------------------------------------------

def is_mass(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "mass"

# -----------------------------------------------------------------

def is_photometric_density(argument):

    """
    This function ...
    :return:
    """

    qty = parse_quantity(argument, density=True)
    return qty.unit.is_spectral_density

# -----------------------------------------------------------------

def is_angle(argument):

    """
    This function ...
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "angle"

# -----------------------------------------------------------------

def is_temperature(argument):

    """
    This function ...
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "temperature"

# -----------------------------------------------------------------
