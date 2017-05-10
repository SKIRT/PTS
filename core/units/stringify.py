#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.units Contains functions for converting units, angles and quantities to strings

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# -----------------------------------------------------------------

# Output unit->string replacements
output_replacements = OrderedDict()
output_replacements["solMass"] = "Msun"
output_replacements["solLum"] = "Lsun"

# -----------------------------------------------------------------

def stringify_unit(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    from .unit import PhotometricUnit

    # Get parsing type
    if isinstance(unit, PhotometricUnit):
        if unit.density: parsing_type = "photometric_density_unit"
        else: parsing_type = "photometric_unit"
    else: parsing_type = "unit"

    # Return type and stringified unit
    return parsing_type, represent_unit(unit)

# -----------------------------------------------------------------

def represent_unit(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    string = str(unit)

    for key in output_replacements:
        string = string.replace(key, output_replacements[key])

    # Remove whitespace
    string = string.replace(" ", "")

    return string

# -----------------------------------------------------------------

def stringify_quantity(quantity, scientific=False, decimal_places=2, fancy=False, ndigits=None):

    """
    This function ...
    :param quantity:
    :param scientific:
    :param decimal_places:
    :param fancy:
    :param ndigits:
    :return:
    """

    # Stringify the unit
    unit_type, unit_string = stringify_unit(quantity.unit)

    # Determine quantity parsing type
    if unit_type == "photometric_unit": parsing_type = "photometric_quantity"
    elif unit_type == "photometric_density_unit": parsing_type = "photometric_density_quantity"
    elif unit_type == "unit": parsing_type = "quantity"
    else: raise ValueError("Unknown unit type: " + unit_type)

    # Import function that converts real to string
    from ..tools.stringify import str_from_real

    # Return parsing type and stringified quantity
    return parsing_type, str_from_real(quantity.value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits) + " " + unit_string

# -----------------------------------------------------------------

def represent_quantity(quantity, scientific=False, decimal_places=2):

    """
    This function ...
    :param quantity:
    :param scientific:
    :param decimal_places:
    :return:
    """

    from ..tools.stringify import str_from_real
    return str_from_real(quantity.value, scientific=scientific, decimal_places=decimal_places) + " " + represent_unit(quantity.unit)

# -----------------------------------------------------------------
