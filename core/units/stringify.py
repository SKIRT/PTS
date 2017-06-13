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

def stringify_unit(unit, add_physical_type=False):

    """
    This function ...
    :param unit:
    :param add_physical_type:
    :return:
    """

    from .unit import PhotometricUnit

    # Get parsing type
    if isinstance(unit, PhotometricUnit):
        if unit.density: parsing_type = "photometric_density_unit"
        else: parsing_type = "photometric_unit"
    else: parsing_type = "unit"

    # Return type and stringified unit
    return parsing_type, represent_unit(unit, add_physical_type=add_physical_type)

# -----------------------------------------------------------------

def represent_unit(unit, add_physical_type=False):

    """
    This function ...
    :param unit:
    :param add_physical_type:
    :return:
    """

    string = str(unit)

    for key in output_replacements:
        string = string.replace(key, output_replacements[key])

    # Remove whitespace
    string = string.replace(" ", "")

    # Add details if requested
    if add_physical_type: string += " [" + unit.physical_type + "]"

    # Return the string
    return string

# -----------------------------------------------------------------

def stringify_quantity(quantity, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, add_physical_type=False):

    """
    This function ...
    :param quantity:
    :param scientific:
    :param decimal_places:
    :param fancy:
    :param ndigits:
    :param unicode:
    :param add_physical_type:
    :return:
    """

    # Stringify the unit
    unit_type, unit_string = stringify_unit(quantity.unit, add_physical_type=add_physical_type)

    # Determine quantity parsing type
    if unit_type == "photometric_unit": parsing_type = "photometric_quantity"
    elif unit_type == "photometric_density_unit": parsing_type = "photometric_density_quantity"
    elif unit_type == "unit": parsing_type = "quantity"
    else: raise ValueError("Unknown unit type: " + unit_type)

    # Import function that converts real to string
    from ..tools.stringify import str_from_real

    # Return parsing type and stringified quantity
    return parsing_type, str_from_real(quantity.value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits, unicode=unicode) + " " + unit_string

# -----------------------------------------------------------------

def str_from_quantity(quantity, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, add_physical_type=False):

    """
    This function ...
    :param quantity:
    :param scientific:
    :param decimal_places:
    :param fancy:
    :param ndigits:
    :param unicode:
    :param add_physical_type:
    :param kwargs:
    :return:
    """

    ptype, string = stringify_quantity(quantity, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits, unicode=unicode, add_physical_type=add_physical_type)
    return string

# -----------------------------------------------------------------

def represent_quantity(quantity, scientific=False, decimal_places=2, add_physical_type=False):

    """
    This function ...
    :param quantity:
    :param scientific:
    :param decimal_places:
    :param add_physical_type:
    :return:
    """

    from ..tools.stringify import str_from_real
    return str_from_real(quantity.value, scientific=scientific, decimal_places=decimal_places) + " " + represent_unit(quantity.unit, add_physical_type=add_physical_type)

# -----------------------------------------------------------------

def str_from_quantity_range(the_range, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, add_physical_type=False, **kwargs):

    """
    This function ...
    :param the_range: 
    :param scientific: 
    :param decimal_places: 
    :param fancy: 
    :param ndigits: 
    :param unicode:
    :param add_physical_type:
    :param kwargs:
    :return: 
    """

    min_str = str_from_quantity(the_range.min, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits, unicode=unicode, add_physical_type=add_physical_type)
    max_str = str_from_quantity(the_range.max, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits, unicode=unicode, add_physical_type=add_physical_type)

    return min_str + " > " + max_str

# -----------------------------------------------------------------
