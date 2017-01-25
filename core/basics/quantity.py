#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.quantity Contains the PhotometricQuantity class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from .unit import stringify_unit, represent_unit, parse_unit

# -----------------------------------------------------------------

def parse_quantity(argument, density=False):

    """
    This function ...
    :param argument:
    :param density:
    :return:
    """

    # NEW IMPLEMENTATION
    units = ""
    number = 1.0
    while argument:
        try:
            number = float(argument)
            break
        except ValueError:
            units = argument[-1:] + units
            argument = argument[:-1]
    if units == "": raise ValueError("Unit is not specified")

    # Create quantity
    return number * parse_unit(units.strip(), density=density)

# -----------------------------------------------------------------

def stringify_quantity(quantity):

    """
    This function ...
    :param quantity:
    :return:
    """

    # Stringify the unit
    unit_type, unit_string = stringify_unit(quantity.unit)

    # Determine quantity parsing type
    if unit_type == "photometric_unit": parsing_type = "photometric_quantity"
    elif unit_type == "photometric_density_unit": parsing_type = "photometric_density_quantity"
    elif unit_type == "unit": parsing_type = "quantity"
    else: raise ValueError("Unknown unit type: " + unit_type)

    # Return parsing type and stringified quantity
    return parsing_type, repr(quantity.value) + " " + unit_string

# -----------------------------------------------------------------

def represent_quantity(quantity):

    """
    This function ...
    :param quantity:
    :return:
    """

    return repr(quantity.value) + " " + represent_unit(quantity.unit)

# -----------------------------------------------------------------

# NOT OPERATIONAL YET
class PhotometricQuantity(Quantity):
    
    """
    This class ...
    """
    
    def __new__(cls, value, unit=None, dtype=None, copy=True, order=None,
                subok=False, ndmin=0):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        return super(PhotometricQuantity, cls).__new__(value, unit, dtype, copy, order, subok, ndmin)

# -----------------------------------------------------------------
