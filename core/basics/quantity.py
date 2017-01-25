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
from .unit import stringify_unit

# -----------------------------------------------------------------

def parse_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try: quantity = PhotometricQuantity(argument)
    except ValueError: quantity = Quantity(argument)
    return quantity
        
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
