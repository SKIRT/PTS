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
