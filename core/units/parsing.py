#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.units.parsing Contains unit, quantity and angle parsing functions

# -----------------------------------------------------------------

# Import astronomical modules
from astropy.units import Unit

# Import astronomical modules
from astropy.units import Quantity
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ..tools import types
from .utils import clean_unit_string

# -----------------------------------------------------------------

def parse_unit(argument, density=False, brightness=False, density_strict=False, brightness_strict=False):

    """
    This function ...
    :param argument:
    :param density:
    :param brightness:
    :param density_strict:
    :param brightness_strict:
    :return:
    """

    from .unit import PhotometricUnit

    try: unit = PhotometricUnit(argument, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)
    except ValueError:
        if types.is_string_type(argument): argument = clean_unit_string(argument)
        unit = Unit(argument)
    return unit

# -----------------------------------------------------------------

def parse_unit_or_quantity(argument, density=False, brightness=False, density_strict=False, brightness_strict=False):

    """
    This function ...
    :param argument:
    :param density:
    :param brightness:
    :param density_strict:
    :param brightness_strict:
    :return:
    """

    # Parse as quantity
    quantity = parse_quantity(argument, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)

    # Return
    if quantity.value == 1.: return quantity.unit
    else: return quantity

# -----------------------------------------------------------------

def is_photometric_unit(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from .unit import PhotometricUnit

    try:
        unit = PhotometricUnit(argument)
        return True
    except ValueError: return False

# -----------------------------------------------------------------

def parse_photometric_unit(argument, density=False, brightness=False, density_strict=False, brightness_strict=False):

    """
    This function ...
    :param argument:
    :param density:
    :param brightness:
    :param density_strict:
    :param brightness_strict:
    :return:
    """

    from .unit import PhotometricUnit

    unit = PhotometricUnit(argument, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)
    return unit

# -----------------------------------------------------------------

def split_value_unit(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    units = ""
    number = 1.0
    while argument:
        try:
            number = float(argument)
            break
        except ValueError:
            units = argument[-1:] + units
            argument = argument[:-1]

    # Return
    return number, units

# -----------------------------------------------------------------

def parse_quantity(argument, density=False, physical_type=None, brightness=False, default_unit=None,
                   density_strict=False, brightness_strict=False, allow_scalar=False):

    """
    This function ...
    :param argument:
    :param density:
    :param physical_type:
    :param brightness:
    :param default_unit: for when numerical value is zero and no unit is in the argument string
    :param density_strict:
    :param brightness_strict:
    :param allow_scalar:
    :return:
    """

    # Argument is a quantity (or derived type)
    if isinstance(argument, Quantity):

        number = argument.value
        unit = argument.unit
        unit = parse_unit(unit, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)

    # Argument is a string
    elif types.is_string_type(argument):

        # Split
        number, units = split_value_unit(argument)
        #print(number, units)

        # Check whether unit is given
        if units == "":

            # Check whether the number is zero and default unit is specified
            if number == 0 and default_unit is not None: unit = parse_unit(default_unit, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)

            #
            elif allow_scalar: return number

            # Otherwise, we cannot make a quantity
            else: raise ValueError("Unit is not specified")

        # Parse the unit
        else: unit = parse_unit(units.strip(), density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)

    # Invalid input
    else: raise ValueError("Argument must be either a string or a quantity (type is " + str(type(argument)) + ")")

    # Check physical type
    if physical_type is not None:
        if unit.physical_type != physical_type: raise ValueError("The quantity was not parsed as a quantity of '" + physical_type + "'")

    # Create quantity
    return number * unit

# -----------------------------------------------------------------

def is_photometric_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    # Argument is a quantity (or derived type)
    if isinstance(argument, Quantity):

        number = argument.value
        unit = argument.unit
        return is_photometric_unit(unit)

    # Argument is a string
    elif types.is_string_type(argument):

        # Split
        number, units = split_value_unit(argument)
        return is_photometric_unit(units)

    # Invalid input
    else: raise ValueError("Argument must be either a string or a quantity")

# -----------------------------------------------------------------

def parse_angle(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    quantity = parse_quantity(argument, physical_type="angle", default_unit="deg")
    return Angle(quantity.value, quantity.unit)

# -----------------------------------------------------------------

def possible_physical_types(argument):

    """
    This function ....
    :param argument:
    :return:
    """

    physical_types = set()

    for density in True, False:
        for brightness in True, False:

            unit = parse_unit_or_quantity(argument, density=density, brightness=brightness)
            physical_types.add(unit.physical_type)

    # Return
    return list(physical_types)

# -----------------------------------------------------------------

def possible_density_flags(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    flags = set()

    for density in True, False:
        unit = parse_unit_or_quantity(argument, density=density)
        flags.add(unit.density)

    # Return
    return list(flags)

# -----------------------------------------------------------------

def possible_brightness_flags(arguments):

    """
    This function ...
    :param arguments:
    :return:
    """

    flags = set()

    for brightness in True, False:
        unit = parse_unit_or_quantity(arguments, brightness=brightness)
        flags.add(unit.brightness)

    # Return
    return list(flags)

# -----------------------------------------------------------------
