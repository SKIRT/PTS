#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.parsing Provides useful functions for parsing strings into a variety of types.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
import warnings
from collections import OrderedDict
from decimal import Decimal
#import numpy as np

# Import the relevant PTS classes and modules
#from ..basics.range import IntegerRange, RealRange, QuantityRange # imported below to accomodate clean python installs
from . import filesystem as fs
from . import types
# Imported below to accomodate clean python installs
#from ..filter.broad import BroadBandFilter
#from ..filter.broad import identifiers as broad_band_identifiers
#from ..filter.narrow import NarrowBandFilter, generate_aliases_ranges, wavelength_range_for_spec
#from ..filter.filter import parse_filter
#from ..basics.errorbar import ErrorBar
#from ..basics.unit import PhotometricUnit, parse_unit
#from ..basics.quantity import parse_quantity

# -----------------------------------------------------------------

def integer_or_real_or_string(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    # First try integer, then real, then string
    try: return integer(argument)
    except ValueError:
        try: return real(argument)
        except ValueError:
            return string(argument)

# -----------------------------------------------------------------

def integer_or_real_or_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try: return integer(argument)
    except ValueError:
        try: return real(argument)
        except ValueError:
            return quantity(argument)

# -----------------------------------------------------------------

def real_or_quantity(argument):

    """
    This fucntion ...
    :param argument:
    :return:
    """

    try: return real(argument)
    except ValueError: return quantity(argument)

# -----------------------------------------------------------------

def real_list_or_quantity_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try: return real_list(argument)
    except ValueError: return quantity_list(argument)

# -----------------------------------------------------------------

def integer_real_and_quantity_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    elements = string_list(argument)

    values = []

    for element in elements:
        value = integer_or_real_or_quantity(element)
        values.append(value)

    return values

# -----------------------------------------------------------------

def bit(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    value = integer(argument)
    if value == 0 or value == 1: return value
    else: raise ValueError("")

# -----------------------------------------------------------------

def digit(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    value = integer(argument)
    if 0 <= value <= 9: return value
    else: raise ValueError("Not a digit: must be integer between 0 and 9")

# -----------------------------------------------------------------

def binary(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    chars = characters(argument)
    return [bit(char) for char in chars]

# -----------------------------------------------------------------

def characters(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    return list(argument)

# -----------------------------------------------------------------

def character(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    chars = character(argument)
    if len(chars) == 1: return chars[0]
    else: raise ValueError("Contains multiple characters")

# -----------------------------------------------------------------

def letter(argument):

    """
    THis function ...
    :param argument: 
    :return: 
    """

    char = character(argument)
    if char.isalpha(): return char
    else: raise ValueError("Is not a letter character")

# -----------------------------------------------------------------

def integer_array(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    import numpy as np
    return np.array(integer_list(argument))

# -----------------------------------------------------------------

def real_array(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    import numpy as np
    return np.array(real_list(argument))

# -----------------------------------------------------------------

def boolean(entry):

    """
    Boolean value (True or False). Allowed: 'True', 'T', 'y', 'yes', 'False', 'n', 'no'
    :param entry:
    :return:
    """

    lowercase = entry.lower().strip()

    if lowercase == "true" or lowercase == "y" or lowercase == "yes" or lowercase == "t": return True
    elif lowercase == "false" or lowercase == "n" or lowercase == "no" or lowercase == "f": return False
    else: raise ValueError("Invalid boolean specification: " + entry)

# -----------------------------------------------------------------

def decimal(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return Decimal(argument)

# -----------------------------------------------------------------

def positive_decimal(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = decimal(argument)
    if value < 0: raise ValueError("Value is smaller than zero")
    return value

# -----------------------------------------------------------------

def negative_decimal(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = decimal(argument)
    if value > 0: raise ValueError("Value is greater than zero")
    return value

# -----------------------------------------------------------------

def integer(argument):

    """
    Integer value
    :param argument:
    :return:
    """

    try: return int(argument)
    except ValueError: # ValueError: invalid literal for int() with base 10:

        from . import numbers

        # Parse as decimal number (to keep precision)
        value = decimal(argument)

        # Check whether integer
        if not numbers.is_integer(value): raise ValueError("Not an integer number")

        # Convert to int and return
        return int(value)

# -----------------------------------------------------------------

def positive_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = integer(argument)
    if value < 0: raise ValueError("Value is smaller than zero")
    return value

# -----------------------------------------------------------------

def negative_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = integer(argument)
    if value > 0: raise ValueError("Value is greater than zero")
    return value

# -----------------------------------------------------------------

def even_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = integer(argument)
    if value % 2 != 0: raise ValueError("Integer is not even")
    return value

# -----------------------------------------------------------------

def even_positive_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = positive_integer(argument)
    if value % 2 != 0: raise ValueError("Integer is not even")

# -----------------------------------------------------------------

def even_negative_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = negative_integer(argument)
    if value % 2 != 0: raise ValueError("Integer is not even")
    return value

# -----------------------------------------------------------------

def odd_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = integer(argument)
    if value % 2 == 0: raise ValueError("Integer is not odd")
    return value

# -----------------------------------------------------------------

def odd_positive_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = positive_integer(argument)
    if value % 2 == 0: raise ValueError("Integer is not odd")
    return value

# -----------------------------------------------------------------

def odd_negative_integer(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = negative_integer(argument)
    if value % 2 == 0: raise ValueError("Integer is not odd")
    return value

# -----------------------------------------------------------------

def real(argument):

    """
    Real (floating-point) value
    :param argument:
    :return:
    """

    return float(argument)

# -----------------------------------------------------------------

def positive_real(argument):

    """
    Positive real (floating-point) value (>=0)
    :param argument:
    :return:
    """

    value = real(argument)
    if value < 0: raise ValueError("Value is smaller than zero")
    return value

# -----------------------------------------------------------------

def negative_real(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = real(argument)
    if value > 0: raise ValueError("Value is greater than zero")
    return value

# -----------------------------------------------------------------

def fraction(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    value = real(argument)
    if value > 1 or value < 0: raise ValueError("Value should be from 0 to 1")
    return value

# -----------------------------------------------------------------

def string(argument):

    """
    String
    :param argument:
    :return:
    """

    return argument

# -----------------------------------------------------------------

def nonempty_string(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    if argument == "": raise ValueError("String is empty")
    return argument

# -----------------------------------------------------------------

def string_no_spaces(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    if " " in argument: raise ValueError("The string cannot contain spaces")
    return argument

# -----------------------------------------------------------------

def nonempty_string_no_spaces(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    argument = string_no_spaces(argument)
    return nonempty_string(argument)

# -----------------------------------------------------------------

def real_range(argument):

    """
    Range of real (floating-point) values
    :param argument:
    :return:
    """

    from ..basics.range import IntegerRange, RealRange, QuantityRange

    min_value, max_value = real_tuple(argument.replace(">", ","))
    return RealRange(min_value, max_value)

# -----------------------------------------------------------------

def integer_range(argument):

    """
    Range of integer values
    :param argument:
    :return:
    """

    from ..basics.range import IntegerRange, RealRange, QuantityRange

    min_value, max_value = integer_tuple(argument.replace(">", ","))
    return IntegerRange(min_value, max_value)

# -----------------------------------------------------------------

def quantity_range(argument):

    """
    Range of (Astropy) quantities
    :param argument:
    :return:
    """

    from ..basics.range import QuantityRange
    min_quantity, max_quantity = quantity_tuple(argument.replace(">", ","))
    return QuantityRange(min_quantity, max_quantity)

# -----------------------------------------------------------------

def photometric_quantity_range(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..basics.range import QuantityRange
    min_quantity, max_quantity = photometric_quantity_tuple(argument.replace(">", ","))
    return QuantityRange(min_quantity, max_quantity)

# -----------------------------------------------------------------

def photometric_density_quantity_range(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..basics.range import QuantityRange
    min_quantity, max_quantity = photometric_density_quantity_tuple(argument.replace(">", ","))
    return QuantityRange(min_quantity, max_quantity)

# -----------------------------------------------------------------

def range(argument):

    """
    Integer, real or quantity range
    :param argument:
    :return:
    """

    try: range = integer_range(argument)
    except ValueError:
        try: range = real_range(argument)
        except ValueError:
            try: range = quantity_range(argument)
            except ValueError: raise ValueError("Not a valid range")

    return range

# -----------------------------------------------------------------

def directory_path(argument):

    """
    Converts a relative path or directory name to an absolute directory path, and checks whether this
    directory exists
    :param argument:
    :return:
    """

    path = fs.absolute_path(argument)
    if not fs.is_directory(path): raise ValueError("Is not a directory: " + path)
    return path

# -----------------------------------------------------------------

def directorypath_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return [directory_path(path) for path in string_list(argument)]

# -----------------------------------------------------------------

def path(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    path = fs.absolute_path(argument)
    if not fs.is_file(path) and not fs.is_directory(path): raise ValueError("Is not a file or directory: " + path)
    return path

# -----------------------------------------------------------------

def path_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return [path(p) for p in string_list(argument)]

# -----------------------------------------------------------------

def file_path(argument):

    """
    Converts a relative path or filename to an absolute filepath, and checks whether this file exists
    :param argument:
    :return:
    """

    path = fs.absolute_path(argument)
    if not fs.is_file(path): raise ValueError("Is not a file: " + path)
    return path

# -----------------------------------------------------------------

def filepath_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return [file_path(path) for path in string_list(argument)]

# -----------------------------------------------------------------

def string_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = argument.split(",")
        return a, b
    except: raise ValueError("Tuple must be of format a,b")

# -----------------------------------------------------------------

def integer_tuple(argument):

    """
    Tuple of integer values
    :param argument:
    :return:
    """

    try:
        a, b = map(int, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a,b")

# -----------------------------------------------------------------

def real_tuple(argument):

    """
    Tuple of real (floating-point) values
    :param argument:
    :return:
    """

    try:
        a, b = map(float, argument.split(","))
        return a, b
    except: raise ValueError("Tuple must be of format a,b")

# -----------------------------------------------------------------

def quantity_tuple(argument):

    """
    Tuple of (Astropy) quantities
    :param argument:
    :return:
    """

    try:
        a, b = map(quantity, argument.split(","))
        return a, b
    except: raise ValueError("Tuple must be of format a unit_a, b unit_b")

# -----------------------------------------------------------------

def angle_tuple(argument):

    """
    Tuple of angles
    :param argument:
    :return:
    """

    try:
        a, b = map(angle, argument.split(","))
        return a, b
    except: raise ValueError("Tuple must be of format a unit_a, b unit_b")

# -----------------------------------------------------------------

def quantity_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    quantities = []
    for item in string_list(argument): quantities.append(quantity(item))
    return quantities

# -----------------------------------------------------------------

def ascending_quantity_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from . import sequences
    quantities = quantity_list(argument)
    if not sequences.is_ascending(quantities): raise ValueError("List is not ascending")
    return quantities

# -----------------------------------------------------------------

def descending_quantity_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from . import sequences
    quantities = quantity_list(argument)
    if not sequences.is_descending(quantities): raise ValueError("List is not descending")
    return quantities

# -----------------------------------------------------------------

def mixed_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return tuple(argument.split(","))

# -----------------------------------------------------------------

def quantity_vector(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    tuple_string = argument.split("(")[1].split(")")[0]
    x, y = tuple_string.split(", ")

    x = quantity(x[2:])
    y = quantity(y[2:])

    from ...magic.basics.vector import Vector
    return Vector(x, y)

# -----------------------------------------------------------------

def string_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    if argument == "": return []
    else:
        parts = argument.split(",")
        return [part.strip() for part in parts]

# -----------------------------------------------------------------

def string_or_string_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try: return string_list(argument)
    except ValueError: return string(argument)

# -----------------------------------------------------------------

def mixed_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    if argument == "": return []
    else:
        #return [eval(value) for value in argument.split(",")]
        values = []
        for value in argument.split(","):
            try: value = eval(value)
            except: warnings.warn("Could not evaluate all elements of the mixed list '" + str(argument) + "', so treating as strings")
            values.append(value)
        return values

# -----------------------------------------------------------------

def duration(argument):

    """
    Duration in seconds from hh:mm:ss format
    :param argument:
    :return:
    """

    # Calculate the walltime in seconds
    hours, minutes, seconds = argument.split(':')
    duration = int(hours)*3600 + int(minutes)*60 + int(seconds)

    # Return the duration in seconds
    return duration

# -----------------------------------------------------------------

def integer_list(argument):

    """
    A list of integer values, based on a string denoting a certain range (e.g. '3-9') or a
    set of integer values seperated by commas ('2,14,20')
    :param argument:
    :return:
    """

    if argument == "": return []

    if "-" in argument and "," in argument:

        parts = argument.split(",")
        #print("PARTS", parts)
        total_int_list = []
        for part in parts: total_int_list += integer_list(part)
        return total_int_list

    # Split the string
    splitted = argument.split('-')

    #print("SPLITTED", splitted)

    if len(splitted) == 0: raise ValueError("No range given")
    elif len(splitted) == 1:

        splitted = splitted[0].split(",")

        #print("SPLITTED 2", splitted)

        # Check if the values are valid
        for value in splitted:
            if not value.isdigit(): raise ValueError("Argument contains unvalid characters")

        # Only leave unique values: NO!!!!!!??? WHY WAS THIS EVER HERE???
        #return list(set([int(value) for value in splitted]))
        return [integer(value) for value in splitted]

    elif len(splitted) == 2:

        if not (splitted[0].isdigit() and splitted[1].isdigit()): ValueError("Not a valid integer range")
        return range(int(splitted[0]), int(splitted[1])+1)

    else: raise ValueError("Values must be seperated by commas or by a '-' in the case of a range")

# -----------------------------------------------------------------

def ascending_integer_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from . import sequences
    integers = integer_list(argument)
    if not sequences.is_ascending(integers): raise ValueError("List is not ascending")
    return integers

# -----------------------------------------------------------------

def descending_integer_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from . import sequences
    integers = integer_list(argument)
    if not sequences.is_descending(integers): raise ValueError("List is not descending")
    return integers

# -----------------------------------------------------------------

def real_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return [real(value) for value in string_list(argument)]

# -----------------------------------------------------------------

def ascending_real_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from . import sequences
    reals = real_list(argument)
    #print(reals)
    if not sequences.is_ascending(reals): raise ValueError("List is not ascending")
    return reals

# -----------------------------------------------------------------

def descending_real_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from . import sequences
    reals = real_list(argument)
    if not sequences.is_descending(reals): raise ValueError("List is not descending")
    return reals

# -----------------------------------------------------------------

def ordered_dictionary(argument):

    """
    Thisf unction ...
    :param argument:
    :return:
    """

    if argument == "": return OrderedDict()
    else:
        d = eval("OrderedDict([(" + argument.replace(",", "), (").replace(":", ", ") + ")])")
        return d

# -----------------------------------------------------------------

def dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    if argument == "": return dict()
    else:
        d = eval("{" + argument + "}")
        if not isinstance(d, dict): raise ValueError("Not a proper specification of a dictionary")
        return d

# -----------------------------------------------------------------

def string_string_dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        if not types.is_string_type(value): raise ValueError("All values must be strings")
    return d

# -----------------------------------------------------------------

def string_string_list_dictionary(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        if not types.is_string_sequence(value): raise ValueError("All values must be string sequences")
    return d

# -----------------------------------------------------------------

def string_integer_dictionary(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        if not types.is_integer_type(value): raise ValueError("All values must be integers")
    return d

# -----------------------------------------------------------------

def string_real_dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        if not types.is_real_type(value): raise ValueError("All values must be real")
    return d

# -----------------------------------------------------------------

def filter_real_dictionary(argument):

    """
    Thisf unction ...
    :param argument:
    :return:
    """

    d = dictionary(argument)
    new = dict()
    for key, value in d.items():
        fltr = filter(key)
        new[fltr] = value

    # Return the dictionary
    return new

# -----------------------------------------------------------------

def string_filepath_dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        d[key] = file_path(value) # check if parsing as filepath succeeds
    return d

# -----------------------------------------------------------------

def filepath_or_string_filepath_dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try: return string_filepath_dictionary(argument)
    except ValueError: return file_path(argument)

# -----------------------------------------------------------------

def string_unit_dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        d[key] = unit(value)
    return d

# -----------------------------------------------------------------

def string_photometricunit_dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        d[key] = photometric_unit(value)
    return d

# -----------------------------------------------------------------

def string_tuple_dictionary(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    d = dictionary(argument)
    for key, value in d.items():
        if not types.is_string_type(key): raise ValueError("All keys must be strings")
        if not types.is_tuple(value): raise ValueError("All values must be tuples")
    return d

# -----------------------------------------------------------------

def unit(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_unit
    return parse_unit(argument) # can be photometric, but doesn't need to be

# -----------------------------------------------------------------

def photometric_unit(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.unit import PhotometricUnit
    return PhotometricUnit(argument)

# -----------------------------------------------------------------

def photometric_density_unit(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.unit import PhotometricUnit
    return PhotometricUnit(argument, density=True, density_strict=True)

# -----------------------------------------------------------------

def photometric_brightness_unit(argument):

    """
    THis function ...
    :param argument: 
    :return: 
    """

    from ..units.unit import PhotometricUnit
    return PhotometricUnit(argument, brightness=True, brightness_strict=True)

# -----------------------------------------------------------------

def length_unit(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_unit
    unit = parse_unit(argument)
    if unit.physical_type != "length": raise ValueError("Not a length unit")
    else: return unit

# -----------------------------------------------------------------

def quantity(argument):

    """
    An Astropy quantity
    """

    from ..units.parsing import parse_quantity
    return parse_quantity(argument)

# -----------------------------------------------------------------

def photometric_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    return parse_quantity(argument)

# -----------------------------------------------------------------

def photometric_quantity_tuple(argument):

    """
    This fucntion ...
    :param argument:
    :return:
    """

    try:
        a, b = map(photometric_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a photometric_unit_a, b photometric_unit_b")

# -----------------------------------------------------------------

def photometric_density_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    return parse_quantity(argument, density=True)

# -----------------------------------------------------------------

def photometric_density_quantity_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(photometric_density_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a photometric_density_unit_a, b photometric_density_unit_b")

# -----------------------------------------------------------------

def length_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    qty = parse_quantity(argument)
    if qty.unit.physical_type != "length": raise ValueError("Not a length")
    return qty

# -----------------------------------------------------------------

def length_quantity_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(length_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a length_unit_a, b length_unit_b")

# -----------------------------------------------------------------

def temperature_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    qty = parse_quantity(argument)
    if qty.unit.physical_type != "temperature": raise ValueError("Not a temperature")
    return qty

# -----------------------------------------------------------------

def temperature_quantity_tuple(argument):

    """
    This fucntion ...
    :param argument:
    :return:
    """

    try:
        a, b = map(temperature_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a temperature_unit_a, temperature_unit_b")

# -----------------------------------------------------------------

def mass_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    qty = parse_quantity(argument)
    if qty.unit.physical_type != "mass": raise ValueError("Not a mass")
    return qty

# -----------------------------------------------------------------

def mass_quantity_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(mass_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a mass_unit_a, mass_unit_b")

# -----------------------------------------------------------------

def mass_density_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    qty = parse_quantity(argument)
    if qty.unit.physical_type != "mass density": raise ValueError("Not a mass density")
    return qty

# -----------------------------------------------------------------

def mass_density_quantity_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(mass_density_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a mass_density_unit_a, b mass_density_unit_b")

# -----------------------------------------------------------------

def mass_surface_density_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    qty = parse_quantity(argument)

    # Check
    b = qty.unit / "m"
    if b.physical_type != "mass density": raise ValueError("Not a mass surface density")

    # Return
    return qty

# -----------------------------------------------------------------

def mass_surface_density_quantity_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(mass_surface_density_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a mass_surface_density_unit_a, b mass_surface_density_unit_b")

# -----------------------------------------------------------------

def data_quantity(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..units.parsing import parse_quantity
    qty = parse_quantity(argument)
    if qty.unit.physical_type != "data quantity": raise ValueError("Not a data quantity")
    return qty

# -----------------------------------------------------------------

def data_quantity_tuple(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    try:
        a, b = map(data_quantity, argument.split(','))
        return a, b
    except: raise ValueError("Tuple must be of format a data_unit_a, b data_unit_b")

# -----------------------------------------------------------------

def angle(argument):

    """
    An Astropy Angle
    :param argument:
    :return:
    """

    from ..units.parsing import parse_angle
    return parse_angle(argument)

# -----------------------------------------------------------------

def errorbar(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..basics.errorbar import ErrorBar

    upper = None
    if ">" in argument: lower, upper = quantity_tuple(argument.replace(">", ","))
    else: lower = quantity(argument)

    # Create error bar
    return ErrorBar(lower, upper)

# -----------------------------------------------------------------

def photometric_errorbar(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..basics.errorbar import ErrorBar

    upper = None
    if ">" in argument: lower, upper = photometric_quantity_tuple(argument)
    else: lower = photometric_quantity(argument)

    # Create error bar
    return ErrorBar(lower, upper)

# -----------------------------------------------------------------

def pixel_limits(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    lst = integer_list(argument)
    assert len(lst) == 4
    return lst

# -----------------------------------------------------------------

def calibration_error(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ...magic.misc.calibration import CalibrationError
    return CalibrationError.from_string(argument)

# -----------------------------------------------------------------

def url(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    regex = re.compile(
        r'^(?:http|ftp)s?://'  # http:// or https://
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain...
        r'localhost|'  # localhost...
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # ...or ip
        r'(?::\d+)?'  # optional port
        r'(?:/?|[/?]\S+)$', re.IGNORECASE)

    if not regex.match(argument): raise ValueError("Invalid URL")
    else: return argument

# -----------------------------------------------------------------

def url_list(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    return [url(element) for element in string_list(argument)]

# -----------------------------------------------------------------

def image_path(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    path = file_path(argument)
    if path.endswith("fits"): raise ValueError("Unrecognized file type")
    else: return path

# -----------------------------------------------------------------

def filter(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..filter.filter import parse_filter
    return parse_filter(argument)

# -----------------------------------------------------------------

def filter_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return [filter(arg) for arg in string_list(argument)]

# -----------------------------------------------------------------

def broad_band_filter(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..filter.broad import BroadBandFilter
    return BroadBandFilter(argument)

# -----------------------------------------------------------------

def narrow_band_filter(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..filter.narrow import NarrowBandFilter
    return NarrowBandFilter(argument)

# -----------------------------------------------------------------

def broad_band_filter_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    filters = []
    for item in string_list(argument): filters.append(broad_band_filter(item))
    return filters

# -----------------------------------------------------------------

def narrow_band_filter_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    filters = []
    for item in string_list(argument): filters.append(narrow_band_filter(item))
    return filters

# -----------------------------------------------------------------

def lazy_filter_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..filter.broad import BroadBandFilter
    from ..filter.broad import identifiers as broad_band_identifiers
    from ..filter.narrow import NarrowBandFilter, generate_aliases_ranges, wavelength_range_for_spec
    from ..filter.filter import parse_filter

    filters = []
    for arg in string_list(argument):

        try:

            # Try to parse the filter
            fltr = parse_filter(arg)
            filters.append(fltr)

        # If parsing directly failes
        except ValueError:

            # Try matching with broad bands
            for spec in broad_band_identifiers:

                identifier = broad_band_identifiers[spec]

                if "instruments" in identifier:
                    if arg in identifier.instruments:
                        filters.append(BroadBandFilter(spec))
                        continue # this filter matches
                if "observatories" in identifier:
                    if arg in identifier.observatories:
                        filters.append(BroadBandFilter(spec))
                        continue # this filter matches

            # Try matching with narrow bands defined by wavelength ranges
            for spec, alias in generate_aliases_ranges():

                if alias not in argument: continue

                # Get wavelength range
                wavelength_range = wavelength_range_for_spec(spec)

                # Create two filters, one for the minimum and one for the maximum wavelength
                fltr_min = NarrowBandFilter(wavelength_range.min, name=alias + " min")
                fltr_max = NarrowBandFilter(wavelength_range.max, name=alias + " max")

                filters.append(fltr_min)
                filters.append(fltr_max)

    # Return the list of filters
    return filters

# -----------------------------------------------------------------

def lazy_broad_band_filter_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..filter.broad import BroadBandFilter
    from ..filter.broad import identifiers as broad_band_identifiers

    filters = []
    for arg in string_list(argument):

        try:
            # Try to parse the filter
            fltr = BroadBandFilter(arg)
            filters.append(fltr)

        except ValueError:

            # Try matching with broad bands
            for spec in broad_band_identifiers:

                identifier = broad_band_identifiers[spec]

                #print(spec)

                if "instruments" in identifier:
                    if arg in identifier.instruments:
                        filters.append(BroadBandFilter(spec))
                        continue  # this filter matches
                if "observatories" in identifier:
                    if arg in identifier.observatories:
                        filters.append(BroadBandFilter(spec))
                        continue  # this filter matches

    # Return the filters
    return filters

# -----------------------------------------------------------------

def lazy_narrow_band_filter_list(argument):

    """
    This fucntion ...
    :param argument:
    :return:
    """

    from ..filter.narrow import NarrowBandFilter, generate_aliases_ranges, wavelength_range_for_spec

    filters = []
    for arg in string_list(argument):

        try:
            # Try to parse the filter
            fltr = NarrowBandFilter(arg)
            filters.append(fltr)

        except ValueError:

            # Try matching with narrow bands defined by wavelength ranges
            for spec, alias in generate_aliases_ranges():

                if alias not in argument: continue

                # Get wavelength range
                wavelength_range = wavelength_range_for_spec(spec)

                # Create two filters, one for the minimum and one for the maximum wavelength
                fltr_min = NarrowBandFilter(wavelength_range.min, name=alias + " min")
                fltr_max = NarrowBandFilter(wavelength_range.max, name=alias + " max")

                filters.append(fltr_min)
                filters.append(fltr_max)

    # Return the filters
    return filters

# -----------------------------------------------------------------

def pixelcoordinate(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ...magic.basics.coordinate import PixelCoordinate
    x, y = real_tuple(argument)
    return PixelCoordinate(x, y)

# -----------------------------------------------------------------

def skycoordinate(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ...magic.basics.coordinate import SkyCoordinate
    ra, dec = quantity_tuple(argument)
    return SkyCoordinate(ra=ra, dec=dec)

# -----------------------------------------------------------------

def sed_entry(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    fltr, flux, flux_error = argument.split("::")

    fltr = filter(fltr)
    flux = photometric_quantity(flux)
    flux_error = photometric_errorbar(flux_error)

    return fltr, flux, flux_error

# -----------------------------------------------------------------

def sed_entry_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    entries = []
    for item in string_list(argument): entries.append(sed_entry(item))
    return entries

# -----------------------------------------------------------------

def instrument_frame(argument):

    """
    This function ....
    :param argument:
    :return:
    """

    from ...modeling.basics.instruments import InstrumentFrame
    return InstrumentFrame(**dictionary(argument))

# -----------------------------------------------------------------

def instrument_frame_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return [instrument_frame(arg) for arg in string_list(argument)]

# -----------------------------------------------------------------

def parts(argument):

    """
    This function ...
    :param argument: 
    :return: 
    """

    import numpy as np
    values = real_array(argument)
    relative_values = list(values / np.sum(values))
    return relative_values

# -----------------------------------------------------------------

def weights(argument):

    """
    This fucntion ...
    :param argument: 
    :return: 
    """

    import numpy as np
    ratios = np.array(parts(argument))
    the_weights = list(ratios * len(ratios))
    return the_weights

# -----------------------------------------------------------------

def pixel_shape(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ...magic.basics.vector import PixelShape
    shape = integer_tuple(argument)
    return PixelShape.from_xy_tuple(shape)

# -----------------------------------------------------------------

def colour(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    from ..basics.colour import parse_colour
    return parse_colour(argument)

# -----------------------------------------------------------------

def percentage(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return real(argument) / 100.

# -----------------------------------------------------------------

def integer_extent(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    x, y = integer_tuple(argument)
    from ...magic.basics.vector import IntegerExtent
    return IntegerExtent(x, y)

# -----------------------------------------------------------------

def real_extent(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    x, y = real_tuple(argument)
    from ...magic.basics.vector import RealExtent
    return RealExtent(x, y)

# -----------------------------------------------------------------

def angle_extent(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    x, y = angle_tuple(argument)
    from ...magic.basics.vector import AngleExtent
    return AngleExtent(x, y)

# -----------------------------------------------------------------

def quantity_extent(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    x, y = quantity_tuple(argument)
    from ...magic.basics.vector import QuantityExtent
    return QuantityExtent(x, y)

# -----------------------------------------------------------------

def sky_extent(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    x, y = angle_tuple(argument)
    from ...magic.basics.stretch import SkyExtent
    return SkyExtent(x, y)

# -----------------------------------------------------------------

def physical_extent(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    x, y = length_quantity_tuple(argument)
    from ...magic.basics.stretch import PhysicalExtent
    return PhysicalExtent(x, y)

# -----------------------------------------------------------------

def parallelization(argument, default_nthreads_per_core=1):

    """
    This function ...
    :param argument:
    :param default_nthreads_per_core:
    :return:
    """

    from ..simulation.parallelization import Parallelization

    # Parse
    splitted = argument.split(":")
    if len(splitted) < 2: raise ValueError("Invalid input: must be 'ncores:nprocesses[:nthreads_per_core]'")
    if len(splitted) > 3: raise ValueError("Invalid input: must be 'ncores:nprocesses[:nthreads_per_core]'")

    ncores = integer(splitted[0])
    nprocesses = integer(splitted[1])

    # Get number of threads per core (hyperthreading)
    if len(splitted) == 3: nthreads_per_core = integer(splitted[2])
    else: nthreads_per_core = default_nthreads_per_core

    # Create and return the parallelization
    return Parallelization(ncores, nthreads_per_core, nprocesses)

# -----------------------------------------------------------------
