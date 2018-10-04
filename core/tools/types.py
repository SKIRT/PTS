#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.types Checking types of objects.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from types import NoneType
import collections
try:
    HAS_NP = True
    import numpy as np
except ImportError: HAS_NP = False

# -----------------------------------------------------------------

# IMPLEMENT INTELLIGENT META-PROGRAMMING TYPING THINGS:
# SEE: https://www.youtube.com/watch?v=sPiWg5jSoZI (around 1:10 ->  1:30 etc.) !!

# -----------------------------------------------------------------

if HAS_NP: boolean_types = [bool, np.bool, np.bool_]
else: boolean_types = [bool]

# -----------------------------------------------------------------

def is_boolean_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in boolean_types:
        if type(value) == tp: return True
    return False

# -----------------------------------------------------------------

if HAS_NP: integer_types = [int, long, np.int32, np.int64, np.uint32, np.uint64]
else: integer_types = [int, long]

# -----------------------------------------------------------------

def is_integer_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in integer_types:
        if type(value) == tp: return True
    return False

# -----------------------------------------------------------------

if HAS_NP: real_types = [float, np.float32, np.float64]
else: real_types = [float]

# -----------------------------------------------------------------

def is_real_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in real_types:
        if type(value) == tp: return True
    return False

# -----------------------------------------------------------------

def is_real_or_integer(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_real_type(value) or is_integer_type(value)

# -----------------------------------------------------------------

if HAS_NP: string_types = [basestring, str, np.string_]
else: string_types = [basestring, str]

# -----------------------------------------------------------------

def is_string_type(value):

    """
    This function ...
    :param value:
    :return:
    """

    for tp in string_types:

        # Literal
        if type(value) == tp: return True

        # Derived from
        if isinstance(value, tp): return True

    return False

# -----------------------------------------------------------------

def is_dictionary(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    return isinstance(value, dict)

# -----------------------------------------------------------------

def is_dictionary_of_dictionaries(value, passive=False):

    """
    This function ...
    :param value:
    :param passive:
    :return: 
    """

    if not is_dictionary(value): raise ValueError("Not a dictionary")
    if len(value) == 0: raise ValueError("Empty dictionary")

    if passive:

        for key in value:
            # from the moment one element is also a dict, return True
            if is_dictionary(value[key]): return True
            else: return False
        return False

    else:

        # Loop over the keys
        for key in value:
            if not is_dictionary(value[key]): return False
        return True

# -----------------------------------------------------------------

def is_sequence(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    #print(type(value))

    #return isinstance(value, collections.Iterable) and not is_string_type(value) and not is_dictionary(value)
    return isinstance(value, collections.Sequence) and not is_string_type(value) and not is_tuple(value) # NEW

# -----------------------------------------------------------------

def is_list(value):

    """
    This function ...
    :param value:
    :return:
    """

    return isinstance(value, list)

# -----------------------------------------------------------------

def is_string_sequence(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    if not is_sequence(value): return False
    for element in value:
        if not is_string_type(element): return False
    return True

# -----------------------------------------------------------------

def is_integer_sequence(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_sequence(value): return False
    for element in value:
        if not is_integer_type(element): return False
    return True

# -----------------------------------------------------------------

def is_real_sequence(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_sequence(value): return False
    for element in value:
        if not is_real_type(element): return False
    return True

# -----------------------------------------------------------------

def get_dtype_name(dtype):

    """
    This function ...
    :param name:
    :return:
    """

    if dtype.name.startswith("string"): return "string"
    elif dtype.name.startswith("float"): return "real"
    elif dtype.name.startswith("int"): return "integer"
    elif dtype.name.startswith("bool"): return "boolean"
    else: raise ValueError("Unknown type: '" + dtype.name + "'")

# -----------------------------------------------------------------

def is_string_array(value):

    """
    Thisn function ...
    :param value:
    :return:
    """

    if not is_array_like(value): return False

    # WRONG! loops over rows if 2D!
    #for element in value:
    #    if not is_string_type(element): return False
    #return True

    return get_dtype_name(value.dtype) == "string"

# -----------------------------------------------------------------

def is_integer_array(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_array_like(value): return False

    # WRONG! loops over rows if 2D!
    #for element in value:
    #    if not is_integer_type(element): return False
    #return True

    return get_dtype_name(value.dtype) == "integer"

# -----------------------------------------------------------------

def is_real_array(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_array_like(value): return False

    # WRONG! loops over rows if 2D!
    #for element in value:
    #    if not is_real_type(element): return False
    #return True

    return get_dtype_name(value.dtype) == "real"

# -----------------------------------------------------------------

def is_real_or_integer_array(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_real_array(value) or is_integer_array(value)

# -----------------------------------------------------------------

def is_boolean_array(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_array_like(value): return False

    return get_dtype_name(value.dtype) == "boolean"

# -----------------------------------------------------------------

def is_string_column(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_astropy_column(value): return False
    for element in value:
        if not is_string_type(element): return False
    return True

# -----------------------------------------------------------------

def is_integer_column(value):

    """
    Thisn function ...
    :param value:
    :return:
    """

    if not is_astropy_column(value): return False
    for element in value:
        if not is_integer_type(element): return False
    return True

# -----------------------------------------------------------------

def is_real_column(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_astropy_column(value): return False
    for element in value:
        if not is_real_type(element): return False
    return True

# -----------------------------------------------------------------

def is_tuple(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    return isinstance(value, tuple)

# -----------------------------------------------------------------

def is_string_tuple(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_tuple(value): return False
    for element in value:
        if not is_string_type(element): return False
    return True

# -----------------------------------------------------------------

def is_integer_tuple(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not is_tuple(value): return False
    for element in value:
        if not is_integer_type(element): return False
    return True

# -----------------------------------------------------------------

def is_none(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    return isinstance(value, NoneType)

# -----------------------------------------------------------------

def is_unit(value):

    """
    This fucntion ...
    :param value:
    :return:
    """

    #from astropy.units import Unit
    from astropy.units import UnitBase # also IrreducibleUnit etc.
    return isinstance(value, UnitBase)

# -----------------------------------------------------------------

def is_quantity(value):

    """
    This function ...
    :param value:
    :return:
    """

    from astropy.units import Quantity
    return isinstance(value, Quantity)

# -----------------------------------------------------------------

def is_angle(value):

    """
    This function ...
    :param value:
    :return:
    """

    return hasattr(value, "unit") and value.unit.physical_type == "angle"

# -----------------------------------------------------------------

def is_length_unit(value):

    """
    Thisn function ...
    :param value:
    :return:
    """

    return is_unit(value) and value.physical_type == "length"

# -----------------------------------------------------------------

def is_length_quantity(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_quantity(value) and value.unit.physical_type == "length"

# -----------------------------------------------------------------

def is_array_like(value):

    """
    This function ...
    :param value:
    :return:
    """

    return hasattr(value, "__array__") and hasattr(value, "__getitem__") and can_get_item(value) and not hasattr(value, "unit") and not is_boolean_type(value) and not is_integer_type(value) and not is_real_type(value) and not is_string_type(value)

# -----------------------------------------------------------------

def is_astropy_column(value):

    """
    This function ...
    :param value:
    :return:
    """

    return type(value).__name__ == "MaskedColumn" or type(value).__name__ == "Column"

# -----------------------------------------------------------------

def can_get_item(value):

    """
    This function ...
    :param value:
    :return:
    """

    #print(value, type(value))

    try:
        length = len(value)
    except TypeError: return False

    if len(value) == 0: return True
    else:
        try:
            item = value[0]
            return True
        except IndexError: return False

# -----------------------------------------------------------------

def is_sequence_or_array(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_sequence(value) or is_array_like(value)

# -----------------------------------------------------------------

def is_sequence_or_tuple(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_sequence(value) or is_tuple(value)

# -----------------------------------------------------------------

def is_integer_sequence_or_tuple(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_integer_sequence(value) or is_integer_tuple(value)

# -----------------------------------------------------------------

def is_frame(value):

    """
    This function ...
    :param value:
    :return:
    """

    from ...magic.core.frame import Frame
    return isinstance(value, Frame)

# -----------------------------------------------------------------

def is_image(value):

    """
    This function ...
    :param value:
    :return:
    """

    from ...magic.core.image import Image
    return isinstance(value, Image)

# -----------------------------------------------------------------

def is_datacube(value):

    """
    This function ...
    :param value:
    :return:
    """

    from ...magic.core.datacube import DataCube
    return isinstance(value, DataCube)

# -----------------------------------------------------------------

def is_2d_array(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_array_like(value) and value.ndim == 2

# -----------------------------------------------------------------

def is_3d_array(value):

    """
    This function ...
    :param value:
    :return:
    """

    return is_array_like(value) and value.ndim == 3

# -----------------------------------------------------------------`
