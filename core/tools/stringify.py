#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.stringify Provides useful functions for converting objects of various types to strings.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy

# Import the relevant PTS classes and modules
from . import types
from . import introspection
from . import sequences
from . import strings
from . import numbers

# -----------------------------------------------------------------

#def tostr(value, scientific=None, decimal_places=2, fancy=False, ndigits=None, delimiter=", ", **kwargs):
def tostr(value, **kwargs):

    """
    This function ...
    :param value:
    :param kwargs:
    :return: 
    """

    # Get the 'scientific' flag
    scientific = kwargs.get("scientific", None)

    # Set default number of decimal places
    decimal_places = 2

    # Set scientific flag flexibly, if scientific flag was not passed explicitly
    if scientific is None:

        # Integer value
        if types.is_integer_type(value) or (types.is_real_type(value) and numbers.is_integer(value)):

            # Convert to be certain (if from float)
            value = int(value)

            #if -1e4 <= value <= 1e4: scientific = False
            if -999 < value < 999: scientific = False
            else: scientific = True

            # No decimals for integers
            decimal_places = 0

        # Real value
        elif types.is_real_type(value):

            #if -1e4 <= value <= 1e4: scientific = False
            if -999.99 < value < 999.99: scientific = False
            else: scientific = True

        # Quantity
        elif introspection.lazy_isinstance(value, "Quantity", "astropy.units", return_false_if_fail=True):

            if -999.99 < value.value < 999.99: scientific = False
            else: scientific = True

        elif introspection.lazy_isinstance(value, "QuantityRange", "pts.core.basics.range", return_false_if_fail=True):

            if -999.99 < value.min.value and value.max.value < 999.99: scientific = False
            else: scientific = True

        elif introspection.lazy_isinstance(value, "RealRange", "pts.core.basics.range", return_false_if_fail=True):

            if -999.99 < value.min and value.max < 999.99: scientific = False
            else: scientific = True

        elif introspection.lazy_isinstance(value, "IntegerRange", "pts.core.basics.range", return_false_if_fail=True):

            if -999 < value.min and value.max < 999: scientific = False
            else: scientific = True

        # Other
        else: scientific = False

        # Set the options
        kwargs["scientific"] = scientific
        kwargs["decimal_places"] = decimal_places

    # Stringify
    #return stringify(value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits, delimiter=delimiter, **kwargs)[1].strip()
    return stringify(value, **kwargs)[1].strip()

# -----------------------------------------------------------------

#def stringify(value, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, delimiter=",", **kwargs):
def stringify(value, **kwargs):

    """
    This function ...
    :param value:
    :param kwargs:
    :return:
    """

    #print("strdelimiter", delimiter)
    #print("strkwargs", kwargs)

    # List or derived from list
    if isinstance(value, list): return stringify_list(value, **kwargs)

    # Dictionary
    if isinstance(value, dict): return stringify_dict(value, **kwargs)

    # Array or derived from Array, but not quantity!
    #elif isinstance(value, np.ndarray) and not isinstance(value, Quantity):
    #elif introspection.try_importing_module("numpy", True) and (isinstance(value, np.ndarray) and not hasattr(value, "unit")):
    # WE ALSO TEST IF THIS IS NOT A NUMPY INTEGER, FLOAT OR BOOLEAN (because they have a __array__ attribute)
    elif hasattr(value, "__array__") and hasattr(value, "__getitem__") and can_get_item(value) and not hasattr(value, "unit") and not types.is_boolean_type(value) and not types.is_integer_type(value) and not types.is_real_type(value) and not types.is_string_type(value): return stringify_array(value, **kwargs)

    # Column or masked masked column
    elif type(value).__name__ == "MaskedColumn" or type(value).__name__ == "Column": return stringify_array(value, **kwargs)

    # Tuple or derived from tuple
    elif isinstance(value, tuple): return stringify_tuple(value, **kwargs)

    # All other
    #else: return stringify_not_list(value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits, unicode=unicode, **kwargs)
    else: return stringify_not_list(value, **kwargs)

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

def stringify_list(value, **kwargs):

    """
    This function ...
    :param value:
    :param kwargs:
    :return: 
    """

    #print("kwargs", kwargs)

    #if len(value) == 0: raise ValueError("Cannot stringify an empty list")
    if len(value) == 0: return "list", ""

    # If delimiter is passed for stringifying the values in the list
    value_kwargs = copy.copy(kwargs)
    if "value_delimiter" in value_kwargs: value_kwargs["delimiter"] = value_kwargs.pop("value_delimiter")

    strings = []
    ptype = None
    ptypes = set()
    for entry in value:

        #print("Herekwargs", kwargs)

        #parsetype, val = stringify_not_list(entry)
        parsetype, val = stringify(entry, **kwargs)

        if ptype is None: ptype = parsetype
        elif ptype != parsetype:
            #raise ValueError("Nonuniform list")
            ptype = "mixed"

        # Add the parse type
        ptypes.add(parsetype)

        # Add the string
        strings.append(val)

    from ..basics.configuration import parent_type
    from .logging import log

    # Investigate the different ptypes
    parent_types = [parent_type(type_name) for type_name in ptypes]
    #print("Parent types:", parent_types)
    if sequences.all_equal(parent_types) and parent_types[0] is not None: ptype = parent_types[0]
    elif ptype == "mixed": log.warning("Could not determine a common type for '" + stringify(parent_types)[1] + "'")

    # Get delimiter for list
    delimiter = kwargs.pop("delimiter", ",")

    # Return the type and the string
    if ptype.endswith("list"):
        top_delimiter = delimiter + " "
        return ptype + "_list", top_delimiter.join(strings)
    else: return ptype + "_list", delimiter.join(strings)

# -----------------------------------------------------------------

def stringify_dict(value, **kwargs):

    """
    This function ...
    :param value:
    :param kwargs:
    :return: 
    """

    #if len(value) == 0: raise ValueError("Cannot stringify an empty dictionary")
    if len(value) == 0: return "dictionary", ""

    keytype = None
    ptype = None
    parts = []

    keytypes = set()
    ptypes = set()

    # Only for stringifying the values
    value_kwargs = copy.copy(kwargs)
    # If delimiter is passed for stringifying the values in the list
    if "value_delimiter" in value_kwargs: value_kwargs["delimiter"] = value_kwargs.pop("value_delimiter")

    # Loop over the dictionary keys
    for key in value:

        ktype, kstring = stringify(key, **kwargs)

        # Add key type
        keytypes.add(ktype)

        # Check key type
        if keytype is None: keytype = ktype
        elif keytype != ktype: keytype = "mixed"

        v = value[key]

        # Stringify the value
        vtype, vstring = stringify(v, **value_kwargs)

        # Add value type
        ptypes.add(vtype)

        # Check value type
        if ptype is None: ptype = vtype
        elif ptype != vtype: ptype = "mixed"

        # Determine line
        if ptype == "integer" or ptype == "real" or ptype == "boolean": string = "'" + kstring + "': " + vstring
        else: string = "'" + kstring + "': '" + vstring + "'"

        # Add line
        parts.append(string)

    from ..basics.configuration import parent_type
    from .logging import log

    # Investigate the different keytypes
    parent_key_types = [parent_type(type_name) for type_name in keytypes]
    #print("Parent key types:", parent_key_types)
    if sequences.all_equal(parent_key_types) and parent_key_types[0] is not None: ptype = parent_key_types[0]
    elif keytype == "mixed": log.warning("Could not determine a common type for '" + stringify(parent_key_types)[1] + "'")

    # Investigate the different value types
    parent_value_types = [parent_type(type_name) for type_name in ptypes]
    #print("Parent value types:", parent_value_types)
    if sequences.all_equal(parent_value_types) and parent_value_types[0] is not None: ptype = parent_value_types[0]
    elif ptype == "mixed": log.warning("Could not determine a common type for '" + stringify(parent_value_types)[1] + "'")

    # Get delimiter
    delimiter = kwargs.pop("delimiter", ",")

    # Return
    return keytype + "_" + ptype + "_dictionary", delimiter.join(parts)

# -----------------------------------------------------------------

def stringify_array(value, **kwargs):

    """
    This function ...
    :param value:
    :param kwargs:
    :return: 
    """

    # Get delimiter
    delimiter = kwargs.pop("delimiter", ",")

    ptype, val = stringify_not_list(value[0], **kwargs)
    return ptype + "_array", delimiter.join([repr(el) for el in value])

    #ptype, val = stringify_not_list(value[0])
    #return ptype + "_array", ",".join([repr(el) for el in value])

# -----------------------------------------------------------------

def stringify_tuple(value, **kwargs):

    """
    This function ...
    :param value:
    :param kwargs:
    :return: 
    """

    value_kwargs = copy.copy(kwargs)
    if "value_delimiter" in value_kwargs: value_kwargs["delimiter"] = value_kwargs.pop("value_delimiter")

    #print("kwargs", kwargs)

    strings = []
    ptype = None
    for entry in value:

        parsetype, val = stringify_not_list(entry, **kwargs)

        if ptype is None:
            ptype = parsetype
        elif ptype != parsetype:
            raise ValueError("Nonuniform tuple")

        strings.append(val)

    # Get delimiter
    delimiter = kwargs.pop("delimiter", ",")

    # Return
    return ptype + "_tuple", delimiter.join(strings)

# -----------------------------------------------------------------

#def stringify_not_list(value, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, **kwargs):
def stringify_not_list(value, **kwargs):

    """
    This function does stringify, but not for iterables
    :param value:
    :param kwargs:
    :return:
    """

    #print("okwargs", kwargs)

    # Standard
    if types.is_boolean_type(value): return "boolean", str_from_bool(value, **kwargs)
    elif types.is_integer_type(value): return "integer", str_from_integer(value, **kwargs)
    elif types.is_real_type(value): return "real", str_from_real(value, **kwargs)
    elif types.is_string_type(value): return "string", value
    elif types.is_none(value): return "None", "None"

    # Special
    elif introspection.lazy_isinstance(value, "UnitBase", "astropy.units"): return introspection.lazy_call("stringify_unit", "pts.core.units.stringify", value, **kwargs)
    elif introspection.lazy_isinstance(value, "Quantity", "astropy.units"): return introspection.lazy_call("stringify_quantity", "pts.core.units.stringify", value, **kwargs)
    elif introspection.lazy_isinstance(value, "Angle", "astropy.coordinates"): return "angle", str_from_angle(value, **kwargs)
    elif introspection.lazy_isinstance(value, "RealRange", "pts.core.basics.range"): return "real_range", str_from_real_range(value, **kwargs)
    elif introspection.lazy_isinstance(value, "IntegerRange", "pts.core.basics.range"): return "integer_range", str_from_integer_range(value, **kwargs)
    elif introspection.lazy_isinstance(value, "QuantityRange", "pts.core.basics.range"): return "quantity_range", introspection.lazy_call("str_from_quantity_range", "pts.core.units.stringify", value, **kwargs)
    elif introspection.lazy_isinstance(value, "SkyCoordinate", "pts.magic.basics.coordinate"): return "skycoordinate", repr(value.ra.value) + " " + str(value.ra.unit) + "," + repr(value.dec.value) + " " + str(value.dec.unit)
    elif introspection.lazy_isinstance(value, "SkyStretch", "pts.magic.basics.stretch"): return "skystretch", repr(value.ra.value) + " " + str(value.ra.unit) + "," + repr(value.dec.value) + " " + str(value.dec.unit)
    elif introspection.lazy_isinstance(value, "Filter", "pts.core.filter.filter"): return introspection.lazy_call("stringify_filter", "pts.core.filter.filter", value, **kwargs)

    # Other
    #elif introspection.isinstance(Instrument):

    # Unrecognized
    else: raise ValueError("Unrecognized type: " + str(type(value)))

# -----------------------------------------------------------------

def str_from_dictionary(dictionary, **kwargs):

    """
    This function ...
    :param dictionary:
    :param kwargs:
    :return:
    """

    parts = []
    for key in dictionary:
        value = dictionary[key]
        vtype, vstring = stringify(value, **kwargs)
        string = key + ": " + vstring
        parts.append(string)
    return ",".join(parts)

# -----------------------------------------------------------------

def stringify_string_fancy(string, **kwargs):

    """
    This function ...
    :param string:
    :return:
    """

    width = kwargs.pop("width", 100)
    lines_prefix = kwargs.pop("lines_prefix", "")

    from textwrap import wrap
    return "string", lines_prefix + ("\n" + lines_prefix).join(wrap(string, width))

# -----------------------------------------------------------------

def stringify_list_fancy(lst, **kwargs):

    """
    This function ...
    :param lst:
    :param kwargs:
    :return:
    """

    width = kwargs.pop("width", 100)
    delimiter = kwargs.pop("delimiter", ", ")
    lines_prefix = kwargs.pop("lines_prefix", "")

    from textwrap import wrap

    ptype, string = stringify(lst)
    return ptype, lines_prefix + ("\n" + lines_prefix).join(wrap(string.replace(",", delimiter), width))

# -----------------------------------------------------------------

def stringify_paths(paths, **kwargs):

    """
    This function ...
    :param paths:
    :param kwargs:
    :return:
    """

    # Get options
    base = kwargs.pop("basse", None)

    if base is None: return "path_list", stringify_list(paths)[1]
    else:

        from . import filesystem as fs
        absolute_base = fs.absolute_path(base)

        # Return the type and the relative paths as a string list
        return "string_list", stringify_list([fs.absolute_path(path).split(absolute_base)[1] for path in paths])[1]

# -----------------------------------------------------------------

def str_from_integer(integer, **kwargs):

    """
    This function ...
    :param integer:
    :param kwargs:
    :return:
    """
    
    # Get settings
    scientific = kwargs.pop("scientific", False)
    decimal_places = kwargs.pop("decimal_places", 2)
    fancy = kwargs.pop("fancy", False)
    ndigits = kwargs.pop("ndigits", None)
    unicode = kwargs.pop("unicode", False)

    # Check input
    if ndigits is not None and ndigits < 1: raise ValueError("Number of digits cannot be smaller than 1")

    if scientific:
        if fancy:
            if ndigits is not None:
                power = len(str(integer)) - 1
                digits = []
                str_rounded = str(integer)
                for i in range(ndigits):
                    digit = str_rounded[i]
                    digits.append(digit)
                if unicode: return digits[0].decode("utf8") + u"." + u"".join(digits[1:]) + u" " + strings.multiplication + u" 10" + strings.superscript(power) # DOESN'T WORK??
                else: return digits[0] + "." + "".join(digits[1:]) + " x 10^" + str(power)
            else:
                result = "{:.0e}".format(integer).replace("+", "").replace("e0", "e")
                power = int(result.split("e")[1])
                if unicode: result = result.split("e")[0].decode("utf8") + u" " + strings.multiplication + u" 10" + strings.superscript(power) # DOESN'T WORK
                else: result = result.split("e")[0] + " x 10^" + str(power)
                return result
        else:
            if ndigits is not None: decimal_places = ndigits - 1
            #return "{:.0e}".format(integer).replace("+", "").replace("e0", "e")
            return ("{:." + str(decimal_places) + "e}").format(float(integer)).replace("+", "").replace("e0", "e")
    else: return str(integer)

# -----------------------------------------------------------------

#def str_from_integer_range(the_range, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, **kwargs):
def str_from_integer_range(the_range, **kwargs):

    """
    Thi function ...
    :param the_range:
    :param kwargs:
    :return:
    """

    min_str = str_from_integer(the_range.min, **kwargs)
    max_str = str_from_integer(the_range.max, **kwargs)

    return min_str + " > " + max_str

# -----------------------------------------------------------------

#def str_from_real(real, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, **kwargs):
def str_from_real(real, **kwargs):

    """
    This function ...
    :param real:
    :param kwargs:
    :return:
    """

    # Get kwargs
    scientific = kwargs.pop("scientific", False)
    decimal_places = kwargs.pop("decimal_places", 2)
    fancy = kwargs.pop("fancy", False)
    ndigits = kwargs.pop("ndigits", None)
    unicode = kwargs.pop("unicode", False)

    # Check input
    if ndigits is not None and ndigits < 1: raise ValueError("Number of digits cannot be smaller than 1")

    #print(real)

    if scientific:
        if fancy:
            if ndigits is not None:
                power = len(str(real).split(".")[0]) - 1
                digits = []
                rounded = numbers.round_to_n_significant_digits(real, ndigits)
                str_rounded = str(rounded)
                #print(str_rounded)
                #if "." in str_rounded: enditeration = ndigits + 1
                #else: enditeration = ndigits
                if "." in str_rounded: str_rounded = "".join(str_rounded.split("."))
                for i in range(ndigits):
                    digit = str_rounded[i]
                    #if digit == ".": continue # happens if rounded does stil contain dot
                    digits.append(digit)
                #print("digits", digits)
                if unicode: return digits[0].decode("utf8") + u"." + u"".join(digits[1:]) + u" " + strings.multiplication + u" 10" + strings.superscript(power).decode("utf8") # DOESN'T WORK??
                else: return digits[0] + "." + "".join(digits[1:]) + " x 10^" + str(power)
            else:
                result = ("{:." + str(decimal_places) + "e}").format(real).replace("+", "").replace("e0", "e")
                power = int(result.split("e")[1])
                #result = result.split("e")[0].decode("utf8") + u" " + strings.multiplication + u" 10" + strings.superscript(power).decode("utf8")
                #result = result.split("e")[0].decode("utf8") + u" " + strings.multiplication + u" 10" + strings.superscript(power).decode("utf8")
                if unicode: result = result.split("e")[0].decode("utf8") + u" " + u"x" + u" 10" + strings.superscript(power).decode("utf8") # SOMETHING LIKE THIS?? DOESN'T WORK??
                else: result = result.split("e")[0] + " x 10^" + str(power)
                return result
        else:
            if ndigits is not None: decimal_places = ndigits - 1
            return ("{:." + str(decimal_places) + "e}").format(real).replace("+", "").replace("e0", "e")

    else: return repr(real)

# -----------------------------------------------------------------

#def str_from_real_range(the_range, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False, **kwargs):
def str_from_real_range(the_range, **kwargs):

    """
    This function ...
    :param the_range:
    :param kwargs:
    :return:
    """

    min_str = str_from_real(the_range.min, **kwargs)
    max_str = str_from_real(the_range.max, **kwargs)

    return min_str + " > " + max_str

# -----------------------------------------------------------------

def yes_or_no(boolean, **kwargs):

    """
    This function ...
    :param boolean:
    :param kwargs:
    :return:
    """

    # Get options
    short = kwargs.pop("short", False)

    answer = "yes" if boolean else "no"
    if short: return answer[0]
    else: return answer

# -----------------------------------------------------------------

def str_from_bool(boolean, **kwargs):

    """
    This function ...
    :param boolean:
    :param kwargs:
    :return:
    """

    # Get options
    lower = kwargs.pop("lower", False)

    if lower: return str(boolean).lower()
    else: return str(boolean)

# -----------------------------------------------------------------

def str_from_angle(angle, **kwargs):

    """
    This function ...
    :param angle:
    :param kwargs:
    :return:
    """

    return str_from_real(angle.value, **kwargs) + " " + str(angle.unit).replace(" ", "")

# -----------------------------------------------------------------
