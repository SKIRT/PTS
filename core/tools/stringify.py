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

# Import the relevant PTS classes and modules
from . import types
from . import introspection
from . import sequences
from . import strings
from . import numbers

# -----------------------------------------------------------------

def tostr(value, scientific=None, decimal_places=2, fancy=False, ndigits=None):

    """
    This function ...
    :param value: 
    :param scientific:
    :param decimal_places:
    :param fancy:
    :param ndigits:
    :return: 
    """

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

        # Other
        else: scientific = False

    # Stringify
    return stringify(value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits)[1].strip()

# -----------------------------------------------------------------

def stringify(value, scientific=False, decimal_places=2, fancy=False, ndigits=None):

    """
    This function ...
    :param value:
    :param scientific:
    :param decimal_places:
    :param fancy:
    :param ndigits:
    :return:
    """

    # List or derived from list
    if isinstance(value, list): return stringify_list(value)

    # Dictionary
    if isinstance(value, dict): return stringify_dict(value)

    # Array or derived from Array, but not quantity!
    #elif isinstance(value, np.ndarray) and not isinstance(value, Quantity):
    #elif introspection.try_importing_module("numpy", True) and (isinstance(value, np.ndarray) and not hasattr(value, "unit")):
    # WE ALSO TEST IF THIS IS NOT A NUMPY INTEGER, FLOAT OR BOOLEAN (because they have a __array__ attribute)
    elif hasattr(value, "__array__") and hasattr(value, "__getitem__") and can_get_item(value) and not hasattr(value, "unit") and not types.is_boolean_type(value) and not types.is_integer_type(value) and not types.is_real_type(value) and not types.is_string_type(value): return stringify_array(value)

    # Column or masked masked column
    elif type(value).__name__ == "MaskedColumn" or type(value).__name__ == "Column": return stringify_array(value)

    # Tuple or derived from tuple
    elif isinstance(value, tuple): return stringify_tuple(value)

    # All other
    else: return stringify_not_list(value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits)

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

def stringify_list(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    #if len(value) == 0: raise ValueError("Cannot stringify an empty list")
    if len(value) == 0: return "list", ""

    strings = []
    ptype = None
    ptypes = set()
    for entry in value:

        #parsetype, val = stringify_not_list(entry)
        parsetype, val = stringify(entry)

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

    # Return the type and the string
    return ptype + "_list", ",".join(strings)

# -----------------------------------------------------------------

def stringify_dict(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    #if len(value) == 0: raise ValueError("Cannot stringify an empty dictionary")
    if len(value) == 0: return "dictionary", ""

    keytype = None
    ptype = None
    parts = []

    keytypes = set()
    ptypes = set()

    for key in value:

        ktype, kstring = stringify(key)

        # Add key type
        keytypes.add(ktype)

        # Check key type
        if keytype is None: keytype = ktype
        elif keytype != ktype: keytype = "mixed"

        v = value[key]

        vtype, vstring = stringify(v)

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

    # Return
    return keytype + "_" + ptype + "_dictionary", ",".join(parts)

# -----------------------------------------------------------------

def stringify_array(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    ptype, val = stringify_not_list(value[0])
    return ptype + "_array", ",".join([repr(el) for el in value])

    #ptype, val = stringify_not_list(value[0])
    #return ptype + "_array", ",".join([repr(el) for el in value])

# -----------------------------------------------------------------

def stringify_tuple(value):

    """
    This function ...
    :param value: 
    :return: 
    """

    strings = []
    ptype = None
    for entry in value:

        parsetype, val = stringify_not_list(entry)

        if ptype is None:
            ptype = parsetype
        elif ptype != parsetype:
            raise ValueError("Nonuniform tuple")

        strings.append(val)

    return ptype + "_tuple", ",".join(strings)

# -----------------------------------------------------------------

def stringify_not_list(value, scientific=False, decimal_places=2, fancy=False, ndigits=None):

    """
    This function does stringify, but not for iterables
    :param value:
    :param scientific:
    :param decimal_places:
    :param fancy:
    :param ndigits:
    :return:
    """

    # Standard
    if types.is_boolean_type(value): return "boolean", str_from_bool(value)
    elif types.is_integer_type(value): return "integer", str_from_integer(value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits)
    elif types.is_real_type(value): return "real", str_from_real(value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits)
    elif types.is_string_type(value): return "string", value
    elif types.is_none(value): return "None", "None"

    # Special
    elif introspection.lazy_isinstance(value, "UnitBase", "astropy.units"): return introspection.lazy_call("stringify_unit", "pts.core.units.stringify", value)
    elif introspection.lazy_isinstance(value, "Quantity", "astropy.units"): return introspection.lazy_call("stringify_quantity", "pts.core.units.stringify", value, scientific=scientific, decimal_places=decimal_places, fancy=fancy, ndigits=ndigits)
    elif introspection.lazy_isinstance(value, "Angle", "astropy.coordinates"): return "angle", str_from_angle(value)
    elif introspection.lazy_isinstance(value, "RealRange", "pts.core.basics.range"): return "real_range", repr(value)
    elif introspection.lazy_isinstance(value, "IntegerRange", "pts.core.basics.range"): return "integer_range", repr(value)
    elif introspection.lazy_isinstance(value, "QuantityRange", "pts.core.basics.range"): return "quantity_range", repr(value)
    elif introspection.lazy_isinstance(value, "SkyCoordinate", "pts.magic.basics.coordinate"): return "skycoordinate", repr(value.ra.value) + " " + str(value.ra.unit) + "," + repr(value.dec.value) + " " + str(value.dec.unit)
    elif introspection.lazy_isinstance(value, "SkyStretch", "pts.magic.basics.stretch"): return "skystretch", repr(value.ra.value) + " " + str(value.ra.unit) + "," + repr(value.dec.value) + " " + str(value.dec.unit)
    elif introspection.lazy_isinstance(value, "NarrowBandFilter", "pts.core.filter.narrow"): return "narrow_band_filter", str(value)
    elif introspection.lazy_isinstance(value, "BroadBandFilter", "pts.core.filter.broad"): return "broad_band_filter", str(value)

    # Other
    #elif introspection.isinstance(Instrument):

    # Unrecognized
    else: raise ValueError("Unrecognized type: " + str(type(value)))

# -----------------------------------------------------------------

def str_from_dictionary(dictionary):

    """
    This function ...
    :param dictionary:
    :return:
    """

    parts = []
    for key in dictionary:
        value = dictionary[key]
        vtype, vstring = stringify(value)
        string = key + ": " + vstring
        parts.append(string)
    return ",".join(parts)

# -----------------------------------------------------------------

def stringify_string_fancy(string, width=100, lines_prefix=""):

    """
    This function ...
    :param string:
    :param width:
    :param lines_prefix:
    :return:
    """

    from textwrap import wrap
    return "string", lines_prefix + ("\n" + lines_prefix).join(wrap(string, width))

# -----------------------------------------------------------------

def stringify_list_fancy(lst, width=100, delimiter=", ", lines_prefix=""):

    """
    This function ...
    :param lst:
    :param width:
    :param delimiter:
    :param lines_prefix:
    :return:
    """

    from textwrap import wrap

    ptype, string = stringify(lst)
    return ptype, lines_prefix + ("\n" + lines_prefix).join(wrap(string.replace(",", delimiter), width))

# -----------------------------------------------------------------

def stringify_paths(paths, base=None):

    """
    This function ...
    :param paths:
    :param delimiter:
    :param base:
    :return:
    """

    if base is None: return "path_list", stringify_list(paths)[1]
    else:

        from . import filesystem as fs
        absolute_base = fs.absolute_path(base)

        # Return the type and the relative paths as a string list
        return "string_list", stringify_list([fs.absolute_path(path).split(absolute_base)[1] for path in paths])

# -----------------------------------------------------------------

def str_from_integer(integer, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False):

    """
    This function ...
    :param integer:
    :param scientific:
    :param fancy:
    :param ndigits:
    :param unicode:
    :return:
    """

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

def str_from_real(real, scientific=False, decimal_places=2, fancy=False, ndigits=None, unicode=False):

    """
    This function ...
    :param real:
    :param scientific:
    :param decimal_places:
    :param fancy:
    :param ndigits:
    :return:
    """

    # Check input
    if ndigits is not None and ndigits < 1: raise ValueError("Number of digits cannot be smaller than 1")

    if scientific:
        if fancy:
            if ndigits is not None:
                power = len(str(real).split(".")[0]) - 1
                digits = []
                rounded = numbers.round_to_n_significant_digits(real, ndigits)
                str_rounded = str(rounded)
                for i in range(ndigits):
                    digit = str_rounded[i]
                    digits.append(digit)
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

def yes_or_no(boolean, short=False):

    """
    This function ...
    :param boolean:
    :param short:
    :return:
    """

    answer = "yes" if boolean else "no"
    if short: return answer[0]
    else: return answer

# -----------------------------------------------------------------

def str_from_bool(boolean, lower=False):

    """
    This function ...
    :param boolean:
    :param lower:
    :return:
    """

    if lower: return str(boolean).lower()
    else: return str(boolean)

# -----------------------------------------------------------------

def str_from_angle(angle):

    """
    This function ...
    :param angle:
    :return:
    """

    return repr(angle.value) + " " + str(angle.unit).replace(" ", "")

# -----------------------------------------------------------------
