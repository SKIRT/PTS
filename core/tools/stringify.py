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
from types import NoneType

# Import the relevant PTS classes and modules
from . import types
from . import introspection

# -----------------------------------------------------------------

def stringify(value, scientific=False, decimal_places=2):

    """
    This function ...
    :param value:
    :param scientific:
    :param decimal_places:
    :return:
    """

    # List or derived from list
    if isinstance(value, list):

        strings = []
        ptype = None
        for entry in value:

            #parsetype, val = stringify_not_list(entry)
            parsetype, val = stringify(entry)

            if ptype is None: ptype = parsetype
            elif ptype != parsetype:
                #raise ValueError("Nonuniform list")
                ptype = "mixed"

            strings.append(val)

        return ptype + "_list", ",".join(strings)

    # Array or derived from Array, but not quantity
    #elif isinstance(value, np.ndarray) and not isinstance(value, Quantity):
    #elif introspection.try_importing_module("numpy", True) and (isinstance(value, np.ndarray) and not hasattr(value, "unit")):
    elif hasattr(value, "__array__") and not hasattr(value, "unit"):

        ptype, val = stringify_not_list(value[0])
        return ptype + "_array", ",".join([repr(el) for el in value])

    # Tuple or derived from tuple
    elif isinstance(value, tuple):

        strings = []
        ptype = None
        for entry in value:

            parsetype, val = stringify_not_list(entry)

            if ptype is None: ptype = parsetype
            elif ptype != parsetype: raise ValueError("Nonuniform tuple")

            strings.append(val)

        return ptype + "_tuple", ",".join(strings)

    # All other
    else: return stringify_not_list(value, scientific=scientific, decimal_places=decimal_places)

# -----------------------------------------------------------------

def stringify_not_list(value, scientific=False, decimal_places=2):

    """
    This function ...
    :param value:
    :param scientific:
    :param decimal_places:
    :return:
    """

    #from pts.magic.basics.coordinate import SkyCoordinate
    #from pts.magic.basics.stretch import SkyStretch

    #from astropy.units import UnitBase, Quantity
    #from astropy.coordinates import Angle

    #from ..filter.broad import BroadBandFilter
    #from ..filter.narrow import NarrowBandFilter

    #from ..basics.range import RealRange, IntegerRange, QuantityRange

    #from ..basics.unit import stringify_unit
    #from ..basics.quantity import stringify_quantity

    if types.is_boolean_type(value): return "boolean", str_from_bool(value)

    elif types.is_integer_type(value): return "integer", str_from_integer(value, scientific=scientific)

    elif types.is_real_type(value): return "real", str_from_real(value, scientific=scientific, decimal_places=decimal_places)

    elif isinstance(value, basestring): return "string", value

    elif isinstance(value, NoneType): return "None", "None"

    elif introspection.try_importing_class("UnitBase", "astropy.units", True) and isinstance(value, UnitBase): return stringify_unit(value)

    elif introspection.try_importing_class("Quantity", "astropy.units", True) and isinstance(value, Quantity): return stringify_quantity(value)

    elif introspection.try_importing_class("Angle", "astropy.coordinates", True) and isinstance(value, Angle): return "angle", str_from_angle(value)

    elif introspection.try_importing_class("RealRange", "pts.core.basics.range", True) and isinstance(value, RealRange): return "real_range", repr(value)
    elif introspection.try_importing_class("IntegerRange", "pts.core.basics.range", True) and isinstance(value, IntegerRange): return "integer_range", repr(value)
    elif introspection.try_importing_class("QuantityRange", "pts.core.basics.range", True) and isinstance(value, QuantityRange): return "quantity_range", repr(value)

    elif introspection.try_importing_class("SkyCoordinate", "pts.magic.basics.coordinate", True) and isinstance(value, SkyCoordinate): return "skycoordinate", repr(value.ra.value) + " " + str(value.ra.unit) + "," + repr(value.dec.value) + " " + str(value.dec.unit)

    elif introspection.try_importing_class("SkyStretch", "pts.magic.basics.stretch", True) and isinstance(value, SkyStretch): return "skystretch", repr(value.ra.value) + " " + str(value.ra.unit) + "," + repr(value.dec.value) + " " + str(value.dec.unit)

    elif introspection.try_importing_class("NarrowBandFilter", "pts.core.filter.narrow", True) and isinstance(value, NarrowBandFilter): return "narrow_band_filter", str(value)

    elif introspection.try_importing_class("BroadBandFilter", "pts.core.filter.broad", True) and isinstance(value, BroadBandFilter): return "broad_band_filter", str(value)

    else: raise ValueError("Unrecognized type: " + str(type(value)))

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

def str_from_integer(integer, scientific=False):

    """
    This function ...
    :param integer:
    :param scientific:
    :return:
    """

    if scientific: return "{:.0e}".format(integer).replace("+", "").replace("e0", "e")
    else: return str(integer)

# -----------------------------------------------------------------

def str_from_real(real, scientific=False, decimal_places=2):

    """
    This function ...
    :param real:
    :param scientific:
    :param decimal_places:
    :return:
    """

    if scientific: return ("{:." + str(decimal_places) + "e}").format(real).replace("+", "").replace("e0", "e")
    else: return repr(real)

# -----------------------------------------------------------------

def yes_or_no(boolean, short=False):

    """
    This function ...
    :param boolean:
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

    return "angle", repr(angle.value) + " " + str(angle.unit).replace(" ", "")

# -----------------------------------------------------------------
