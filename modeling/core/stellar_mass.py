#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.stellar_mass Contains stellar mass estimators.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...magic.core.frame import Frame
from ...core.units.parsing import parse_unit as u
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

def hubble_stage_to_type(stage, add_subtype=False):

    """
    # SOURCE: https://en.wikipedia.org/wiki/Galaxy_morphological_classification
    This function ...
    :param stage:
    :param add_subtype:
    :return:
    """

    if stage < -3.5:

        if add_subtype:

            if stage < -5.5: subtype = "E-"
            elif stage < -4.5: subtype = "E"
            else: subtype = "E+"
            return "E", subtype

        else: return "E"

    elif stage < -0.5:

        if add_subtype:

            if stage < -2.5: subtype = "S0-"
            elif stage < -1.5: subtype = "S00"
            else: subtype = "S0+"
            return "S0", subtype

        else: return "S0"

    elif stage < 0.5:

        if add_subtype: return "S0/a", None
        else: return "S0/a"

    elif stage < 1.5:

        if add_subtype: return "Sa", None
        else: return "Sa"

    elif stage < 2.5:

        if add_subtype: return "Sab", None
        else: return "Sab"

    elif stage < 3.5:

        if add_subtype: return "Sb", None
        else: return "Sb"

    elif stage < 4.5:

        if add_subtype: return "Sbc", None
        else: return "Sbc"

    elif stage < 7.5:

        if add_subtype:
            if stage < 5.5: subtype = "Sc"
            elif stage < 6.5: subtype = "Scd"
            else: subtype = "Sd"
            return "Sc", subtype
        else: return "Sc"

    elif stage < 8.5:

        if add_subtype: return "Sc/Irr", None
        else: return "Sc/Irr"

    else:

        if add_subtype:
            if stage < 9.5: subtype = "Sm"
            elif stage < 10.5: subtype = "Im"
            else: subtype = None
            return "Irr", subtype
        else: return "Irr"

# -----------------------------------------------------------------

def hubble_type_to_stage(hubble_class):

    """
    This function ...
    # SOURCE: https://en.wikipedia.org/wiki/Galaxy_morphological_classification
    :param hubble_class:
    :return:
    """

    # E
    if hubble_class == "E": return -5

    # SUBDIVISION (VAUCOULEURS)
    elif hubble_class == "cE": return -6
    elif hubble_class == "E+": return -4

    # S0
    elif hubble_class == "S0": return -2

    # SUBDIVISION (VAUCOULEURS)
    elif hubble_class == "S0-": return -3
    elif hubble_class == "S00": return -2
    elif hubble_class == "S0+": return -1

    # S0/a
    elif hubble_class == "S0/a": return 0

    # Sa
    elif hubble_class == "Sa": return 1

    # Sab
    elif hubble_class == "Sab": return 2
    # synonym
    elif hubble_class == "Sa-b": return 2

    # Sb
    elif hubble_class == "Sb": return 3

    # Sbc
    elif hubble_class == "Sbc": return 4
    elif hubble_class == "Sb-c": return 4

    # Sc
    elif hubble_class == "Sc": return 5.5
    elif hubble_class == "Scd": return 6
    elif hubble_class == "Sc-d": return 6
    elif hubble_class == "Sd": return 7

    # Sc/Irr
    elif hubble_class == "Sc/Irr": return 8
    elif hubble_class == "Sc-Irr": return 8
    elif hubble_class == "Sdm": return 8

    # Higher
    elif hubble_class == "Sm": return 9
    elif hubble_class == "Irr": return 9.5
    elif hubble_class == "Im": return 10
    elif hubble_class == "I": return 9.5

    # Invalid
    else: raise ValueError("Unknown hubble class: '" + hubble_class + "'")

# -----------------------------------------------------------------

# OLIVER 2010
# (M∗/M⊙)/[νLν(3.6)/L⊙] to be 38.4, 40.8, 27.6, 35.3, 18.7 and 26.7, for types E, Sab, Sbc, Scd, Sdm and sb
# measuring the 3.6 μm monochromatic luminosity in total solar units, not in units of the Sun’s monochromatic 3.6 μm lu- minosity

oliver_stellar_mass_factors = OrderedDict()
oliver_stellar_mass_factors["E"] = 38.4
oliver_stellar_mass_factors["Sab"] = 40.8
oliver_stellar_mass_factors["Sb"] = 26.7
oliver_stellar_mass_factors["Sbc"] = 27.6
oliver_stellar_mass_factors["Scd"] = 35.3
oliver_stellar_mass_factors["Sdm"] = 18.7

# -----------------------------------------------------------------

def get_oliver_stellar_mass_factor(hubble_type, hubble_subtype=None):

    """
    Thisf unction ...
    :param hubble_type:
    :param hubble_subtype:
    :return:
    """

    # Get the factor
    if hubble_type not in oliver_stellar_mass_factors:
        if hubble_subtype is not None:
            if hubble_subtype not in oliver_stellar_mass_factors: raise ValueError("Hubble type '" + hubble_type + "' or '" + hubble_subtype + "' not supported")
            else: factor = oliver_stellar_mass_factors[hubble_subtype]
        else: raise ValueError("Hubble type '" + hubble_type + "' not supported")
    else: factor = oliver_stellar_mass_factors[hubble_type]

    # Return the factor
    return factor

# -----------------------------------------------------------------

def oliver_stellar_mass(i1_luminosity, hubble_type, hubble_subtype=None, distance=None):

    """
    This function ...
    :param i1_luminosity:
    :param hubble_type:
    :param hubble_subtype:
    :param distance:
    :return:
    """

    from ..core.data import Data3D

    # Get the I1 wavelength
    i1_wavelength = parse_filter("IRAC I1").wavelength

    # Get the factor
    oliver_factor = get_oliver_stellar_mass_factor(hubble_type, hubble_subtype=hubble_subtype)

    # Frame
    if isinstance(i1_luminosity, Frame):

        converted = i1_luminosity.converted_to("Lsun", density=True, wavelength=i1_wavelength, distance=distance)
        converted *= oliver_factor
        converted.unit = "Msun"
        return converted

    # 3D data
    elif isinstance(i1_luminosity, Data3D):

        factor = i1_luminosity.unit.conversion_factor("Lsun", density=True, wavelength=i1_wavelength, distance=distance)
        factor *= oliver_factor
        return i1_luminosity.converted_by_factor(factor, "Msun", new_name="Mstar", new_description="Stellar mass (Oliver)")

    # Photometric quantity
    else: return i1_luminosity.to("Lsun", density=True, wavelength=i1_wavelength, distance=distance).value * oliver_factor * u("Msun")
    
# -----------------------------------------------------------------
