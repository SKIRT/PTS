#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.sfr Contains SFR estimators.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from ...core.filter.filter import parse_filter
from ...magic.core.frame import Frame
from ...core.units.parsing import parse_unit as u
from ...core.units.unit import get_converted_value
from ...core.tools import types

# -----------------------------------------------------------------

# Kennicutt & Evans 2012
kennicutt_evans_logc = 43.35

# Salim 20??
salim = 1.08e-28

# Hao 2011
hao = 4.55e-44

# Calzetti 2010
calzetti_a = 1.70e-43
calzetti_b = 2.03e-44
calzetti_exponent = 0.048

# Calzetti 2007
calzetti_local = 1.31e-38
calzetti_local_exponent = 0.885

# Kennicutt 1998
kennicutt = 1.4e-28 # FUV
kennicutt_tir = 4.5e-44 # TIR

# -----------------------------------------------------------------

def converted_by_factor(luminosity, conversion, new_unit, required_unit, unit=None, wavelength=None, distance=None,
                        name=None, description=None):

    """
    This function ...
    :param luminosity:
    :param factor:
    :param new_unit:
    :param required_unit:
    :param unit:
    :param wavelength:
    :param distance:
    :param name:
    :param description:
    :return:
    """

    from ..core.data import Data3D

    # Frame
    if isinstance(luminosity, Frame):

        if unit is not None: raise ValueError("Cannot specify unit")

        converted = luminosity.converted_to(required_unit, wavelength=wavelength, distance=distance)
        converted *= conversion
        converted.unit = new_unit
        return converted

    # 3D data
    elif isinstance(luminosity, Data3D):

        if unit is not None: raise ValueError("Cannot specify unit")

        factor = luminosity.unit.conversion_factor(required_unit, wavelength=wavelength, distance=distance)
        factor *= conversion
        return luminosity.converted_by_factor(factor, new_unit, new_name=name, new_description=description)

    # Photometric quantity
    elif types.is_quantity(luminosity):
        if unit is not None: raise ValueError("Cannot specify unit")
        return luminosity.to(required_unit, wavelength=wavelength, distance=distance).value * conversion * new_unit

    # Array
    elif types.is_array_like(luminosity):
        if unit is None: raise ValueError("Unit is not specified")
        factor = unit.conversion_factor(required_unit, wavelength=wavelength, distance=distance)
        factor *= conversion
        return luminosity * factor

    # Invalid
    else: raise ValueError("Invalid type for 'luminosity'")

# -----------------------------------------------------------------

def estimate_fuv_to_sfr(estimator, factor, required_unit, fuv_luminosity, unit=None, distance=None):

    """
    This function ...
    :param estimator:
    :param factor:
    :param required_unit:
    :param fuv_luminosity:
    :param unit:
    :param distance:
    :return:
    """

    # Get the FUV wavelength
    fuv_wavelength = parse_filter("GALEX FUV").wavelength
    new_unit = u("Msun/yr")
    required_unit = u(required_unit, density=True)

    # Set name and description
    name = "SFR"
    description = "star formation rate (FUV [" + estimator + "])"

    # Return converted to SFR (various different types)
    return converted_by_factor(fuv_luminosity, factor, new_unit, required_unit, unit=unit, wavelength=fuv_wavelength, distance=distance, name=name, description=description)

# -----------------------------------------------------------------

def estimate_tir_to_sfr(estimator, factor, required_unit, tir_luminosity, unit=None, distance=None):

    """
    This function ...
    :param estimator:
    :param factor:
    :param required_unit:
    :param fuv_luminosity:
    :param unit:
    :param distance:
    :return:
    """

    # Get the
    new_unit = u("Msun/yr")
    required_unit = u(required_unit)

    # Set name and description
    name = "SFR"
    description = "star formation rate (TIR [" + estimator + "])"

    # Return converted to SFR (various different types)
    return converted_by_factor(tir_luminosity, factor, new_unit, required_unit, unit=unit, distance=distance, name=name, description=description)

# -----------------------------------------------------------------

def kennicutt_evans_fuv_to_sfr(fuv_luminosity, unit=None, distance=None):

    """
    This function ...
    :param fuv_luminosity:
    :param unit:
    :param distance:
    :return:
    """

    # Calculate factor
    calibration = 1./ 10**kennicutt_evans_logc

    # Calculate and return
    return estimate_fuv_to_sfr("Kennicutt & Evans", calibration, "erg/s", fuv_luminosity, unit=unit, distance=distance)

# -----------------------------------------------------------------

def kennicutt_fuv_to_sfr(fuv_luminosity, unit=None, distance=None):

    """
    This function ...
    :param fuv_luminosity:
    :param unit:
    :param distance:
    :return:
    """

    return estimate_fuv_to_sfr("Kennicutt", kennicutt, "erg/s/Hz", fuv_luminosity, unit=unit, distance=distance)

# -----------------------------------------------------------------

def salim_fuv_to_sfr(fuv_luminosity, unit=None, distance=None):

    """
    This function ...
    :param fuv_luminosity:
    :param unit:
    :param distance:
    :return:
    """

    return estimate_fuv_to_sfr("Salim", salim, "erg/s/Hz", fuv_luminosity, unit=unit, distance=distance)

# -----------------------------------------------------------------

def hao_fuv_to_sfr(fuv_luminosity, unit=None, distance=None):

    """
    This function ...
    :param fuv_luminosity:
    :param unit:
    :param distance:
    :return:
    """

    return estimate_fuv_to_sfr("Hao", hao, "erg/s", fuv_luminosity, unit=unit, distance=distance)

# -----------------------------------------------------------------

def calzetti_24um_to_sfr(luminosity, unit=None, distance=None):

    """
    SFR = 1.70 × 10−43L(24) × [2.03 × 10−44L(24)] 0.048
    :param luminosity:
    :param unit:
    :param distance:
    :return:
    """

    # Condition:
    # if L(24) ≥ 5 x 1043​erg s−1

    # Set name and description
    name = "SFR"
    description = "star formation rate (24 micron [Calzetti])"

    # Set new unit
    new_unit = u("Msun/yr")

    # Get wavelength
    wavelength = parse_filter("MIPS 24mu").wavelength

    # Set required unit
    required_unit = u("erg/s", density=True)

    from ..core.data import Data3D

    # Frame
    if isinstance(luminosity, Frame):

        # LOCAL

        if unit is not None: raise ValueError("Cannot specify unit")

        # Get converted frame
        converted = luminosity.converted_to(required_unit, wavelength=wavelength, distance=distance)

        # Return
        #new_array = calzetti_a * converted.data * (calzetti_b * converted.data)**calzetti_exponent
        new_array = calzetti_local * converted.data**calzetti_local_exponent
        frame = Frame(new_array, name=name, description=description, unit=new_unit, wcs=luminosity.wcs,
                     pixelscale=luminosity.pixelscale, psf_filter=luminosity.psf_filter, fwhm=luminosity.fwhm,
                      distance=luminosity.distance)
        if distance is not None: frame.distance = distance
        return frame

    # 3D data
    elif isinstance(luminosity, Data3D):

        # LOCAL

        if unit is not None: raise ValueError("Cannot specify unit")

        # Get converted data
        converted = luminosity.converted_to(required_unit, wavelength=wavelength, distance=distance)

        # Calculate and return
        #new_array = calzetti_a * converted.values * (calzetti_b * converted.values)**calzetti_exponent
        new_array = calzetti_local * converted.values**calzetti_local_exponent
        return Data3D.from_other(luminosity, values=new_array, unit=new_unit, name=name, description=description)

    # Photometric quantity
    elif types.is_quantity(luminosity):

        # GLOBAL

        # Check
        if unit is not None: raise ValueError("Cannot specify unit")

        # Get converted value
        value = get_converted_value(luminosity, required_unit, distance=distance, wavelength=wavelength)

        # Calculate and return
        new_value = calzetti_a * value * (calzetti_b * value)**calzetti_exponent
        return Quantity(new_value, new_unit)

    # Array
    elif types.is_array_like(luminosity):

        # GLOBAL

        # Check
        if unit is None: raise ValueError("Unit is not specified")

        # Create converted array
        factor = unit.conversion_factor(required_unit, wavelength=wavelength, distance=distance)
        converted = luminosity * factor

        # Return
        return calzetti_a * converted * (calzetti_b * converted)**calzetti_exponent

    # Invalid
    else: raise ValueError("Invalid type for 'luminosity'")

# -----------------------------------------------------------------

def kennicutt_tir_to_sfr(luminosity, unit=None, distance=None):

    """
    This function ...
    :param luminosity:
    :param unit:
    :param distance:
    :return:
    """

    return estimate_tir_to_sfr("Kennicutt", kennicutt_tir, "erg/s", luminosity, unit=unit, distance=distance)

# -----------------------------------------------------------------
