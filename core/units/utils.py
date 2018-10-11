#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.units.utils Contains unit utilitiy functions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict, OrderedDict

# Import astronomical modules
from astropy.table import Table
from astropy import constants
from astropy.units import Unit, CompositeUnit

# -----------------------------------------------------------------

_c = 2.99792458e8     # light speed in m/s
_AU = 1.49597871e11   # astronomical unit in m
_pc = 3.08567758e16   # parsec in m
_Lsun = 3.839e26      # solar bolometric luminosity in W (without solar neutrino radiation)
_Msun = 1.9891e30     # solar mass in kg
_arcsec2 = 2.350443053909789e-11  # solid angle of 1 square arc second in steradian

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# Flux zero point for AB magnitudes
ab_mag_zero_point = 3631. * Unit("Jy")

# 2MASS F_0 (in Jy)
f_0_2mass = {"2MASS.J": 1594.0, "2MASS.H": 1024.0, "2MASS.Ks": 666.7}

# Extended source photometrical correction coefficients
photometrical_correction_irac = {"IRAC.I1": 0.91, "IRAC.I2": 0.94, "IRAC.I3": 0.66, "IRAC.I4": 0.74}

# Wise magnitude zero points (those in the header are rounded ??)
#m_0_wise = {"WISE.W1": 20.73, "WISE.W2": 19.567, "WISE.W3": 17.600, "WISE.W4": 12.980}

# WISE F_0 (in W / cm2 / um) (from http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/figures/sec4_3gt4.gif)
f_0_wise = {"WISE.W1": 8.1787e-15, "WISE.W2": 2.4150e-15, "WISE.W3": 6.5151e-17, "WISE.W4": 5.0901e-18}

# Absolute calibration magnitudes for WISE bands
absolute_calibration_wise = {"WISE.W1": 0.034, "WISE.W2": 0.041, "WISE.W3": -0.030, "WISE.W4": 0.029}

# Correction factor for the calibration discrepancy between WISE photometric standard blue stars and red galaxies
calibration_discrepancy_wise_w4 = 0.92

# GALEX conversion factors from count/s to flux
#  - FUV: Flux [erg sec-1 cm-2 Angstrom-1] = 1.40 x 10-15 x CPS
#  - NUV: Flux [erg sec-1 cm-2 Angstrom-1] = 2.06 x 10-16 x CPS
galex_conversion_factors = {"GALEX.FUV": 1.40e-15, "GALEX.NUV": 2.06e-16}

# 1 Jy = 1e-23 erg/s/cm/Hz (erg / [s * cm * Hz])
ergscmHz_to_Jy = 1e23

# Effective wavelengths for GALEX
effective_wavelengths = {"GALEX.FUV": 1528.0 * Unit("Angstrom"), "GALEX.NUV": 2271.0 * Unit("Angstrom")}

# Conversion between flux density in SI units and Jansky
jansky_to_si = 1e-26  # 1 Jy in W / [ m2 * Hz]
si_to_jansky = 1e26  # 1 W / [m2 * Hz] in Jy

# -----------------------------------------------------------------

def ab_to_jansky(ab_magnitude):

    """
    This function ...
    :param ab_magnitude:
    :return:
    """

    return ab_mag_zero_point.to("Jy").value * np.power(10.0, -2./5. * ab_magnitude)

# -----------------------------------------------------------------

def jansky_to_ab(jansky_fluxes):

    """
    This function ...
    :param jansky_fluxes:
    :return:
    """

    return -5./2. * np.log10(jansky_fluxes / ab_mag_zero_point.to("Jy").value)

# -----------------------------------------------------------------

# Construct table for conversion between AB and Vega magnitude scale

# Band lambda_eff "mAB - mVega" MSun(AB) MSun(Vega)
# U	3571 0.79 6.35 5.55
# B	4344 -0.09 5.36 5.45
# V	5456 0.02 4.80 4.78
# R	6442 0.21 4.61 4.41
# I	7994 0.45 4.52 4.07
# J	12355 0.91 4.56 3.65
# H	16458 1.39 4.71 3.32
# Ks 21603 1.85	5.14 3.29
# u	3546 0.91 6.38 5.47
# g	4670 -0.08 5.12 5.20
# r	6156 0.16 4.64 4.49
# i	7472 0.37 4.53 4.16
# z	8917 0.54 4.51 3.97

# From: http://cosmo.nyu.edu/mb144/kcorrect/linear.ps
# or http://astroweb.case.edu/ssm/ASTR620/alternateabsmag.html
band_column = ["U", "B", "V", "R", "I", "J", "H", "Ks", "u", "g", "r", "i", "z"]
lambda_column = [3571., 4344., 5456., 6442., 7994., 12355., 16458., 21603., 3546., 4670., 6156., 7472., 8917.]
mab_mvega_column = [0.79, -0.09, 0.02, 0.21, 0.45, 0.91, 1.39, 1.85, 0.91, -0.08, 0.16, 0.37, 0.54]
msun_ab_column = [6.35, 5.36, 4.80, 4.61, 4.52, 4.56, 4.71, 5.14, 6.38, 5.12, 4.64, 4.53, 4.51]
msun_vega_column = [5.55, 5.45, 4.78, 4.41, 4.07, 3.65, 3.32, 3.29, 5.47, 5.20, 4.49, 4.16, 3.97]

data = [band_column, lambda_column, mab_mvega_column, msun_ab_column, msun_vega_column]
names = ["Band", "Effective wavelength", "mAB - mVega", "MSun(AB)", "MSun(Vega"]

# Create the table
#ab_vega_conversion_table = tables.new(data, names)

# Create a new table from the data
ab_vega_conversion_table = Table(data=data, names=names, masked=True)

# -----------------------------------------------------------------

def vega_to_ab(vega_magnitude, band):

    """
    This function ...
    :param vega_magnitude:
    :param band:
    :return:
    """

    from ...core.tools import tables

    # Get the index of the row corresponding to the specified band
    band_index = tables.find_index(ab_vega_conversion_table, band)

    # Return the magnitude in the AB scale
    return vega_magnitude + ab_vega_conversion_table["mAB - mVega"][band_index]

# -----------------------------------------------------------------

def ab_to_vega(ab_magnitude, band):

    """
    This function ...
    :param ab_magnitude:
    :param band:
    :return:
    """

    from ...core.tools import tables

    # Get the index of the row corresponding to the specified band
    band_index = tables.find_index(ab_vega_conversion_table, band)

    # Return the magnitude in the Vega system
    return ab_magnitude - ab_vega_conversion_table["mAB - mVega"][band_index]

# -----------------------------------------------------------------

def photometry_2mass_mag_to_jy(magnitude, band):

    """
    This function ...
    :param magnitude:
    :param band:
    :return:
    """

    return f_0_2mass["2MASS."+band] * 10.**(-magnitude/2.5)

# -----------------------------------------------------------------

# Get magnitudes (asinh magnitudes: Lupton et al. (1999))
# whereby the logarithmic magnitude scale transitions to a linear scale in flux density f at low S/N:
#
# m = −2.5/ln(10) * [ asinh((f/f0)/(2b)) + ln(b) ]
#
# f0 = 3631 Jy, the zero point of the AB flux scale
# The quantity b for the co-addition is given in Table 2, along with the asinh magnitude associated with a zero-flux object.

# From THE SEVENTH DATA RELEASE OF THE SLOAN DIGITAL SKY SURVEY (Kevork N. Abazajian, 2009)
# Table 2
# Asinh Magnitude Softening Parameters for the Co-Addition
# Band    b            Zero-Flux Magnitude    m
#                      (m(f/f0 = 0))          (f/f0 = 10b)
# -----------------------------------------------------------------
# u      1.0 × 10−11   27.50                  24.99
# g      0.43 × 10−11  28.42                  25.91
# r      0.81 × 10−11  27.72                  25.22
# i      1.4 × 10−11   27.13                  24.62
# z      3.7 × 10−11   26.08                  23.57

band_column = ["u", "g", "r", "i", "z"]
b_column = [1.0e-11, 0.43e-11, 0.81e-11, 1.4e-11, 3.7e-11]
zeroflux_column = [27.50, 28.42, 27.72, 27.13, 26.08]
m_column = [24.99, 25.91, 25.22, 24.62, 23.57]

data = [band_column, b_column, zeroflux_column, m_column]
names = ["Band", "b", "Zero-flux magnitude", "m"]
#sdss_asinh_parameters_table = tables.new(data, names)

# -----------------------------------------------------------------

nanomaggy = 3.613e-6 * Unit("Jy")
nanomaggy_string = "(3.613e-6 Jy)"

# -----------------------------------------------------------------

# Input string->unit replacements
input_replacements = OrderedDict()
input_replacements["COUNTS"] = "count"
input_replacements["PIXEL"] = "pix"
input_replacements["DN"] = "count"
input_replacements["SEC"] = "second"
input_replacements["nanomaggy"] = nanomaggy_string
input_replacements["nmaggy"] = nanomaggy_string
input_replacements["nmaggy"] = nanomaggy_string
input_replacements["nmgy"] = nanomaggy_string
input_replacements["nMgy"] = nanomaggy_string
input_replacements["nanomaggies"] = nanomaggy_string
input_replacements["solLum"] = "Lsun" # first convert solLum to Lsun so that UM IS NOT REPLACED BY MICRON, GIVEN RISE TO SOLLMICRON !!! (THIS IS ALSO WHY THE INPUT_REPLACEMENTS IS AN ORDERED DICT)
input_replacements["um"] = "micron"
input_replacements["kiB"] = "kB"
input_replacements["GiB"] = "GB"
input_replacements["MiB"] = "MB"
input_replacements["TiB"] = "TB"

# -----------------------------------------------------------------

dimensionless_unit_strings = ["dimensionless", "fraction", "dimless", "none", "None", "no", "no_unit", "scalar", "n.a."]

# -----------------------------------------------------------------

def get_physical_type(unit_string):

    """
    This function ...
    :param unit_string:
    :return:
    """

    if "[" in unit_string and "]" in unit_string: return unit_string.split("[")[1].split("]")[0].strip()
    else: return None

# -----------------------------------------------------------------

def interpret_physical_type(physical_type):

    """
    This function ...
    :param physical_type:
    :return:
    """

    # Get
    density = physical_type.endswith("density")
    brightness = "surface brightness" in physical_type

    # Return
    return density, brightness

# -----------------------------------------------------------------

def clean_unit_string(string):

    """
    This function ...
    :param string:
    :return:
    """

    # Remove things within brackets
    if "[" in string and "]" in string:
        string = string.split("[")[0]
        # If nothing remains, that is weird
        if string.strip() == "": raise ValueError("Invalid unit")

    # Remove hash words (#density and #brightness)
    from ..tools import strings
    string = strings.remove_hash_words(string)

    for key in input_replacements:
        string = string.replace(key, input_replacements[key])
    if string.count("(") == 1 and string.count(")") == 1 and string.startswith("(") and string.endswith(")"): string = string[1:-1]

    # Return
    if string in dimensionless_unit_strings: return ""
    else: return string

# -----------------------------------------------------------------

def physical_base_types(unit, power=None):

    """
    This function ...
    :param unit:
    :param power:
    :return:
    """

    if power is None: return [base.physical_type for base in unit.bases]
    else: return [base.physical_type for base, pwr in zip(unit.bases, unit.powers) if pwr == power]

# -----------------------------------------------------------------

def physical_base_types_units_and_powers(unit, power=None):

    """
    This function ...
    :param unit:
    :param power:
    :return:
    """

    result = []

    if power is None:
        for base, pwr in zip(unit.bases, unit.powers):
            result.append((base.physical_type, base, pwr))
    else:
        for base, pwr in zip(unit.bases, unit.powers):
            if pwr != power: continue
            result.append((base.physical_type, base, pwr))

    return result

# -----------------------------------------------------------------

def physical_base_types_units_and_powers_as_dict(unit, power=None):

    """
    This function ...
    :param unit:
    :param power:
    :return:
    """

    result = defaultdict(list)

    if power is None:

        for base, pwr in zip(unit.bases, unit.powers):
            result[base.physical_type].append((base, pwr))

    else:

        for base, pwr in zip(unit.bases, unit.powers):
            if pwr != power: continue
            result[base.physical_type].append((base, pwr))

    return result

# -----------------------------------------------------------------

def occurences_for_type(unit, physical_type):

    """
    This function ...
    :param unit:
    :param physical_type:
    :return:
    """

    occurences = []

    for base, power in zip(unit.bases, unit.powers):

        if base.physical_type == physical_type:

            occurences.append((base, power))

    return occurences

# -----------------------------------------------------------------

def reduce_unit(unit):

    """
    This function ...
    :param unit: 
    :return: 
    """

    #print("Reducing unit '" + str(unit) + "' ...")

    base_types = physical_base_types_units_and_powers_as_dict(unit, power=1)

    #print("BASE TYPES", base_types)

    # Flag that is set to True when things as Jy are encountered, which allows us to say that we are not dealing with intrinsic surface brightnesses, but with fluxes
    cannot_be_intrinsic_brightness = False

    # If we have a combination of a flux density and a frequency
    if "spectral flux density" in base_types and "frequency" in base_types:

        #print("1")

        # If they occur more than once, we have something weird
        if len(base_types["spectral flux density"]) == 1 and len(base_types["frequency"]) == 1:

            # Assert power is one
            assert base_types["spectral flux density"][0][1] == 1

            # We have a flux, so no intrinsic brightness
            cannot_be_intrinsic_brightness = True

            # Replace the flux density unit with the decomposed version
            represents = base_types["spectral flux density"][0][0].represents

            #print("REPRESENTS", represents, type(represents), represents.physical_type)
            #print("REPRESENTS x2", represents.represents)

            #flux_density_base_types = physical_base_types_units_and_powers_as_dict(represents)
            #print("REPRESENTS x2", flux_density_base_types)

            # ACCOMODATING THE FACT THAT FOR FOR EXAMPLE MJy, .REPRESENTS RESULTS IN 1e6 Jy, INSTEAD OF IMMEDIATELY
            # 1e6 W / Hz / m2 !!!!
            if represents.physical_type == "spectral flux density":

                factor = represents.scale
                bases = represents.bases
                powers = represents.powers
                #represents = represents.represents

                #print(factor)
                #print(bases)
                #print(powers)

                # If the number of bases is NOT one, we have the 'normal' case, where Jy is decomposed into W, Hz and m
                # and factor, bases and powers will be 1e-26, [Unit("W"), Unit("Hz"), Unit("m")], and [1, -1, -2]
                if len(bases) == 1:

                    # This is almost impossible, but check anyway
                    if powers[0] != 1: raise RuntimeError("Something went wrong")

                    # Decompose again
                    base = bases[0]
                    base_represents = base.represents

                    # Set represents
                    #represents = base_represents
                    #represents._scale *= factor

                    represents = CompositeUnit(base_represents.scale * factor, base_represents.bases, base_represents.powers)

                # Not one: do nothing
                else: pass

            #print("REPRESENTS", represents)

            # Replace
            unit /= base_types["spectral flux density"][0][0]
            unit *= represents

        # Invalid
        elif len(base_types["spectral flux density"]) > 1: raise ValueError("Invalid photometric unit: multiple occurences of spectral flux base unit")
        elif len(base_types["frequency"]) > 1: raise ValueError("Invalid photometric unit: multiple occurences of frequency base unit")

    # If we have a combination of a flux density and a length to the power of 2
    if "spectral flux density" in base_types and "length" in physical_base_types(unit, power=2):

        #print("2")

        # if they occur more than once, we have something weird
        if len(base_types["spectral flux density"]) == 1:

            # Assert power is one
            assert base_types["spectral flux density"][0][1] == 1

            # Can still be intrinsic brightness if distance dependence is done away with and another is added (e.g. * m2 / pc2)

            # Replace the flux density unit with the decomposed version
            represents = base_types["spectral flux density"][0][0].represents

            # Replace
            unit /= base_types["spectral flux density"][0][0]
            unit *= represents

        # Invalid
        else: raise ValueError("Invalid photometric unit: multiple occurences of spectral flux base unit")

    # Return the unit
    return unit, cannot_be_intrinsic_brightness

# -----------------------------------------------------------------

def analyse_unit(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    #print("Analysing unit '" + str(unit) + "' ...")

    # Initialize different parts of the unit
    #scale_factor = unit.scale
    base_unit = Unit("")
    wavelength_unit = Unit("")
    frequency_unit = Unit("")
    length_unit = Unit("")
    solid_angle_unit = Unit("")

    # Special extra unit that is certainly about the physical scale (intrinsic brightness)
    intrinsic_scale_unit = None

    #print("BEFORE", unit)
    #print(unit)

    # Reduce (for example, convert 'Jy * Hz' to '1e-26 W / m2')
    unit, cannot_be_intrinsic_brightness = reduce_unit(unit)
    scale_factor = unit.scale

    #print(unit)
    #print("AFTER", unit)

    # Look in the bases whether different physical types occur twice
    for physical_type in ["frequency", "time", "length", "solid angle"]:

        occurences = occurences_for_type(unit, physical_type)

        #print(occurences)

        # Powers are
        if len(occurences) == 2:

            base_a = occurences[0][0]
            base_b = occurences[1][0]

            power_a = occurences[0][1]
            power_b = occurences[1][1]

            if power_a == -power_b:

                # Eliminate this physical type from the unit

                # Calculate the ratio of the two units
                ratio = (base_a / base_b).to("")

                # To the power
                factor = ratio**power_a

                # Eliminate first base unit
                unit *= base_a**(-power_a)
                unit *= base_b**(-power_b)

                # Correct with factor
                unit = CompositeUnit(factor * unit.scale, unit.bases, unit.powers)

    #print(unit)
    #print(unit.bases)

    # Loop over the bases
    for base, power in zip(unit.bases, unit.powers):

        #print(base, power)

        if power > 0:

            # We have spectral flux density (e.g. Jy) as a base unit
            if base.physical_type == "spectral flux density":

                # Decompose the unit
                sf, bu, wu, fu, du, su, _, _ = analyse_unit(base.represents)

                scale_factor *= sf
                base_unit *= bu

                # Wavelength unit
                if wu != "":

                    if fu != "": raise ValueError("Wavelength unit and frequency unit cannot both be defined")
                    if wavelength_unit != "": raise ValueError("Encountered two wavelength units: " + str(wavelength_unit) + " and " + str(wu))
                    wavelength_unit = wu

                # Frequency unit
                if fu != "":

                    if wu != "": raise ValueError("Wavelength unit and frequency unit cannot both be defined")
                    if frequency_unit != "": raise ValueError("Encountered two frequency units: " + str(frequency_unit) + " and " + str(fu))

                    frequency_unit = fu

                # Distance unit
                if du != "":

                    if length_unit != "": raise ValueError("Encountered two length units: " + str(length_unit) + " and " + str(du))
                    length_unit = du

                    # We cannot have an intrinsic brightness (in other words, the length unit IS A DISTANCE LENGTH UNIT, NOT AN INTRINSIC LENGTH SCALE)
                    cannot_be_intrinsic_brightness = True

                # Solid angle unit
                if su != "":

                    if solid_angle_unit != "": raise ValueError("Encountered two solid angle units: " + str(solid_angle_unit) + " and " + str(su))
                    solid_angle_unit = su

            # Unit of power
            elif base.physical_type == "power":

                if power != 1: raise ValueError("Found a power of " + str(power) + " for a unit of radiative power")
                if base_unit != "": raise ValueError("Found a unit of power but base unit already defined by '" + str(base_unit) + "'")
                base_unit = base

            # Unit of energy
            elif base.physical_type == "energy":

                if power != 1: raise ValueError("Found a power of " + str(power) + " for a unit of energy")
                if base_unit != "": raise ValueError("Found a unit of energy but base unit already defined by '" + str(base_unit) + "'")
                base_unit = base

            # Unknown
            elif base.physical_type == "unknown": base_unit *= base ** power

            # Else
            else: raise ValueError("Not a photometric unit: found " + base.physical_type + "^" + str(power) + " dimension as a base")

        # Time unit
        elif base.physical_type == "time":

            # Check the power
            if power != -1: raise ValueError("Found a unit of time but not as inversely proportional to the base unit (instead the power is " + str(power) + ")")

            # Set
            base_unit *= base ** power

        # Length unit
        elif base.physical_type == "length":

            # No power
            if power == -1:

                # Check whether not already defined
                if wavelength_unit != "": raise ValueError("Encountered two wavelength units: " + str(wavelength_unit) + " and " + str(base))

                # Set
                wavelength_unit = base

            # Power of 2
            elif power == -2:

                # Check whether already defined
                if length_unit != "":

                    # Check whether intrinsic scale unit is not yet defined
                    if intrinsic_scale_unit is not None: raise ValueError("Encountered two intrinsic scale units: " + str(intrinsic_scale_unit) + " and " + str(base**2))

                    # Check which one is the bigger unit (e.g. 'pc' > 'm')
                    # PREVIOUS LENGTH UNIT WAS BIGGER THAN NEW LENGTH UNIT
                    if length_unit > base**2:

                        # Set long length unit as intrinsic scale unit (e.g. 'pc'), and the other (e.g. 'm' as unit for distance (flux = L / distance**2)
                        intrinsic_scale_unit = length_unit
                        length_unit = base**2

                    # PREVIOUS LENGTH UNIT WAS SMALLER THAN NEW LENGTH UNIT
                    else: intrinsic_scale_unit = base**2

                # So far only length unit
                else:

                    # Set
                    length_unit = base ** 2

            # Power of 3
            elif power == -3:

                # Check whether wavelength unit not already defined
                if wavelength_unit != "": raise ValueError("Encountered two wavelength units: " + str(wavelength_unit) + " and " + str(base))

                # Split into wavelength and length unit
                wavelength_unit = base
                length_unit = base ** 2

            # Power of 4
            elif power == -4:

                # Check whether not yet defined
                if length_unit != "": raise ValueError("Encountered two length units: " + str(length_unit) + " and " + str(base**2))
                if intrinsic_scale_unit is not None: raise ValueError("Encountered two intrinsic scale units: " + str(intrinsic_scale_unit) + " and " + str(base**2))

                # Split into length (for distance) and intrinsic scale unit
                length_unit = base ** 2
                intrinsic_scale_unit = base ** 2

            # Invalid power
            else: raise ValueError("Not a photometric unit: found length^" + str(power) + " dimension")

        # Frequency unit
        elif base.physical_type == "frequency":

            # Check whether not already defined
            if frequency_unit != "": raise ValueError("Encountered two frequency units: " + str(frequency_unit) + " and " + str(base))

            # Check power
            if power != -1: raise ValueError("Found a unit of frequency but not as inversely proportional to the base unit (instead the power is " + str(power) + ")")

            # Set
            frequency_unit = base

        # Solid angle unit
        elif base.physical_type == "solid angle":

            # Check whether not already defined
            if solid_angle_unit != "": raise ValueError("Encountered two solid angle units: " + str(solid_angle_unit) + " and " + str(base))

            # Check power
            if power != -1: raise ValueError("Found a unit of solid angle but not as inversely proportional to the base unit (instead the power is " + str(power) + ")")

            # Set
            solid_angle_unit = base

        # Angle unit
        elif base.physical_type == "angle":

            # Check whether not already defined
            if solid_angle_unit != "": raise ValueError("Found an angle unit but the solid angle unit is already defined: " + str(solid_angle_unit))

            # Check power
            if power != -2: raise ValueError("Found an angle unit but not as squared inversely proportional to the base unit (instead the power is " + str(power) + ")")

            # Set
            solid_angle_unit = base ** 2

        # Not recognized
        else: raise ValueError("Not a photometric unit: found " + str(base) + "^" + str(power))

    # Check if wavelength and frequency unit are not both defined
    if wavelength_unit != "" and frequency_unit != "": raise ValueError("Not a photometric unit: found wavelength^-1 and frequency^-1 dimensions")

    # Check whether a base unit is found
    if base_unit is None or base_unit == "": raise ValueError("Not a photometric unit: found no unit of energy or luminosity")

    # Return
    return scale_factor, base_unit, wavelength_unit, frequency_unit, length_unit, solid_angle_unit, cannot_be_intrinsic_brightness, intrinsic_scale_unit

# -----------------------------------------------------------------

def divide_units_reverse(unit_a, unit_b):

    """
    This function ...
    :param unit_a:
    :param unit_b:
    :return:
    """

    # Re-evaluate everything, cannot be a photometric quantityb anymore
    return CompositeUnit(1, [unit_b, unit_a], [1, -1], _error_check=False)

# -----------------------------------------------------------------
