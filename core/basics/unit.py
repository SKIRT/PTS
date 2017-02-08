#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.unit Contains the PhotometricUnit class.

# -----------------------------------------------------------------

# Import standard modules
import math
import copy
import warnings
import numpy as np

# Import astronomical modules
from astropy.units import Unit, UnitBase, CompositeUnit, spectral, Quantity
from astropy import constants
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...magic.basics.pixelscale import Pixelscale
from .quantity import PhotometricQuantity, parse_quantity

# -----------------------------------------------------------------

_c = 2.99792458e8     # light speed in m/s
_AU = 1.49597871e11   # astronomical unit in m
_pc = 3.08567758e16   # parsec in m
_Lsun = 3.839e26      # solar bolometric luminosity in W (without solar neutrino radiation)
_Msun = 1.9891e30     # solar mass in kg
_arcsec2 = 2.350443053909789e-11  # solid angle of 1 square arc second in steradian

# -----------------------------------------------------------------

nanomaggy = 3.613e-6 * Unit("Jy")
nanomaggy_string = "(3.613e-6 Jy)"

# -----------------------------------------------------------------

# Input string->unit replacements
input_replacements = dict()
input_replacements["DN"] = "count"
input_replacements["SEC"] = "second"
input_replacements["nanomaggy"] = nanomaggy_string
input_replacements["nmaggy"] = nanomaggy_string
input_replacements["nmaggy"] = nanomaggy_string
input_replacements["nmgy"] = nanomaggy_string
input_replacements["nMgy"] = nanomaggy_string
input_replacements["nanomaggies"] = nanomaggy_string
input_replacements["um"] = "micron"

# -----------------------------------------------------------------

# Output unit->string replacements
output_replacements = dict()
output_replacements["solMass"] = "Msun"
output_replacements["solLum"] = "Lsun"

# -----------------------------------------------------------------

# LUMINOSITY: W
# WAVELENGTH LUMINOSITY DENSITY: W/micron
# FREQUENCY LUMINOSITY DENSITY: W/Hz
# NEUTRAL LUMINOSITY DENSITY: W

# FLUX: W/m2
# WAVELENGTH FLUX DENSITY: W/m2/micron
# FREQUENCY FLUX DENSITY: W/m2/Hz
# NEUTRAL FLUX DENSITY: W/m2

# INTENSITY: W/sr
# WAVELENGTH INTENSITY DENSITY: W/sr/micron
# FREQUENCY INTENSITY DENSITY: W/sr/Hz
# NEUTRAL INTENSITY DENSITY: W/m2

# SURFACE BRIGHTNESS: W/m2/sr
# WAVELENGTH SURFACE BRIGHTNESS DENSITY: W/m2/sr/micron
# FREQUENCY SURFACE BRIGHTNESS DENSITY: W/m2/sr/micron
# NEUTRAL SURFACE BRIGHTNESS DENSITY: W/m2/sr

# -----------------------------------------------------------------

# SKIRT:
# 'W/m': 'wavelengthluminositydensity',
#'W/micron': 'wavelengthluminositydensity',
#'Lsun/micron': 'wavelengthluminositydensity',
#'W/Hz': 'frequencyluminositydensity',
#'erg/s/Hz': 'frequencyluminositydensity',
#'Lsun/Hz': 'frequencyluminositydensity',
#'W/m2': 'neutralfluxdensity',
#'W/m2/sr': 'neutralsurfacebrightness',
#'W/m2/arcsec2': 'neutralsurfacebrightness',
#'W/m3': 'wavelengthfluxdensity',
#'W/m2/micron': 'wavelengthfluxdensity',
#'W/m3/sr': 'wavelengthsurfacebrightness',
#'W/m2/micron/sr': 'wavelengthsurfacebrightness',
#'W/m2/micron/arcsec2': 'wavelengthsurfacebrightness',
#'W/m2/Hz': 'frequencyfluxdensity',
#'Jy': 'frequencyfluxdensity',
#'mJy': 'frequencyfluxdensity',
#'MJy': 'frequencyfluxdensity',
#'erg/s/cm2/Hz': 'frequencyfluxdensity',
#'W/m2/Hz/sr': 'frequencysurfacebrightness',
#'W/m2/Hz/arcsec2': 'frequencysurfacebrightness',
#'Jy/sr': 'frequencysurfacebrightness',
#'Jy/arcsec2': 'frequencysurfacebrightness',
#'MJy/sr': 'frequencysurfacebrightness',
#'MJy/arcsec2': 'frequencysurfacebrightness'

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

def clean_unit_string(string):

    """
    This function ...
    :param string:
    :return:
    """

    for key in input_replacements:
        string = string.replace(key, input_replacements[key])
    if string.count("(") == 1 and string.count(")") == 1 and string.startswith("(") and string.endswith(")"): string = string[1:-1]
    return string

# -----------------------------------------------------------------

def parse_unit(argument, density=False):

    """
    This function ...
    :param argument:
    :param density:
    :return:
    """

    try: unit = PhotometricUnit(argument, density=density)
    except ValueError:
        if isinstance(argument, basestring): argument = clean_unit_string(argument)
        unit = Unit(argument)
    return unit

# -----------------------------------------------------------------

def stringify_unit(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    # Get parsing type
    if isinstance(unit, PhotometricUnit):
        if unit.density: parsing_type = "photometric_density_unit"
        else: parsing_type = "photometric_unit"
    else: parsing_type = "unit"

    # Return type and stringified unit
    return parsing_type, represent_unit(unit)

# -----------------------------------------------------------------

def represent_unit(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    string = str(unit)

    for key in output_replacements:
        string = string.replace(key, output_replacements[key])

    return string

# -----------------------------------------------------------------

def analyse_unit(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    scale_factor = unit.scale
    base_unit = Unit("")
    wavelength_unit = Unit("")
    frequency_unit = Unit("")
    distance_unit = Unit("")
    solid_angle_unit = Unit("")

    for base, power in zip(unit.bases, unit.powers):

        if power > 0:

            if base.physical_type == "spectral flux density":

                sf, bu, wu, fu, du, su = analyse_unit(base.represents)

                scale_factor *= sf
                base_unit *= bu

                if wu != "":

                    assert fu == ""
                    assert wavelength_unit == ""

                    wavelength_unit = wu

                if fu != "":

                    assert wu == ""
                    assert frequency_unit == ""

                    frequency_unit = fu

                if du != "":

                    assert distance_unit == ""
                    distance_unit = du

                if su != "":

                    assert su == ""
                    solid_angle_unit = su

            elif base.physical_type == "power":
                assert power == 1
                assert base_unit == ""
                base_unit = base

            elif base.physical_type == "energy":
                assert power == 1
                assert base_unit == ""
                base_unit = base

            elif base.physical_type == "unknown":
                base_unit *= base ** power

            else: raise ValueError("Not a photometric unit: found " + base.physical_type + "^" + str(power) + " dimension as a base")

        elif base.physical_type == "time":
            assert power == -1
            base_unit *= base ** power
        elif base.physical_type == "length":
            if power == -1:
                wavelength_unit = base
            elif power == -2:
                distance_unit = base ** 2
            elif power == -3:
                wavelength_unit = base
                distance_unit = base ** 2
            else: raise ValueError("Not a photometric unit: found length^" + str(power) + " dimension")
        elif base.physical_type == "frequency":
            assert power == -1, power
            frequency_unit = base
        elif base.physical_type == "solid angle":
            assert power == -1, power
            solid_angle_unit = base
        elif base.physical_type == "angle":
            assert solid_angle_unit == ""
            assert power == -2, power
            solid_angle_unit = base ** power

    # Check if wavelength and frequency unit are not both defined
    if wavelength_unit != "" and frequency_unit != "": raise ValueError("Not a photometric unit: found wavelength^-1 and frequency^-1 dimensions")

    if base_unit is None or base_unit == "": raise ValueError("Not a photometric unit: found no unit of energy or luminosity")

    # Return
    return scale_factor, base_unit, wavelength_unit, frequency_unit, distance_unit, solid_angle_unit

# -----------------------------------------------------------------

class PhotometricUnit(CompositeUnit):

    """
    This function ...
    """

    __slots__ = [] # make the class objects immutable

    # -----------------------------------------------------------------

    def __init__(self, unit, density=False):

        """
        The constructor ...
        :param unit:
        :param density:
        """

        # Already a photometric unit
        if isinstance(unit, PhotometricUnit):

            self.density = unit.density
            self.scale_factor = unit.scale_factor
            self.base_unit = unit.base_unit
            self.wavelength_unit = unit.wavelength_unit
            self.frequency_unit = unit.frequency_unit
            self.distance_unit = unit.distance_unit
            self.solid_angle_unit = unit.solid_angle_unit

        # Regular unit
        else:

            # Clean unit string
            if isinstance(unit, basestring): unit = clean_unit_string(unit)

            # Parse the unit
            try: unit = Unit(unit)
            except ValueError: raise ValueError("Unit is not recognized")

            # Remove 'per pixel' from the unit
            if "pix" in str(unit): unit *= "pix"

            # Set whether it represents a density
            self.density = density

            # Analyse the unit
            self.scale_factor, self.base_unit, self.wavelength_unit, self.frequency_unit, self.distance_unit, self.solid_angle_unit = analyse_unit(unit)

            # If the wavelength unit is not None or the frequency unit is not None, we have a spectral density
            if self.wavelength_unit is not None and self.wavelength_unit != "": self.density = True
            if self.frequency_unit is not None and self.frequency_unit != "": self.density = True

        # Call the constructor of the base class
        super(PhotometricUnit, self).__init__(unit.scale, unit.bases, unit.powers)

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def __pow__(self, p):

        """
        This function ...
        :param p:
        :return:
        """

        # If the power is one, return a copy of this photometric unit
        if p == 1: return self.copy()

        # Not a photometric unit anymore
        else: return CompositeUnit(1, [self], [p])

    # -----------------------------------------------------------------

    def __div__(self, other):

        """
        This function ...
        :param m:
        :return:
        """

        # If other is a string
        if isinstance(other, basestring): other = parse_unit(other)

        # Divided by another unit
        if isinstance(other, UnitBase): return divide_units(self, other)

        # Divided by a quantity
        elif hasattr(other, "unit"):

            # Get the new unit
            new_unit = divide_units(self, other.unit)

            # Create a quantity
            quantity = parse_quantity(other.value * new_unit)

            # Return the quantity
            return quantity

        # Divided by a number
        else: return PhotometricQuantity(1./other, self)

    # -----------------------------------------------------------------

    def __rdiv__(self, other):

        """
        This function ...
        :param m:
        :return:
        """

        # If the other is a string
        if isinstance(other, basestring): other = parse_unit(other)

        # Another unit divided by this unit
        if isinstance(other, UnitBase): return divide_units_reverse(self, other)

        # Quantity divided by this unit
        elif hasattr(other, "unit"):

            # Get the new unit
            new_unit = divide_units_reverse(self, other)

            # Create a quantity
            quantity = parse_quantity(other.value * new_unit)

            # Return the quantity
            return quantity

        # Regular number divided by this unit
        else: return parse_quantity(Quantity(other, 1./self))

    # -----------------------------------------------------------------

    __truediv__ = __div__

    # -----------------------------------------------------------------

    __rtruediv__ = __rdiv__

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # If the other is a string
        if isinstance(other, basestring): other = parse_unit(other)

        # This unit is multiplied with another unit
        if isinstance(other, UnitBase): return multiply_units(self, other)

        # This unit is multiplied with a quantity
        elif hasattr(other, "unit"):

            # Get the new unit
            new_unit = multiply_units(self, other.unit)

            # Create a quantity
            quantity = parse_quantity(other.value * new_unit)

            # Return the quantity
            return quantity

        # Regular number multiplied by this unit
        else: return PhotometricQuantity(other, self)

    # -----------------------------------------------------------------

    __rmul__ = __mul__

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Try to parse as a photometric unit
        try: other = PhotometricUnit(other)
        except ValueError: raise ValueError("The other unit is not a photometric unit")

        # Use implementation in base class
        return super(PhotometricUnit, self).__eq__(other)

    # -----------------------------------------------------------------

    @property
    def is_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.base_physical_type == "luminosity"

    # -----------------------------------------------------------------

    @property
    def is_flux(self):

        """
        This function ...
        :return:
        """

        return self.base_physical_type == "flux"

    # -----------------------------------------------------------------

    @property
    def is_intensity(self):

        """
        This function ...
        :return:
        """

        return self.base_physical_type == "intensity"

    # -----------------------------------------------------------------

    @property
    def is_surface_brightness(self):

        """
        This function ...
        :return:
        """

        return self.base_physical_type == "surface brightness"

    # -----------------------------------------------------------------

    @property
    def base_physical_type(self):

        """
        This function ...
        :return:
        """

        if self.base_unit.physical_type == "power":

            if self.distance_unit != "":
                if self.solid_angle_unit != "": base = "surface brightness"
                else: base = "flux"
            else:
                if self.solid_angle_unit != "": base = "intensity"
                else: base = "luminosity"

        elif self.base_unit.physical_type == "frequency":
            base = "detection rate"

        else: base = "detections"

        # Return the base type
        return base

    # -----------------------------------------------------------------

    @property
    def spectral_density_type(self):

        """
        This function ...
        :return:
        """

        # Wavelength densities
        if self.is_wavelength_density: return "wavelength"

        # Frequency densities
        elif self.is_frequency_density: return "frequency"

        # Neutral density or regular power/flux/surfacebrightness
        elif self.is_neutral_density: return "neutral"

        # Not a spectral density
        else: return None

    # -----------------------------------------------------------------

    @property
    def is_spectral_density(self):

        """
        This function ...
        :return:
        """

        return self.density

    # -----------------------------------------------------------------

    @property
    def is_wavelength_density(self):

        """
        This function ...
        :return:
        """

        if self.wavelength_unit != "":
            assert self.frequency_unit == ""
            return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def is_frequency_density(self):

        """
        This function ...
        :return:
        """

        if self.frequency_unit != "":
            assert self.wavelength_unit == ""
            return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def is_neutral_density(self):

        """
        This function ...
        :return:
        """

        if self.wavelength_unit == "" and self.frequency_unit == "": return self.density
        else: return False

    # -----------------------------------------------------------------

    @property
    def physical_type(self):

        """
        This function ...
        :return:
        """

        # Get the spectral density type
        density_type = self.spectral_density_type

        # Return the type
        if density_type is None: return self.base_physical_type
        else: return " ".join([density_type, self.base_physical_type, "density"])

    # -----------------------------------------------------------------

    def conversion_factor(self, to_unit, density=False, wavelength=None, frequency=None, distance=None, solid_angle=None,
                          fltr=None, pixelscale=None):

        """
        This function ...
        :param to_unit:
        :param density:
        :param wavelength:
        :param frequency:
        :param distance:
        :param solid_angle:
        :param fltr:
        :param pixelscale:
        :return:
        """

        # Parse "to unit"
        to_unit = PhotometricUnit(to_unit, density=density)

        # Determine wavelength and frequency
        if wavelength is not None:
            if frequency is not None: raise ValueError("Either frequency or wavelength can be specified")
            frequency = wavelength.to("Hz", equivalencies=spectral())
        elif frequency is not None:
            wavelength = frequency.to("micron", equivalencies=spectral())
        elif fltr is not None:
            wavelength = fltr.pivot
            frequency = frequency.to("Hz", equivalencies=spectral())

        # Same type
        if self.physical_type == to_unit.physical_type:
            #factor = self.scale_factor * self / self.to_unit
            factor = self / to_unit
            return factor.scale

        # Convert
        if isinstance(solid_angle, basestring):
            from ..tools import parsing
            solid_angle = parsing.quantity(solid_angle)

        # Convert
        if isinstance(distance, basestring):
            from ..tools import parsing
            distance = parsing.quantity(distance)

        # Convert
        if isinstance(frequency, basestring):
            from ..tools import parsing
            frequency = parsing.quantity(frequency)

        # Convert
        if isinstance(wavelength, basestring):
            from ..tools import parsing
            wavelength = parsing.quantity(wavelength)

        # Convert
        if isinstance(pixelscale, basestring):
            from ..tools import parsing
            pixelscale = parsing.quantity(pixelscale)

        # Make PixelScale instance
        if isinstance(pixelscale, Quantity): pixelscale = Pixelscale(pixelscale)
        elif isinstance(pixelscale, Pixelscale): pass
        elif pixelscale is None: pass
        else: raise ValueError("Don't know what to do with pixelscale of type " + str(type(pixelscale)))

        # If solid angle is None, convert pixelscale to solid angle (of one pixel)
        if solid_angle is None and pixelscale is not None: solid_angle = pixelscale.solid_angle

        # Neutral density
        if self.is_neutral_density:

            if to_unit.is_wavelength_density: new_unit = self / wavelength
            elif to_unit.is_frequency_density: new_unit = self / frequency
            elif to_unit.is_neutral_density: new_unit = self
            else: raise ValueError("Cannot convert from spectral density to integrated quantity") # asked to convert to not a spectral density

        # Wavelength density
        elif self.is_wavelength_density:

            if to_unit.is_neutral_density: new_unit = self * wavelength
            elif to_unit.is_frequency_density: new_unit = self * wavelength / frequency
            elif to_unit.is_wavelength_density: new_unit = self
            else: raise ValueError("Cannot convert from spectral density to integrated quantity")

        # Frequency density
        elif self.is_frequency_density:

            if to_unit.is_neutral_density: new_unit = self * frequency
            elif to_unit.is_frequency_density: new_unit = self
            elif to_unit.is_wavelength_density: new_unit = self * frequency / wavelength
            else: raise ValueError("Cannot convert from spectral density to integrated quantity")

        # Not a spectral density
        else:

            if to_unit.is_neutral_density: raise ValueError("Cannot convert from integrated quantity to spectral density")
            elif to_unit.is_frequency_density: raise ValueError("Cannot convert from integrated quantity to spectral density")
            elif to_unit.is_wavelength_density: raise ValueError("Cannot convert from integrated quantity to spectral density")
            else: new_unit = self

        # Same base type
        if self.base_physical_type == to_unit.base_physical_type:

            # Determine factor
            factor = new_unit.to(to_unit)
            return factor

        # Different base type, luminosity
        elif self.base_physical_type == "luminosity":

            if to_unit.base_physical_type == "flux":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Divide by 4 pi distance**2
                new_unit /= (4.0 * math.pi * distance**2)

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "intensity":

                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Divide by solid angle
                new_unit /= solid_angle

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "surface brightness":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)
                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Divide by 4 pi distance**2 and solid angle
                new_unit /= (4.0 * math.pi * distance**2 * solid_angle)

                # Determine factor
                factor = new_unit.to(to_unit).value

            # Invalid
            else: raise RuntimeError("We shouldn't reach this part")

            # Return the conversion factor
            return factor

        # Different base type, flux
        elif self.base_physical_type == "flux":

            if to_unit.base_physical_type == "luminosity":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Multiply by 4 pi distance**2
                new_unit *= (4.0 * math.pi * distance ** 2)

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "surface brightness":

                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Divide by solid angle
                new_unit /= solid_angle

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "intensity":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)
                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Divide by solid angle
                # Multiply by 4 pi distance**2
                new_unit /= solid_angle * (4.0 * math.pi * distance ** 2)

                # Determine factor
                factor = new_unit.to(to_unit).value

            # Invalid
            else: raise RuntimeError("We shouldn't reach this part")

            # Return the conversion factor
            return factor

        # Different base type, intensity
        elif self.base_physical_type == "intensity":

            if to_unit.base_physical_type == "luminosity":

                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Multiply by solid angle
                new_unit *= solid_angle

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "flux":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)
                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Multiply by solid angle
                # Divide by 4 pi distance**2
                new_unit *= solid_angle / (4.0 * math.pi * distance ** 2)

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "surface brightness":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Divide by 4 pi distance**2
                new_unit /= (4.0 * math.pi * distance ** 2)

                # Determine factor
                factor = new_unit.to(to_unit).value

            # Invalid
            else: raise RuntimeError("We shouldn't reach this part")

            # Return the conversion factor
            return factor

        # Different base type, surface brightness
        else:

            if to_unit.base_physical_type == "luminosity":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)
                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Multiply by 4 pi distance**2 and solid angle
                new_unit *= (4.0 * math.pi * distance ** 2 * solid_angle)

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "flux":

                if solid_angle is None: raise ValueError("Solid angle should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Multiply by solid angle
                new_unit *= solid_angle

                # Determine factor
                factor = new_unit.to(to_unit).value

            elif to_unit.base_physical_type == "intensity":

                if distance is None: raise ValueError("Distance should be specified for conversion from " + self.physical_type + " to " + to_unit.physical_type)

                # Multiply by 4 pi distance**2
                new_unit *= (4.0 * math.pi * distance ** 2)

                # Determine factor
                factor = new_unit.to(to_unit).value

            # Invalid
            else: raise RuntimeError("We shouldn't reach this part")

            # Return the conversion factor
            return factor

# -----------------------------------------------------------------

def is_wavelength(unit):

    """
    This fucntion ...
    :return:
    """

    return unit.physical_type == "length"

# -----------------------------------------------------------------

def is_frequency(unit):

    """
    This function ...
    :return:
    """

    return unit.physical_type == "frequency"

# -----------------------------------------------------------------

def is_time(unit):

    """
    This function ...
    :return:
    """

    return unit.physical_type == "time"

# -----------------------------------------------------------------

def is_inverse_wavelength(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    if len(unit.bases) != 1: return False
    return unit.bases[0].physical_type == "length" and unit.powers[0] == -1

# -----------------------------------------------------------------

def is_inverse_frequency(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    if len(unit.bases) != 1: return False
    return (unit.bases[0].physical_type == "time" and unit.powers[0] == 1) or (unit.bases[0].physical_type == "frequency" and unit.powers[0] == -1)

# -----------------------------------------------------------------

def contains_wavelength(unit):

    """
    This function ...
    :return:
    """

    # Loop over the bases
    for base, power in zip(unit.bases, unit.powers):
        if power == 1 and is_wavelength(base): return True

    return False

# -----------------------------------------------------------------

def contains_frequency(unit):

    """
    This function ...
    :return:
    """

    # Loop over the bases
    for base, power in zip(unit.bases, unit.powers):
        if power == 1 and is_frequency(base): return True
        if power == -1 and is_time(base): return True

    return False

# -----------------------------------------------------------------

def contains_inverse_frequency(unit):

    """
    This function ...
    """

    # Loop over the bases
    for base, power in zip(unit.bases, unit.powers):
        if power == 1 and is_time(base): return True
        if power == -1 and is_frequency(base): return True
    return False

# -----------------------------------------------------------------

def contains_inverse_wavelength(unit):

    """
    This function ...
    :param unit:
    :return:
    """

    # Loop over the bases
    for base, power in zip(unit.bases, unit.powers):
        if power == -1 and is_wavelength(base): return True
    return False

# -----------------------------------------------------------------

def make_composite_multiplication(unit_a, unit_b):

    """
    This function ...
    :param unit_a:
    :param unit_b:
    :return:
    """

    return CompositeUnit(1, [unit_a, unit_b], [1, 1])

# -----------------------------------------------------------------

def multiply_units(unit_a, unit_b):

    """
    This function ...
    :param unit_a:
    :param unit_b:
    :return:
    """

    # If the other unit is dimensionless
    if unit_b == "": return unit_a.copy()

    # If the other unit is dimensionless with a certain scale
    elif unit_b.physical_type == "dimensionless" and unit_b.scale != 1: return PhotometricUnit(CompositeUnit(unit_a.scale * unit_b.scale, unit_a.bases, unit_a.powers), density=unit_a.is_spectral_density)

    # Spectral density
    if unit_a.is_spectral_density:

        # If this is a wavelength density
        if unit_a.is_wavelength_density:

            # From wavelength spectral density to neutral spectral density
            if is_wavelength(unit_b): return PhotometricUnit(make_composite_multiplication(unit_a, unit_b), density=True)
            elif contains_wavelength(unit_b): return parse_unit(make_composite_multiplication(unit_a, unit_b), density=True)
            else: return parse_unit(make_composite_multiplication(unit_a, unit_b))

        # If this is a frequency density
        elif unit_a.is_frequency_density:

            # From frequency spectral density to neutral spectral density
            if is_frequency(unit_b): return PhotometricUnit(make_composite_multiplication(unit_a, unit_b), density=True)
            elif contains_wavelength(unit_b): return parse_unit(make_composite_multiplication(unit_a, unit_b), density=True)
            else: return parse_unit(make_composite_multiplication(unit_a, unit_b))

        # Neutral density
        else:

            # From netural
            if is_inverse_wavelength(unit_b): return PhotometricUnit(make_composite_multiplication(unit_a, unit_b), density=True)
            elif is_inverse_frequency(unit_b): return PhotometricUnit(make_composite_multiplication(unit_a, unit_b), density=True)
            elif contains_inverse_wavelength(unit_b): return parse_unit(make_composite_multiplication(unit_a, unit_b), density=True)
            elif contains_inverse_frequency(unit_b): return parse_unit(make_composite_multiplication(unit_a, unit_b), density=True)
            else: return parse_unit(make_composite_multiplication(unit_a, unit_b))

    # Not a spectral density
    else:

        # If unit b is an inverse wavelength
        if is_inverse_wavelength(unit_b):
            unit = PhotometricUnit(make_composite_multiplication(unit_a, unit_b), density=False)
            if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by multiplication with the unit '" + str(unit_b) + ". This may not be the intention.")
            return unit

        # If unit b is an inverse frequency
        elif is_inverse_frequency(unit_b):
            unit = PhotometricUnit(make_composite_multiplication(unit_a, unit_b), density=False)
            if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by multiplication with the unit '" + str(unit_b) + ". This may not be the intention.")
            return unit

        # Try parsing as spectral photometric quantity (density=True), but possibly no photometric quantity
        elif contains_inverse_wavelength(unit_b):
            unit = parse_unit(make_composite_multiplication(unit_a, unit_b), density=False)
            try:
                if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by multiplication with the unit '" + str(unit_b) + ". This may not be the intention.")
            except AttributeError: pass
            return unit

        # Try parsing as spectral photometric quantity(density=True), but possibly no photometric quantity
        elif contains_inverse_frequency(unit_b):
            unit = parse_unit(make_composite_multiplication(unit_a, unit_b), density=False)
            try:
                if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by multiplication with the unit '" + str(unit_b) + ". This may not be the intention.")
            except AttributeError: pass
            return unit

        # Parse regularly
        else: return parse_unit(make_composite_multiplication(unit_a, unit_b))

# -----------------------------------------------------------------

def make_composite_division(unit_a, unit_b):

    """
    This function ...
    :param unit_a:
    :param unit_b:
    :return:
    """

    return CompositeUnit(1, [unit_a, unit_b], [1, -1], _error_check=False)

# -----------------------------------------------------------------

def divide_units(unit_a, unit_b):

    """
    This function ...
    :param unit_a:
    :param unit_b:
    :return:
    """

    # If the other unit is dimensionless
    if unit_b == "": return unit_a.copy()

    # If the other unit is dimensionless with a certain scale
    elif unit_b.physical_type == "dimensionless" and unit_b.scale != 1: return PhotometricUnit(CompositeUnit(unit_a.scale / unit_b.scale, unit_a.bases, unit_a.powers), density=unit_a.is_spectral_density)

    # If we have a spectral density
    if unit_a.is_spectral_density:

        # If this is a wavelength density
        if unit_a.is_wavelength_density:

            if is_inverse_wavelength(unit_b): return PhotometricUnit(make_composite_division(unit_a, unit_b), density=True)
            elif contains_inverse_wavelength(unit_b): return parse_unit(make_composite_division(unit_a, unit_b), density=True)
            else: return parse_unit(make_composite_division(unit_a, unit_b))

        # If this is a frequency density
        if unit_a.is_frequency_density:

            if is_inverse_frequency(unit_b): return PhotometricUnit(make_composite_division(unit_a, unit_b), density=True)
            elif contains_inverse_frequency(unit_b): return parse_unit(make_composite_division(unit_a, unit_b), density=True)
            else: return parse_unit(make_composite_division(unit_a, unit_b))

        # If this is a neutral density
        else:

            if is_wavelength(unit_b): return PhotometricUnit(make_composite_division(unit_a, unit_b), density=True)
            elif is_frequency(unit_b): return PhotometricUnit(make_composite_division(unit_a, unit_b), density=True)
            elif contains_wavelength(unit_b): return parse_unit(make_composite_division(unit_a, unit_b), density=True)
            elif contains_frequency(unit_b): return parse_unit(make_composite_division(unit_a, unit_b), density=True)
            else: return parse_unit(make_composite_division(unit_a, unit_b))

    # Not a spectral density
    else:

        # If unit b is a wavelength
        if is_wavelength(unit_b):
            unit = PhotometricUnit(make_composite_division(unit_a, unit_b), density=False)
            if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by division with the unit '" + str(unit_b) + "'. This may not be the intention.")
            return unit

        # Unit b is a frequency
        elif is_frequency(unit_b):
            unit = PhotometricUnit(make_composite_division(unit_a, unit_b), density=False)
            if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by division with the unit '" + str(unit_b) + "'. This may not be the intention.")

        # Unit b contains a wavelength
        elif contains_wavelength(unit_b):
            unit = parse_unit(make_composite_division(unit_a, unit_b), density=False)
            try:
                if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by division with unit '" + str(unit_b) + "'. This may not be the intention.")
            except AttributeError: pass
            return unit

        elif contains_frequency(unit_b):
            unit = parse_unit(make_composite_division(unit_a, unit_b), density=False)
            try:
                if unit.density: warnings.warn("A " + unit_a.physical_type + " unit is converted to a " + unit.physical_type + " by division with unit '" + str(unit_b) + "' . This may not be the intention.")
            except AttributeError: pass
            return unit

        # Parse regularly
        else: return parse_unit(make_composite_division(unit_a, unit_b))

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
