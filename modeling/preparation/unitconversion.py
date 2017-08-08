#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.unitconversion Contains the UnitConverter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy import constants
from astropy.units import spectral

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import tables
from ...core.basics.log import log
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# Flux zero point for AB magnitudes
ab_mag_zero_point = 3631. * u("Jy")

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
effective_wavelengths = {"GALEX.FUV": 1528.0 * u("Angstrom"), "GALEX.NUV": 2271.0 * u("Angstrom")}

# Conversion between flux density in SI units and Jansky
jansky_to_si = 1e-26  # 1 Jy in W / [ m2 * Hz]
si_to_jansky = 1e26  # 1 W / [m2 * Hz] in Jy

# -----------------------------------------------------------------

def spectral_factor(wavelength, wavelength_unit, frequency_unit):

    """
    This function ...
    :param wavelength:
    :param wavelength_unit:
    :param frequency_unit:
    :return:
    """

    # Convert string units to Unit objects
    if isinstance(wavelength_unit, basestring): wavelength_unit = u(wavelength_unit)
    if isinstance(frequency_unit, basestring): frequency_unit = u(frequency_unit)

    conversion_factor_unit = wavelength_unit / frequency_unit

    # Calculate the conversion factor
    return (wavelength**2 / speed_of_light).to(conversion_factor_unit).value

# -----------------------------------------------------------------

def neutral_fluxdensity_to_jansky(flux_density, wavelength):

    """
    This function converts a flux density in W/m2 (lambda F_lambda = nu F_nu) to Jansky
    :param flux_density:
    :param wavelength:
    :return:
    """

    # Conversion from W / m2 to W / m3
    wavelength_fluxdensity = flux_density / wavelength.to("m").value

    # Conversion from W / m3 to Jy
    return wavelength_fluxdensity_to_jansky(wavelength_fluxdensity, wavelength)

# -----------------------------------------------------------------

def wavelength_fluxdensity_to_jansky(flux_density, wavelength):

    """
    This function converts a flux density in W/m3 (F_lambda) to Jansky
    :param flux_density:
    :param wavelength:
    :return:
    """

    # Conversion from W / m3 to W / [ m2 * Hz]
    frequency_fluxdensity = flux_density * spectral_factor(wavelength, "m", "Hz")

    # Conversion from W / [m2 * Hz] to Jy
    return frequency_fluxdensity_to_jansky(frequency_fluxdensity)

# -----------------------------------------------------------------

def frequency_fluxdensity_to_jansky(flux_density):

    """
    This function converts a flux density in W/m2/Hz (F_nu) to Jansky
    :param flux_density:
    :return:
    """

    # Conversion from W / [m2 * Hz] (per pixel2) to Jy per pixel
    return flux_density * si_to_jansky

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
ab_vega_conversion_table = tables.new(data, names)

# -----------------------------------------------------------------

def vega_to_ab(vega_magnitude, band):

    """
    This function ...
    :param vega_magnitude:
    :param band:
    :return:
    """

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

class UnitConverter(Configurable):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(UnitConverter, self).__init__(*args, **kwargs)

        # The image
        self.image = None

        # The target unit
        self.target_unit = None

        # The conversion factor
        self.conversion_factor = 1.0

    # -----------------------------------------------------------------

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup(image)

        # 2. Determine the conversion factor
        self.convert()

        # 3. Apply the conversion to the image
        self.apply()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set default values for attributes
        self.image = None
        self.target_unit = None
        self.conversion_factor = 1.0

    # -----------------------------------------------------------------

    def setup(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # Call the setup function of the base class
        super(UnitConverter, self).setup()

        # Create a unit object
        self.target_unit = u(self.config.to_unit)

        # Set the image reference
        self.image = image

    # -----------------------------------------------------------------

    def convert(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the conversion factor to set image in " + str(self.target_unit) + " units ...")

        # Skip the unit conversion for images that are already in the right unit
        if self.image.unit == self.target_unit: return

        # Apply a specific conversion for each instrument
        if "GALEX" in self.image.filter.name: self.convert_galex()
        elif "SDSS" in self.image.filter.name: self.convert_sdss()
        elif "656_1" in self.image.filter.name: self.convert_ha()
        elif "2MASS" in self.image.filter.name: self.convert_2mass()
        elif "IRAC" in self.image.filter.name: self.convert_irac()
        elif "WISE" in self.image.filter.name: self.convert_wise()
        elif "Pacs" in self.image.filter.name: self.convert_pacs()
        elif "SPIRE" in self.image.filter.name: self.convert_spire()
        else: raise ValueError("Unknown image: " + self.image.name)

    # -----------------------------------------------------------------

    def apply(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Applying the unit conversion factor to the image ...")

        # Multiply the image (primary, errors and calibration_errors frame) by the conversion factor
        self.image *= self.conversion_factor

        # Set the new unit
        self.image.unit = self.target_unit

    # -----------------------------------------------------------------

    def spectral_factor_angstrom_to_hz(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return spectral_factor(wavelength, "Angstrom", "Hz")

    # -----------------------------------------------------------------

    def spectral_factor_cm_to_hz(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return spectral_factor(wavelength, "cm", "Hz")

    # -----------------------------------------------------------------

    def spectral_factor_micron_to_hz(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return spectral_factor(wavelength, "micron", "Hz")

    # -----------------------------------------------------------------

    def spectral_factor_meter_to_hz(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return spectral_factor(wavelength, "m", "Hz")

    # -----------------------------------------------------------------

    def pixel_factor(self, pixelscale):

        """
        This function ...
        :param pixelscale:
        :return:
        """

        return (1.0/pixelscale**2).to("pix2/sr").value

    # -----------------------------------------------------------------

    def convert_from_ergscmhz(self):

        """
        This function ...
        :return:
        """

        # Conversion from erg / [s * cm2 * Hz] (per pixel) to Jy (= 1e-23 erg / [s * cm2 * Hz]) (per pixel)
        self.conversion_factor *= ergscmHz_to_Jy

        # Convert from Jy per pixel to the target unit (MJy / sr)
        self.convert_from_jy()

    # -----------------------------------------------------------------

    def convert_from_jy(self):

        """
        This function ...
        :return:
        """

        # Conversion from Jy to MJy
        self.conversion_factor *= 1e-6

        # Conversion from MJy (per pixel2) to MJy / sr
        self.conversion_factor *= self.pixel_factor(self.image.primary.average_pixelscale)

    # -----------------------------------------------------------------

    def convert_galex(self):

        """
        This function ...
        :return:
        """

        # Conversion from count/s (per pixel) to erg / [s * cm2 * Angstrom] (per pixel2)
        self.conversion_factor *= galex_conversion_factors[self.image.filter.name]

        # Conversion from erg / [s * cm2 * Angstrom] (per pixel2) to erg / [s * cm2 * Hz] (per pixel)
        wavelength = effective_wavelengths[self.image.filter.name] # -> differs from self.image.frames.primary.filter.effectivewavelength()
        self.conversion_factor *= self.spectral_factor_angstrom_to_hz(wavelength)

        # Conversion from erg / [s * cm2 * Hz] per pixel2 to the target unit (MJy / sr)
        self.convert_from_ergscmhz()

    # -----------------------------------------------------------------

    def convert_sdss(self):

        """
        This function ...
        :return:
        """

        # Unit = nanomaggies
        # 1 nanomaggie = 3.613e-6 Jy

        # Conversion from nanomaggies to Jy (per pixel2)
        self.conversion_factor *= 3.613e-6

        # Conversion from Jy per pixel2 to the target unit (MJy / sr)
        self.convert_from_jy()

    # -----------------------------------------------------------------

    def convert_ha(self):

        """
        This function ...
        :return:
        """

        # What I thought the unit conversion should be:

        # Conversion from erg / [s * cm2] (per pixel2) to erg / [s * cm2 * micron] (per pixel2) --> divide by eff. bandwidth
        #self.conversion_factor *= 1.0 / self.image.filter.effective_bandwidth()

        # Conversion from erg / [s * cm2 * micron] (per pixel2) to erg / [s * cm2 * Hz] (per pixel2)
        #self.conversion_factor *= self.spectral_factor_micron_to_hz(self.image.wavelength)

        # What Ilse says it should be:

        # Get the frequency of Ha
        frequency = self.image.wavelength.to("Hz", equivalencies=spectral())

        # Conversion from erg / [s * cm2] (per pixel2) to erg / [s * cm2 * Hz] (per pixel2)
        self.conversion_factor *= 1.0 / frequency.value

        # Conversion from erg / [s * cm2 * Hz] per pixel2 to the target unit (MJy / sr)
        self.convert_from_ergscmhz()

        # Correct for contribution of NII (Bendo et. al 2015) NII/Halpha = 0.55 for M81 => divide by a factor of 1.55
        self.conversion_factor *= 1.0 / 1.55

    # -----------------------------------------------------------------

    def convert_2mass(self):

        """
        This function ...
        :return:
        """

        # Conversion from dimensionless flux (X=F/F0) to magnitude:
        #   mag(X) = m_0 (magnitude zero point) - 2.5 * log( X )
        # and then back to flux (Jy) with the flux zero-point:
        #   F(X) = F_0 * 10^(-1/2.5 * mag(X) )
        # Both equations come from: m - m_0 = -2.5 * log(F/F0)
        #   First couple (F_0, m_0):
        #    - F_0 "hidden" in data (data is relative to this value)
        #    - m_0: magnitude given in header
        #   Second couple (F_0, m_0):
        #    - F_0: can be found here: http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
        #    - m_0: the F_0 values found on the webpage above are for a magnitude of zero, so m_0 = 0
        # The two equations combine to:
        #  F(X) = F_0 * 10^(-m_0/2.5) * X
        m_0 = self.image.primary.zero_point
        f_0 = f_0_2mass[self.image.filter.name]
        self.conversion_factor *= f_0 * np.power(10.0, -m_0/2.5)

        # Conversion from Jy per pixel2 to the target unit (MJy / sr)
        self.convert_from_jy()

    # -----------------------------------------------------------------

    def convert_irac(self):

        """
        This function ...
        :return:
        """

        # Multiplication with the appropriate extended source photometrical correction coefficients
        self.conversion_factor *= photometrical_correction_irac[self.image.filter.name]

    # -----------------------------------------------------------------

    def convert_wise(self):

        """
        This function ...
        :return:
        """

        # Conversion from dimensionless unit (DN, or count) to flux (in W / [cm2 * micron]) (per pixel)
        m_0 = self.image.primary.zero_point # rounded values in the header ??
        #m_0 = m_0_wise[self.image.filter.name] + absolute_calibration_wise[self.image.filter.name]
        if self.image.filter.name == "WISE.W4": m_0 += calibration_discrepancy_wise_w4
        f_0 = f_0_wise[self.image.filter.name]
        self.conversion_factor *= f_0 * np.power(10.0, -m_0/2.5)

        # Conversion from W / [cm2 * micron] (per pixel2) to W / [cm2 * Hz] (per pixel)
        self.conversion_factor *= self.spectral_factor_micron_to_hz(self.image.wavelength)

        # Conversion from W / [cm2 * Hz] (per pixel2) to W / [m2 * Hz] (per pixel2)
        self.conversion_factor *= 100.0**2

        # Conversion from W / [m2 * Hz] (per pixel2) to Jy per pixel (1 Jy = 1e-26 W / [m2 * Hz])
        self.conversion_factor *= 1e26

        # Conversion from Jy per pixel to the target unit (MJy / sr)
        self.convert_from_jy()

    # -----------------------------------------------------------------

    def convert_pacs(self):

        """
        This function ...
        :return:
        """

        # Conversion from Jy per pixel to the target unit (MJy / sr)
        self.convert_from_jy()

    # -----------------------------------------------------------------

    def convert_spire(self):

        """
        This function ...
        :return:
        """

        # Conversion from Jy per pixel to the target unit (MJy / sr)
        self.convert_from_jy()

# -----------------------------------------------------------------
