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
from astropy import units as u
from astropy import constants

# Import the relevant AstroMagic classes and modules
from ...magic.basics import Mask

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# Flux zero point for AB magnitudes
ab_mag_zero_point = 3631. * u.Unit("Jy")

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
effective_wavelengths = {"GALEX.FUV": 1528.0 * u.Unit("Angstrom"), "GALEX.NUV": 2271.0 * u.Unit("Angstrom")}

# -----------------------------------------------------------------

class UnitConverter(Configurable):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(UnitConverter, self).__init__(config, "modeling")

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
        self.target_unit = u.Unit(self.config.to_unit)

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
        #self.image.frames.primary = self.image.frames.primary * self.conversion_factor
        #if "errors" in self.image.frames: self.image.frames.errors = self.image.frames.errors * self.conversion_factor

        # Set the new unit
        self.image.set_unit(self.target_unit)

    # -----------------------------------------------------------------

    def spectral_factor_angstrom_to_hz(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return (wavelength**2 / speed_of_light).to("Angstrom / Hz").value

    # -----------------------------------------------------------------

    def spectral_factor_cm_to_hz(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return (wavelength**2 / speed_of_light).to("cm / Hz").value

    # -----------------------------------------------------------------

    def spectral_factor_micron_to_hz(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return (wavelength**2 / speed_of_light).to("micron / Hz").value

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
        self.conversion_factor *= self.pixel_factor(self.image.frames.primary.xy_average_pixelscale)

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
        frequency = self.image.wavelength.to("Hz", equivalencies=u.spectral())

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
        m_0 = self.image.frames.primary.zero_point
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
        m_0 = self.image.frames.primary.zero_point # rounded values in the header ??
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
