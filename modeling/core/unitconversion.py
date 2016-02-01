#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.unitconversion Contains the UnitConverter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
import astropy.units as u
from astropy import constants

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# 2MASS F_0 (in Jy)
f_0_2mass = {"2MASS J": 1594, "2MASS H": 1024, "2MASS K": 666.7}

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

    # -----------------------------------------------------------------

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup(image)

        # 2. Convert
        self.convert()

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
        self.log.info("Converting the unit of the image to " + str(self.target_unit) + " if necessary ...")

        # Skip the unit conversion for images that are already in the right unit
        if self.image.unit == self.target_unit: return

        # Apply a specific conversion for each instrument
        if "GALEX" in self.image.name: self.convert_galex()
        elif "SDSS" in self.image.name: self.convert_sdss()
        elif "Ha" in self.image.name: self.convert_ha()
        elif "2MASS" in self.image.name: self.convert_2mass()
        elif "WISE" in self.image.name: self.convert_wise()
        elif "PACS" in self.image.name: self.convert_pacs()
        elif "SPIRE" in self.image.name: self.convert_spire()
        else: raise ValueError("Unkown image: " + self.image.name)

    # -----------------------------------------------------------------

    def convert_galex(self):

        """
        This function ...
        :return:
        """

        #FUV: Flux [erg sec-1 cm-2 Å-1] = 1.40 x 10-15 x CPS
        #NUV: Flux [erg sec-1 cm-2 Å-1] = 2.06 x 10-16 x CPS

        # Get the pixelscale of the image
        pixelscale = self.image.frames.primary.xy_average_pixelscale.value

        # Get the wavelength of the image
        wavelength = self.image.frames.primary.filter.centerwavelength()
        wavelength = (wavelength * u.Unit("micron")).to("AA")

        # Speed of light in Angstrom per seconds
        c = speed_of_light.to("AA/s")
        spectralfactor = wavelength.value**2 / c.value
        pixelfactor = (206264.806247 / pixelscale)**2
        factor = 1.4e-15 * spectralfactor * 1e17 * pixelfactor

        # Multiply the image (primary and errors frame) by the conversion factor
        self.image *= factor

    # -----------------------------------------------------------------

    def convert_sdss(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def convert_ha(self):

        """
        This function ...
        :return:
        """

        # Get the pixelscale of the image
        pixelscale = self.image.frames.primary.xy_average_pixelscale.value

        pixelfactor = (206264.806247 / pixelscale)**2

        # The wavelength (in meter)
        wavelength = 0.657894736 * 1e-6

        # The speed of light (in m / s)
        c = 299792458

        # The frequency (in Hz)
        frequency = c / wavelength

        # Calculate the conversion factor
        factor = 1e23 * 1e-6 * pixelfactor / frequency

        # Multiply the image (primary and errors frame) by the conversion factor
        self.image *= factor

    # -----------------------------------------------------------------

    def convert_2mass(self):

        """
        This function ...
        :return:
        """

        # Conversion factor to magnitudes
        m_0 = self.image.frames.primary.zero_point

        # Conversion factor back to fluxes
        f_0 = f_0_2mass[self.image.name]

        # Get the pixelscale of the image
        pixelscale = self.image.frames.primary.xy_average_pixelscale.value

        # Calculate the conversion factor
        pixelfactor = (206264.806247 / pixelscale)**2
        factor = 10e-6 * pixelfactor
        factor *= f_0 * np.power(10.0, -m_0/2.5)

        # Multiply the image (primary and errors frame) by the conversion factor
        self.image *= factor

    # -----------------------------------------------------------------

    def convert_wise(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def convert_pacs(self):

        """
        This function ...
        :return:
        """

        # Get the pixelscale of the image
        pixelscale = self.image.frames.primary.xy_average_pixelscale.value

        # Calculate the conversion factor
        pixelfactor = (206264.806247 / pixelscale)**2
        factor = 1e-6 * pixelfactor

        # Multiply the image (primary and errors frame) by the conversion factor
        self.image *= factor

    # -----------------------------------------------------------------

    def convert_spire(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
