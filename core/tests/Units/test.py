#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect
import numpy as np

# Import astronomical modules
from astropy import constants
from astropy.units import spectral

# Import the relevant PTS classes and modules
from pts.core.test.implementation import TestImplementation
from pts.core.tools.logging import log
from pts.core.basics.unit import PhotometricUnit
from pts.core.basics.unit import parse_unit as u
from pts.core.basics.quantity import parse_quantity
from pts.core.tools import filesystem as fs
from pts.core.data.sed import SED, ObservedSED
from pts.core.plot.sed import SEDPlotter
from pts.modeling.preparation.unitconversion import neutral_fluxdensity_to_jansky, si_to_jansky

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "testing the unit conversions"

# -----------------------------------------------------------------

class UnitsTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(UnitsTest, self).__init__(config, interactive)

        self.observed_sed = None
        self.model_sed_1 = None
        self.model_sed_2 = None
        self.model_sed_3 = None
        #self.model_sed_4 = None
        self.mock_sed = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Load the SEDs
        self.load_seds()

        # Simple test
        self.simple()

        # SED test
        self.test_sed()

        # Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(UnitsTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_seds(self):

        """
        This function ...
        :return:
        """

        observed_sed_path = fs.join(this_dir_path, "DustPedia.dat")
        sim_sed_path_1 = fs.join(this_dir_path, "M81_earth_sed1.dat")
        sim_sed_path_2 = fs.join(this_dir_path, "M81_earth_sed2.dat")
        mock_sed_path = fs.join(this_dir_path, "M81_earth_fluxes.dat")

        self.observed_sed = ObservedSED.from_file(observed_sed_path)
        self.model_sed_1 = SED.from_skirt(sim_sed_path_1)
        self.model_sed_2 = SED.from_skirt(sim_sed_path_2)
        self.mock_sed = ObservedSED.from_file(mock_sed_path)

    # -----------------------------------------------------------------

    def simple(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Simple test ...")

        self.test_names()

        self.test_fluxes()

    # -----------------------------------------------------------------

    def test_names(self):

        """
        Tihs function ...
        :return:
        """

        # Create nanomaggy unit
        nanomaggy = PhotometricUnit("nMgy")

        # Conversion factor to Jansky
        factor = nanomaggy.conversion_factor("Jy")

        #print(factor, 3.613e-6)

        assert np.isclose(factor, 3.613e-6)

    # -----------------------------------------------------------------

    def test_fluxes(self):

        """
        This function ...
        :return:
        """

        # Create a quantity that is 1 W/m2 at 200 micron

        wavelength = parse_quantity("200 micron")

        frequency = wavelength.to("Hz", equivalencies=spectral())

        print("")
        print("wavelength:", wavelength)
        print("frequency:", frequency)
        print("")

        flux = parse_quantity("1 W/m2", density=True) # neutral flux density

        print("neutral flux density:", flux)
        print("wavelength flux density:", flux / wavelength)
        print("frequency flux density:", flux / frequency)
        print("")

        converted = flux.value / wavelength.value # to W / (m2 * micron)
        converted = converted / spectral_factor_hz_to_micron(wavelength)
        converted = converted * si_to_jansky * u("Jy") # in Jy

        converted3 = flux.value / frequency.value * si_to_jansky * u("Jy")  # first to W / [m2 * Hz] and then to Jy

        converted_auto = flux.to("Jy", wavelength=wavelength)

        print("")
        print("flux density in Jy:", converted)
        print("flux density in Jy (other way):", converted3)
        print("automatically converted flux density in Jy:", converted_auto)
        print("")

        print("conversion factor:", converted.value / flux.value)
        print("automatically determined conversion factor:", converted_auto.value / flux.value)
        print("")

    # -----------------------------------------------------------------

    def test_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing SED unit conversion ...")

        # Convert manually
        self.convert_sed_manual()

        # Convert automatic
        self.convert_sed_automatic()

    # -----------------------------------------------------------------

    def convert_sed_manual(self):

        """
        This function ...
        :return:
        """

        # 2 different ways should be the same:
        # fluxdensity_ = fluxdensity_jy.to("W / (m2 * micron)", equivalencies=spectral_density(wavelength))
        #fluxdensity = fluxdensity_jy.to("W / (m2 * Hz)").value * spectral_factor_hz_to_micron(wavelength) * u("W / (m2 * micron)")
        # print(fluxdensity_, fluxdensity) # IS OK!
        #fluxdensities.append(fluxdensity.to("W / (m2 * micron)").value)

        wavelengths = self.model_sed_1.wavelengths("micron", asarray=True) # in micron
        fluxes = self.model_sed_1.photometry("W/m2", asarray=True) #self.model_sed_1["Photometry"] # in W/m2 (lambda * F_Lambda)

        #print("wavelengths", wavelengths)
        #print("fluxes", fluxes)

        fluxes = fluxes / wavelengths # in W / (m2 * micron)
        fluxes = [flux / spectral_factor_hz_to_micron(wavelength * u("micron")) for flux, wavelength in zip(fluxes, wavelengths)] # in W / (m2 * Hz)
        fluxes = [flux * si_to_jansky for flux in fluxes]

        self.model_sed_3 = SED.from_arrays(wavelengths, fluxes, wavelength_unit="micron", photometry_unit="Jy")
        #print(self.model_sed_3)

        ## ANALOG WAY:

        #fluxes2 = [neutral_fluxdensity_to_jansky(fluxdensity, wavelength * u("micron")) for fluxdensity, wavelength in zip(self.model_sed_1["Photometry"], wavelengths)]

        #print("Jy fluxes", fluxes2)
        #self.model_sed_4 = SED.from_arrays(wavelengths, fluxes2, wavelength_unit="micron", photometry_unit="Jy")
        #print(self.model_sed_4)

    # -----------------------------------------------------------------

    def convert_sed_automatic(self):

        """
        This function ...
        :return:
        """

        self.model_sed_1.convert_to(photometry_unit="Jy")

        print(self.model_sed_1)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        plotter = SEDPlotter()

        plotter.add_sed(self.observed_sed, "Observation")
        plotter.add_sed(self.model_sed_1, "model 1")
        plotter.add_sed(self.model_sed_2, "model 2")
        plotter.add_sed(self.model_sed_3, "model 3")
        #plotter.add_sed(self.model_sed_4, "model 4")

        plotter.run()

# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    :return:
    """

    pass

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# -----------------------------------------------------------------

def spectral_factor_hz_to_micron(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    wavelength_unit = "micron"
    frequency_unit = "Hz"

    # Convert string units to Unit objects
    if isinstance(wavelength_unit, basestring): wavelength_unit = u(wavelength_unit)
    if isinstance(frequency_unit, basestring): frequency_unit = u(frequency_unit)

    conversion_factor_unit = wavelength_unit / frequency_unit

    # Calculate the conversion factor
    factor = (wavelength ** 2 / speed_of_light).to(conversion_factor_unit).value
    return 1. / factor

# -----------------------------------------------------------------

def spectral_factor_hz_to_meter(wavelength):

    """
    This function ...
    :return:
    """

    wavelength_unit = "m"
    frequency_unit = "Hz"

    # Convert string units to Unit objects
    if isinstance(wavelength_unit, basestring): wavelength_unit = u(wavelength_unit)
    if isinstance(frequency_unit, basestring): frequency_unit = u(frequency_unit)

    conversion_factor_unit = wavelength_unit / frequency_unit

    # Calculate the conversion factor
    factor = (wavelength ** 2 / speed_of_light).to(conversion_factor_unit).value
    return 1./factor

# -----------------------------------------------------------------
