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
from pts.core.basics.log import log
from pts.core.units.unit import PhotometricUnit
from pts.core.units.parsing import parse_unit as u
from pts.core.units.parsing import parse_quantity
from pts.core.tools import filesystem as fs
from pts.core.data.sed import SED, ObservedSED
from pts.core.plot.sed import SEDPlotter
from pts.modeling.preparation.unitconversion import neutral_fluxdensity_to_jansky, si_to_jansky
from pts.core.data.sun import Sun
from pts.core.filter.filter import parse_filter
from pts.core.tools import types

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "testing the unit conversions"

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W") # 3.828×10^26 W

# -----------------------------------------------------------------

class UnitsTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(UnitsTest, self).__init__(*args, **kwargs)

        self.observed_sed = None
        self.model_sed_1 = None
        self.model_sed_2 = None
        self.model_sed_3 = None
        #self.model_sed_4 = None
        self.mock_sed = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Load the SEDs
        self.load_seds()

        # Simple test
        self.simple()

        # SED test
        self.test_sed()

        # Test with solar units
        self.test_solar()

        # Test with surface brightness
        self.test_brightness()

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

        # Inform the user
        log.info("Loading the SEDs ...")

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

        # Test interpretation of names such as Jansky
        self.test_names()

        # Test flux conversions
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

        # Compare
        self.compare_seds()

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

        #print(self.model_sed_1)

    # -----------------------------------------------------------------

    def compare_seds(self):

        """
        This function ...
        :return: 
        """

        for index in range(len(self.model_sed_3)):
            ratio = self.model_sed_1["Photometry"][index] / self.model_sed_3["Photometry"][index]
            #print(self.model_sed_3["Wavelength"][index], self.model_sed_1["Wavelength"][index], self.model_sed_3["Photometry"][index], self.model_sed_1["Photometry"][index], ratio)
            if not np.isclose(ratio, 1.0): print("FAIL")

    # -----------------------------------------------------------------

    def test_solar(self):

        """
        THis function ...
        :return: 
        """

        # Bolometric
        self.test_solar_bolometric()

        # Spectral
        self.test_solar_spectral()

    # -----------------------------------------------------------------

    def test_solar_bolometric(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        #log.info("Loading the frames ...")

        # Get the galaxy distance
        #distance = self.galaxy_properties.distance
        distance = 10. * u("Mpc")

        # Load all the frames and error maps
        #for name in names:

        #frame = self.dataset.get_frame(name)
        #errors = self.dataset.get_errormap(name)

        ## CONVERT TO LSUN

        # Get pixelscale and wavelength
        pixelscale = 1.5 * u("arcsec")
        wavelength = 70. * u("micron")

        ##

        # Get a quantity
        quantity = 1 * u("MJy/sr")

        # Conversion from MJy / sr to Jy / sr
        conversion_factor = 1e6

        # Conversion from Jy / sr to Jy / pix(2)
        conversion_factor *= (pixelscale ** 2).to("sr").value

        # Conversion from Jy / pix to W / (m2 * Hz) (per pixel)
        conversion_factor *= 1e-26

        # Conversion from W / (m2 * Hz) (per pixel) to W / (m2 * m) (per pixel)
        conversion_factor *= (speed_of_light / wavelength ** 2).to("Hz/m").value

        # Conversion from W / (m2 * m) (per pixel) [SPECTRAL FLUX] to W / m [SPECTRAL LUMINOSITY]
        conversion_factor *= (4. * np.pi * distance ** 2).to("m2").value

        # Conversion from W / m [SPECTRAL LUMINOSITY] to W [NEUTRAL SPECTRAL LUMINOSITY]
        conversion_factor *= wavelength.to("m").value

        # Conversion from W to Lsun
        conversion_factor *= 1. / solar_luminosity.to("W").value

        ## CONVERT

        #frame *= conversion_factor
        #frame.unit = "Lsun"
        quantity_manual = quantity * conversion_factor
        quantity_manual.unit = "Lsun"

        print(quantity_manual)

        quantity_automatic = quantity.to("Lsun", density=True, distance=distance, wavelength=wavelength, pixelscale=pixelscale)

        print(quantity_automatic)

        if np.isclose(quantity_manual.value, quantity_automatic.value) and quantity_manual.unit == quantity_automatic.unit: print("SOLAR: OK")
        else: print("SOLAR: FAIL!")

    # -----------------------------------------------------------------

    def test_solar_spectral(self):

        """
        This function ...
        :return: 
        """

        fuv_filter = parse_filter("FUV")
        i1_filter = parse_filter("I1")

        # Solar properties
        sun = Sun()
        #sun_fuv = sun.luminosity_for_filter_as_unit(fuv_filter)  # Get the luminosity of the Sun in the FUV band
        #sun_i1 = sun.luminosity_for_filter_as_unit(i1_filter)  # Get the luminosity of the Sun in the IRAC I1 band

        #print(sun_fuv)
        #print(sun_i1)

        fuv = sun.luminosity_for_filter(fuv_filter, unit="W/Hz")
        i1 = sun.luminosity_for_filter(i1_filter, unit="W/Hz")

        fuv_wav = sun.luminosity_for_wavelength(fuv_filter.wavelength, unit="W/Hz")
        i1_wav = sun.luminosity_for_wavelength(i1_filter.wavelength, unit="W/Hz")

        #print(fuv, fuv_wav)
        #print(i1, i1_wav)

        if np.isclose(fuv.value, fuv_wav.value, rtol=1e-2): print("OK")
        else: print("not OK: " + str(fuv) + " vs " + str(fuv_wav))

        if np.isclose(i1.value, i1_wav.value, rtol=1e-2): print("OK")
        else: print("not OK: " + str(i1) + " vs " + str(i1_wav))

    # -----------------------------------------------------------------

    def test_brightness(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing surface brightness units ...")

        units = [u("Jy"), u("W/micron"), u("Lsun"), u("erg/s/Hz"), u("W/sr"), u("Lsun/pc2", brightness=True), u("W/m2/micron")]

        print("")
        for unit in units:

            print(str(unit))
            print("")
            print(" - symbol: " + unit.symbol)
            print(" - physical type: " + unit.physical_type)
            print(" - base physical type: " + unit.base_physical_type)
            print(" - base unit: " + str(unit.base_unit))
            print(" - density: " + str(unit.density))
            print(" - brightness: " + str(unit.brightness))

            angular_area_unit = unit.corresponding_angular_area_unit
            #print(angular_area_unit)
            intrinsic_area_unit = unit.corresponding_intrinsic_area_unit

            print(" - corresponding angular area unit: " + str(angular_area_unit))
            print(" - corresponding intrinsic area unit: " + str(intrinsic_area_unit))
            print("")

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

def spectral_factor_hz_to_micron(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    wavelength_unit = "micron"
    frequency_unit = "Hz"

    # Convert string units to Unit objects
    if types.is_string_type(basestring): wavelength_unit = u(wavelength_unit)
    if types.is_string_type(frequency_unit): frequency_unit = u(frequency_unit)

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
    if types.is_string_type(wavelength_unit): wavelength_unit = u(wavelength_unit)
    if types.is_string_type(frequency_unit): frequency_unit = u(frequency_unit)

    conversion_factor_unit = wavelength_unit / frequency_unit

    # Calculate the conversion factor
    factor = (wavelength ** 2 / speed_of_light).to(conversion_factor_unit).value
    return 1./factor

# -----------------------------------------------------------------
