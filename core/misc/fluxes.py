#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.fluxes Contains the ObservedFluxCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy.interpolate import interp1d

# Import astronomical modules
from astropy.units import spectral_density, spectral
from astropy import constants

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools.logging import log
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter
from ..filter.filter import parse_filter
from ..data.sed import SED
from ...magic.misc.spire import SPIRE
from ..basics.unit import parse_unit as u
from ..data.sed import ObservedSED
from ..tools import lists
from ..basics.configurable import Configurable
from ..simulation.simulation import createsimulations

# -----------------------------------------------------------------

class ObservedFluxCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ObservedFluxCalculator, self).__init__(config)

        # -- Attributes --

        # The simulation prefix
        self.simulation_prefix = None

        # The paths to the SED files produced by SKIRT
        self.sed_paths = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # Instrument names
        self.instrument_names = None

        # The filters for which the fluxes should be calculated
        self.filters = None

        # The output observed SEDs
        self.mock_seds = dict()

        # The SPIRE instance
        self.spire = SPIRE()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the filters
        self.create_filters()

        # 3. Calculate the observed fluxes
        self.calculate()

        # 4. Write the results
        if self.config.output is not None: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ObservedFluxCalculator, self).setup(**kwargs)

        # Get properties
        if "simulation" in kwargs: simulation = kwargs.pop("simulation")
        elif "simulation_output_path" in kwargs: simulation = createsimulations(kwargs.pop("simulation_output_path"), single=True)
        else: raise ValueError("Either 'simulation' or 'simulation_output_path' must be specified")
        output_path = kwargs.pop("output_path", None)
        filter_names = kwargs.pop("filter_names", None)
        instrument_names = kwargs.pop("instrument_names", None)

        # Obtain the paths to the SED files created by the simulation
        self.sed_paths = simulation.seddatpaths()

        # Get the simulation prefix
        self.simulation_prefix = simulation.prefix()

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Set the instrument names
        self.instrument_names = instrument_names

        # Set output path
        self.config.output = output_path

    # -----------------------------------------------------------------

    def create_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Constructing the filter objects ...")

        # Initialize the list
        self.filters = []

        # Loop over the different filter names
        for filter_name in self.filter_names:

            # Debugging
            log.debug("Constructing the " + filter_name + " filter ...")

            # Create the filter
            fltr = parse_filter(filter_name)

            # Add the filter to the list
            self.filters.append(fltr)

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the observed fluxes ...")

        # Loop over the different SEDs
        for sed_path in self.sed_paths:

            # Get the name of the instrument
            instr_name = instrument_name(sed_path, self.simulation_prefix)

            # If a list of instruments is defined an this instrument is not in this list, skip it
            if self.instrument_names is not None and instr_name not in self.instrument_names: continue

            # Get the name of the SED
            sed_name = fs.name(sed_path).split("_sed")[0]

            # Debugging
            log.debug("Calculating the observed fluxes for the " + sed_name + " SED ...")

            # Create an observed SED for the mock fluxes
            mock_sed = ObservedSED(photometry_unit="Jy")

            # Load the modelled SED
            model_sed = SED.from_skirt(sed_path)

            # Initialize lists
            bb_frequencies = []
            bb_fluxdensities = []

            # Get the wavelengths and flux densities
            wavelengths = model_sed.wavelengths("micron", asarray=True)
            fluxdensities = []
            for wavelength, fluxdensity_jy in zip(model_sed.wavelengths("micron"), model_sed.photometry("Jy")):

                # 2 different ways should be the same:
                #fluxdensity_ = fluxdensity_jy.to("W / (m2 * micron)", equivalencies=spectral_density(wavelength))
                fluxdensity = fluxdensity_jy.to("W / (m2 * Hz)").value * spectral_factor_hz_to_micron(wavelength) * u("W / (m2 * micron)")
                #print(fluxdensity_, fluxdensity) # IS OK!
                fluxdensities.append(fluxdensity.to("W / (m2 * micron)").value)

                if wavelength > 50. * u("micron"):
                    bb_frequencies.append(wavelength.to("Hz", equivalencies=spectral()).value)
                    bb_fluxdensities.append(fluxdensity_jy.to("W / (m2 * Hz)").value) # must be frequency-fluxdensity

            # Convert to array
            fluxdensities = np.array(fluxdensities)  # in W / (m2 * micron)

            # Calculate the derivative of the log(Flux) to the frequency, make a 'spectral index' function
            bb_frequencies = np.array(bb_frequencies)
            bb_fluxdensities = np.array(bb_fluxdensities)
            gradients = np.gradient(np.log10(bb_frequencies), np.log10(bb_fluxdensities))
            spectral_index = interp1d(bb_frequencies, gradients)

            # Loop over the different filters
            for fltr in self.filters:

                # Broad band filter, with spectral convolution
                if isinstance(fltr, BroadBandFilter) and self.config.spectral_convolution:

                    # Debugging
                    log.debug("Calculating the observed flux for the " + str(fltr) + " filter by convolving spectrally ...")

                    # Calculate the flux: flux densities must be per wavelength instead of per frequency!
                    fluxdensity = float(fltr.convolve(wavelengths, fluxdensities)) * u("W / (m2 * micron)")
                    fluxdensity_value = fluxdensity.to("Jy", equivalencies=spectral_density(fltr.pivot)).value # convert back to Jy

                    # For SPIRE, also multiply with Kbeam correction factor
                    if fltr.instrument == "SPIRE":

                        # Calculate the spectral index for the simulated SED at this filter
                        # Use the central wavelength (frequency)
                        central_frequency = fltr.center.to("Hz", equivalencies=spectral()).value
                        spectral_index_filter = spectral_index(central_frequency)

                        # Get the Kbeam factor
                        kbeam = self.spire.get_kbeam_spectral(fltr, spectral_index_filter)

                        # Multiply the flux density
                        fluxdensity_value *= kbeam

                    # Add a point to the mock SED
                    mock_sed.add_point(fltr, fluxdensity_value * u("Jy"))

                # Broad band filter without spectral convolution or narrow band filter
                else:

                    # Debugging
                    log.debug("Getting the observed flux for the " + str(fltr) + " filter ...")

                    # Get the index of the wavelength closest to that of the filter
                    index = lists.find_closest_index(wavelengths, fltr.pivot.to("micron").value)  # wavelengths are in micron

                    # Get the flux density
                    fluxdensity = fluxdensities[index] * u("W / (m2 * micron)")  # flux densities are in W/(m2 * micron)
                    fluxdensity = fluxdensity.to("Jy", equivalencies=spectral_density(fltr.pivot))  # now in Jy

                    # Add a point to the mock SED
                    mock_sed.add_point(fltr, fluxdensity)

                # Narrow band filter
                #elif isinstance(fltr, NarrowBandFilter):

                # Unrecognized
                #else: raise ValueError("Unrecognized filter object: " + str(fltr))

            # Add the complete SED to the dictionary (with the SKIRT SED name as key)
            self.mock_seds[sed_name] = mock_sed

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the mock SEDs
        self.write_mock_seds()

    # -----------------------------------------------------------------

    def write_mock_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mock observed SED ...")

        # Loop over the different flux tables
        for name in self.mock_seds:

            # Determine the path to the output flux table
            path = self.output_path_file(name + "_fluxes.dat")

            # Write out the flux table
            self.mock_seds[name].saveto(path)

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

def instrument_name(sed_path, prefix):

    """
    This function ...
    :param sed_path:
    :param prefix:
    :return:
    """

    return fs.name(sed_path).split("_sed.dat")[0].split(prefix + "_")[1]

# -----------------------------------------------------------------
