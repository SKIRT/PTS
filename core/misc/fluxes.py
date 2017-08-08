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

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..basics.log import log
from ..filter.broad import BroadBandFilter
from ..filter.filter import parse_filter
from ..data.sed import SED
from ...magic.services.spire import SPIRE
from ..units.parsing import parse_unit as u
from ..data.sed import ObservedSED
from ..tools import sequences
from ..basics.configurable import Configurable
from ..simulation.simulation import createsimulations
from ..tools import parsing
from ..tools import types
from ..tools import numbers

# -----------------------------------------------------------------

class ObservedFluxCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ObservedFluxCalculator, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation prefix
        self.simulation_prefix = None

        # The ski file
        self.ski = None

        # The paths to the SED files produced by SKIRT
        self.sed_paths = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # Instrument names
        self.instrument_names = None

        # The errors for the different filters
        self.errors = None

        # The filters for which the fluxes should be calculated
        self.filters = None

        # The output observed SEDs
        self.mock_seds = dict()

        # The SPIRE instance
        self.spire = SPIRE()

        # The SED paths
        self.paths = dict()

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
        errors = kwargs.pop("errors", None)

        # Obtain the paths to the SED files created by the simulation
        self.sed_paths = simulation.seddatpaths()

        # Get the simulation prefix
        self.simulation_prefix = simulation.prefix()

        # The ski file
        self.ski = simulation.ski_file

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Set the instrument names
        self.instrument_names = instrument_names

        # Set the errors
        if errors is not None:
            self.errors = dict()
            # Set the error bars from strings
            for filter_name in errors:
                error = errors[filter_name]
                if types.is_string_type(error): error = parsing.errorbar(error)
                self.errors[filter_name] = error

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

            # Set the conversion info
            conversion_info = dict()
            conversion_info["distance"] = self.ski.get_instrument_distance(instr_name)

            # Debugging
            log.debug("Loading the modelled SED ...")

            # Get the name of the SED
            sed_name = fs.name(sed_path).split("_sed")[0]

            # Load the modelled SED
            model_sed = SED.from_skirt(sed_path)

            # Debugging
            log.debug("Calculating the observed fluxes for the " + sed_name + " SED ...")

            # Convert model sed into observed SED
            mock_sed = create_mock_sed(model_sed, self.filters, self.spire, spectral_convolution=self.config.spectral_convolution, errors=self.errors)

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

            # Debugging
            log.debug("Writing the mock SED '" + name + "' to '" + path + "' ...")

            # Write out the flux table
            self.mock_seds[name].saveto(path)

            # Set the path
            self.paths[name] = path

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

def get_wavelength_and_fluxdensity_arrays(model_sed, conversion_info=None):

    """
    This function ...
    :param model_sed: 
    :param conversion_info:
    :return: 
    """

    # Get arrays
    wavelengths = model_sed.wavelengths("micron", asarray=True)
    fluxdensities = model_sed.photometry("W / (m2 * micron)", asarray=True, conversion_info=conversion_info)

    # Return the wavelengths and fluxdensities
    return wavelengths, fluxdensities

# -----------------------------------------------------------------

def get_wavelength_and_fluxdensity_arrays_old(model_sed):

    """
    This function ...
    :param model_sed: 
    :return: 
    """

    # Get the wavelengths and flux densities
    wavelengths = model_sed.wavelengths("micron", asarray=True)
    fluxdensities = []
    for wavelength, fluxdensity_jy in zip(model_sed.wavelengths("micron"), model_sed.photometry("Jy", conversion_info=conversion_info)):

        # 2 different ways should be the same:
        # fluxdensity_ = fluxdensity_jy.to("W / (m2 * micron)", equivalencies=spectral_density(wavelength))
        # fluxdensity = fluxdensity_jy.to("W / (m2 * Hz)").value * spectral_factor_hz_to_micron(wavelength) * u("W / (m2 * micron)")
        # print(fluxdensity_, fluxdensity) # IS OK!

        fluxdensity = fluxdensity_jy.to("W / (m2 * micron)", wavelength=wavelength)
        fluxdensities.append(fluxdensity.value)

        # Add to data points for black body (to get spectral index)
        if wavelength > 50. * u("micron"):

            bb_frequencies.append(wavelength.to("Hz", equivalencies=spectral()).value)
            bb_fluxdensities.append(fluxdensity_jy.to("W / (m2 * Hz)").value)  # must be frequency-fluxdensity

    # Convert to array
    fluxdensities = np.array(fluxdensities)  # in W / (m2 * micron)

# -----------------------------------------------------------------

def calculate_spectral_indices(model_sed):

    """
    This function ...
    :param model_sed: 
    :return: 
    """

    # Define minimum wavelength
    min_wavelength = 50. * u("micron")

    # Get frequencies and flux densities
    frequencies = model_sed.frequencies("Hz", min_wavelength=min_wavelength, asarray=True)
    fluxdensities = model_sed.photometry("W / (m2 * Hz)", min_wavelength=min_wavelength, asarray=True)

    # Take logarithm
    log_frequencies = np.log10(frequencies)
    log_fluxdensities = np.log10(fluxdensities)

    # Calculate the derivative of the log(Flux) to the frequency, make a 'spectral index' function
    #gradients = np.gradient(log_frequencies, log_fluxdensities) # DOESN'T WORK ANYMORE?
    #gradients = np.diff(log_fluxdensities) / np.diff(log_frequencies)
    new_frequencies, gradients = numbers.derivatives(log_frequencies, log_fluxdensities)

    # Create function
    spectral_indices = interp1d(new_frequencies, gradients)

    # Return the spectral indices
    return spectral_indices

# -----------------------------------------------------------------

def calculate_fluxdensity_convolution(fltr, wavelengths, fluxdensities, spectral_indices, spire, errors=None):

    """
    This function ...
    :param fltr:
    :param wavelengths:
    :param fluxdensities:
    :param spectral_indices:
    :param spire:
    :param errors:
    :return: 
    """

    # Debugging
    log.debug("Calculating the observed flux for the " + str(fltr) + " filter by convolving spectrally ...")

    # Calculate the flux: flux densities must be per wavelength instead of per frequency!
    fluxdensity = float(fltr.convolve(wavelengths, fluxdensities)) * u("W / (m2 * micron)")
    fluxdensity_value = fluxdensity.to("Jy", equivalencies=spectral_density(fltr.pivot)).value  # convert back to Jy

    # For SPIRE, also multiply with Kbeam correction factor
    if fltr.instrument == "SPIRE":

        # Calculate the spectral index for the simulated SED at this filter
        # Use the central wavelength (frequency)
        central_frequency = fltr.center.to("Hz", equivalencies=spectral()).value
        spectral_index_filter = spectral_indices(central_frequency)

        # Get the Kbeam factor
        kbeam = spire.get_kbeam_spectral(fltr, spectral_index_filter)

        # Multiply the flux density
        fluxdensity_value *= kbeam

    # Determine filter name
    filter_name = str(fltr)

    # Add a point to the mock SED
    error = errors[filter_name] if errors is not None and filter_name in errors else None

    # Return the fluxdensity and error
    return fluxdensity_value * u("Jy"), error

# -----------------------------------------------------------------

def calculate_fluxdensity_closest(fltr, wavelengths, fluxdensities, errors=None):

    """
    This function ...
    :param fltr: 
    :param wavelengths: 
    :param fluxdensities: 
    :param errors:
    :return: 
    """

    # Debugging
    log.debug("Getting the observed flux for the " + str(fltr) + " filter ...")

    # Get the index of the wavelength closest to that of the filter
    index = sequences.find_closest_index(wavelengths, fltr.pivot.to("micron").value)  # wavelengths are in micron

    # Get the flux density
    fluxdensity = fluxdensities[index] * u("W / (m2 * micron)")  # flux densities are in W/(m2 * micron)
    fluxdensity = fluxdensity.to("Jy", equivalencies=spectral_density(fltr.pivot))  # now in Jy

    # Determine filter name
    filter_name = str(fltr)

    # Add a point to the mock SED
    error = errors[filter_name] if errors is not None and filter_name in errors else None

    # Return the fluxdensity and error
    return fluxdensity, error

# -----------------------------------------------------------------

def create_mock_sed(model_sed, filters, spire, spectral_convolution=True, errors=None):

    """
    This function ...
    :param model_sed: 
    :param filters:
    :param spire:
    :param spectral_convolution:
    :param errors:
    :return: 
    """

    # Get the wavelengths and flux densities as arrays
    wavelengths, fluxdensities = get_wavelength_and_fluxdensity_arrays(model_sed)

    # Calculate the spectral indices for wavelengths in the dust regime
    spectral_indices = calculate_spectral_indices(model_sed)

    # Create an observed SED for the mock fluxes
    mock_sed = ObservedSED(photometry_unit="Jy")

    # Loop over the different filters
    for fltr in filters:

        # Broad band filter, with spectral convolution
        if isinstance(fltr, BroadBandFilter) and spectral_convolution:

            # Calculate
            fluxdensity, error = calculate_fluxdensity_convolution(fltr, wavelengths, fluxdensities, spectral_indices, spire, errors=errors)

            # Add a data point to the mock SED
            mock_sed.add_point(fltr, fluxdensity, error)

        # Broad band filter without spectral convolution or narrow band filter
        else:

            # Calculate
            fluxdensity, error = calculate_fluxdensity_closest(fltr, wavelengths, fluxdensities, errors=errors)

            # Add a data point to the mock SED
            mock_sed.add_point(fltr, fluxdensity, error)

    # Return the mock observed SED
    return mock_sed

# -----------------------------------------------------------------
