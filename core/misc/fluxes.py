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

# Import astronomical modules
from astropy.units import Unit, spectral_density
from astropy import constants

# Import the relevant PTS classes and modules
from ..tools import tables, filesystem
from ..tools.logging import log
from ..basics.filter import Filter
from ...modeling.core.sed import SED

# -----------------------------------------------------------------

class ObservedFluxCalculator(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(ObservedFluxCalculator, self).__init__()

        # -- Attributes --

        # The paths to the SED files produced by SKIRT
        self.sed_paths = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # The filters for which the fluxes should be calculated
        self.filters = None

        # The flux tables for the different output SEDs
        self.tables = dict()

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None, filter_names=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :param filter_names:
        :return:
        """

        # Obtain the paths to the SED files created by the simulation
        self.sed_paths = simulation.seddatpaths()

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Create the filters
        self.create_filters()

        # Calculate the observed fluxes
        self.calculate()

        # Write the results
        if output_path is not None: self.write(output_path)

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
            filter = Filter.from_string(filter_name)

            # Add the filter to the list
            self.filters.append(filter)

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

            # Get the name of the SED
            sed_name = filesystem.name(sed_path).split("_sed")[0]

            # Debugging
            log.debug("Calculating the observed fluxes for the " + sed_name + " SED ...")

            # Create a flux table for this SED
            names = ["Observatory", "Instrument", "Band", "Wavelength", "Flux"]
            dtypes = ["S10", "S10", "S10", "f8", "f8"]
            data = [[], [], [], [], []]
            table = tables.new(data, names, dtypes=dtypes)
            table["Wavelength"].unit = "micron"
            table["Flux"].unit = "Jy"

            # Load the simulated SED
            sed = SED.from_file(sed_path)

            # Get the wavelengths and flux densities
            wavelengths = sed.wavelengths("micron", asarray=True)
            fluxdensities = []
            for wavelength, fluxdensity_jy in zip(sed.wavelengths("micron"), sed.fluxes("Jy")):

                # 2 different ways should be the same:
                #fluxdensity_ = fluxdensity_jy.to("W / (m2 * micron)", equivalencies=spectral_density(wavelength))
                fluxdensity = fluxdensity_jy.to("W / (m2 * Hz)").value * spectral_factor_hz_to_micron(wavelength) * Unit("W / (m2 * micron)")
                #print(fluxdensity_, fluxdensity) # IS OK!
                fluxdensities.append(fluxdensity.to("W / (m2 * micron)").value)

            fluxdensities = np.array(fluxdensities) # in W / (m2 * micron)

            # Loop over the different filters
            for filter in self.filters:

                # Debugging
                log.debug("Calculating the observed flux for the " + filter.name + " filter ...")

                # Calculate the flux: flux densities must be per wavelength instead of per frequency!
                fluxdensity = float(filter.convolve(wavelengths, fluxdensities)) * Unit("W / (m2 * micron)")
                fluxdensity_value = fluxdensity.to("Jy", equivalencies=spectral_density(filter.pivotwavelength() * Unit("micron"))).value # convert back to Jy

                # Add an entry to the flux table
                table.add_row([filter.observatory, filter.instrument, filter.band, filter.pivotwavelength(), fluxdensity_value])

            # Add the complete table to the dictionary (with the SKIRT SED name as key)
            self.tables[sed_name] = table

    # -----------------------------------------------------------------

    def write(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Writing the observed flux table ...")

        # Loop over the different flux tables
        for name in self.tables:

            # Determine the path to the output flux table
            path = filesystem.join(output_path, name + "_fluxes.dat")

            # Write out the flux table
            tables.write(self.tables[name], path)

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
    if isinstance(wavelength_unit, basestring): wavelength_unit = Unit(wavelength_unit)
    if isinstance(frequency_unit, basestring): frequency_unit = Unit(frequency_unit)

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
    if isinstance(wavelength_unit, basestring): wavelength_unit = Unit(wavelength_unit)
    if isinstance(frequency_unit, basestring): frequency_unit = Unit(frequency_unit)

    conversion_factor_unit = wavelength_unit / frequency_unit

    # Calculate the conversion factor
    factor = (wavelength ** 2 / speed_of_light).to(conversion_factor_unit).value
    return 1./factor

# -----------------------------------------------------------------
