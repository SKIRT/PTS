#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.attenuation Contains the AttenuationCurve class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import interpolate

# Import astronomical modules
from astropy.units import Unit, spectral
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools import tables

# -----------------------------------------------------------------

class AttenuationCurve(object):

    """
    This class ...
    """

    def __init__(self, wavelengths=None, attenuations=None):

        """
        This function ...
        :param wavelengths:
        :param attenuations:
        :return:
        """

        # Column names
        names = ["Wavelength", "Attenuation"]

        # Create the table
        if wavelengths is None or attenuations is None: self.table = Table(names=names, dtype=('f8', 'f8'))
        else: self.table = tables.new([wavelengths, attenuations], names)

        # Set column units
        self.table["Wavelength"].unit = Unit("micron")

    # -----------------------------------------------------------------

    @classmethod
    def from_seds(cls, total, transparent):

        """
        This function ...
        :param total:
        :param transparent:
        :return:
        """

        # Get the wavelengths
        wavelengths = total.wavelengths(unit="micron", add_unit=False)

        # Get the total and transparent fluxes
        total_fluxes = total.fluxes(asarray=True)
        transparent_fluxes = transparent.fluxes(asarray=True)

        # Calculate the attenuations
        attenuations = -2.5 * np.log10(total_fluxes / transparent_fluxes)

        # Create a new AttenuationCurve instance
        return cls(wavelengths, attenuations)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new attenuation curve
        attenuation_curve = cls()

        names = ["Wavelength", "Attenuation", "Band", "Wavelength", "Flux", "Error-", "Error+"]
        observatory_column, instrument_column, band_column, wavelength_column, flux_column, error_min_column, error_plus_column = np.loadtxt(path, unpack=True, dtype=str)
        wavelength_column = wavelength_column.astype(float)
        flux_column = flux_column.astype(float)
        error_min_column = error_min_column.astype(float)
        error_plus_column = error_plus_column.astype(float)
        attenuation_curve.table = tables.new([observatory_column, instrument_column, band_column, wavelength_column, flux_column, error_min_column, error_plus_column], names)
        attenuation_curve.table["Wavelength"].unit = "micron"

        # Return the attenuation curve
        return attenuation_curve

    # -----------------------------------------------------------------

    def add_entry(self, wavelength, attenuation):

        """
        This function ...
        :param wavelength:
        :param attenuation:
        :return:
        """

        self.table.add_row([wavelength.to("micron").value, attenuation])

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Wavelength"], unit=unit)
        else: return tables.column_as_list(self.table["Wavelength"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def attenuations(self, asarray=False):

        """
        This function ...
        :param asarray:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Attenuation"])
        else: return tables.column_as_list(self.table["Attenuation"])

    # -----------------------------------------------------------------

    def attenuation_at(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        interpolated = interpolate.interp1d(self.wavelengths(unit="micron", asarray=True), self.attenuations(asarray=True), kind='linear')
        return interpolated(wavelength.to("micron").value)

    # -----------------------------------------------------------------

    def normalize_at(self, wavelength, value=1.):

        """
        This function ...
        :param wavelength:
        :param value:
        :return:
        """

        attenuation_wavelength = self.attenuation_at(wavelength)
        self.table["Attenuation"] /= attenuation_wavelength * value

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Sort the table by wavelength
        self.table.sort("Wavelength")

        # Write the attenuation data
        tables.write(self.table, path)

# -----------------------------------------------------------------
