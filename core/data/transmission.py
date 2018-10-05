#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.data.transmission Contains the TransmissionCurve class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import interpolate

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..tools import tables, arrays
from ..basics.curve import Curve
from ..units.parsing import parse_unit as u

# -----------------------------------------------------------------

class TransmissionCurve(object):

    """
    This class ...
    """

    def __init__(self, wavelengths=None, transmissions=None, unit="micron"):

        """
        This function ...
        :param wavelengths:
        :param transmissions:
        :param unit:
        :return:
        """

        # Column names
        names = ["Wavelength", "Transmission"]

        # Create the table
        if wavelengths is None or transmissions is None: self.table = Table(names=names, dtype=('f8', 'f8'))
        else: self.table = tables.new([wavelengths, transmissions], names)

        # Set column units
        self.table["Wavelength"].unit = u(unit)

    # -----------------------------------------------------------------

    @classmethod
    def from_filter(cls, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get the wavelengths and transmissions, in micron
        wavelengths = fltr._Wavelengths
        transmissions = fltr._Transmission

        # Make sure that the transmission curves are not floating, but are attached to the x axis
        transmissions[0] = 0.0
        transmissions[-1] = 0.0

        # Create a new TransmissionCurve instance
        return cls(wavelengths, transmissions, unit="micron")

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new transmission curve
        transmission_curve = cls()

        # Load the table
        table = tables.from_file(path, format="ascii.ecsv")

        # Set the table
        transmission_curve.table = table

        # Return the transmission curve
        return transmission_curve

    # -----------------------------------------------------------------

    def add_entry(self, wavelength, transmission):

        """
        This function ...
        :param wavelength:
        :param transmission:
        :return:
        """

        self.table.add_row([wavelength.to("micron").value, transmission])

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return arrays.plain_array(self.table["Wavelength"], unit=unit, array_unit=self.table["Wavelength"].unit)
        else: return arrays.array_as_list(self.table["Wavelength"], unit=unit, add_unit=add_unit, array_unit=self.table["Wavelength"].unit)

    # -----------------------------------------------------------------

    def transmissions(self, asarray=False):

        """
        This function ...
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self.table["Transmission"])
        else: return arrays.array_as_list(self.table["Transmission"])

    # -----------------------------------------------------------------

    def transmission_at(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        interpolated = interpolate.interp1d(self.wavelengths(unit="micron", asarray=True), self.transmissions(asarray=True), kind='linear')
        return interpolated(wavelength.to("micron").value)

    # -----------------------------------------------------------------

    def normalize(self, value=1.0, method="integral"):

        """
        This function ...
        :param value:
        :param method:
        :return:
        """

        if method == "max":

            max_transmission = np.max(self.table["Transmission"])
            factor = value / max_transmission
            self.table["Transmission"] *= factor

        elif method == "integral": raise NotImplementedError("Not implemented yet")
        else: raise ValueError("Invalid option for 'method'")

    # -----------------------------------------------------------------

    def normalize_at(self, wavelength, value=1.):

        """
        This function ...
        :param wavelength:
        :param value:
        :return:
        """

        transmission_wavelength = self.transmission_at(wavelength)
        self.table["Transmission"] *= (value / transmission_wavelength)

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):

        """
        This function ...
        :return:
        """

        index = np.argmin(self.table["Wavelength"])
        return self.table["Wavelength"][index]

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):

        """
        This function ...
        :return:
        """

        index = np.argmax(self.table["Wavelength"])
        return self.table["Wavelength"][index]

    # -----------------------------------------------------------------

    @property
    def min_transmission(self):

        """
        This function ...
        :return:
        """

        index = np.argmin(self.table["Transmission"])
        return self.table["Transmission"][index]

    # -----------------------------------------------------------------

    @property
    def max_transmission(self):

        """
        This function ...
        :return:
        """

        index = np.argmax(self.table["Transmission"])
        return self.table["Transmission"][index]

    # -----------------------------------------------------------------

    @property
    def peak_wavelength(self):

        """
        This function ...
        :return:
        """

        index = np.argmax(self.table["Transmission"])
        return self.table["Wavelength"][index]

    # -----------------------------------------------------------------

    @property
    def mean_wavelength(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    @property
    def median_wavelength(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    @property
    def wavelength_range(self):

        """
        This function ...
        :return:
        """

        from ..basics.range import QuantityRange
        return QuantityRange(self.min_wavelength, self.max_wavelength, unit=self.table["Wavelength"].unit)

    # -----------------------------------------------------------------

    def in_range(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return wavelength in self.wavelength_range

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Sort the table by wavelength
        self.table.sort("Wavelength")

        # Write the transmission data
        tables.write(self.table, path, format="ascii.ecsv")

# -----------------------------------------------------------------
