#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.curve Contains the Curve class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from scipy import interpolate

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .table import SmartTable
from ..tools import tables
from .filter import Filter

# -----------------------------------------------------------------

class Curve(SmartTable):

    """
    This class ...
    """

    column_info = []

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, x_unit=None, y_unit=None, x_name="x", y_name="y", x_description="x values", y_description="y values"):

        """
        This function ...
        :param x_unit:
        :param y_unit:
        :param x_name:
        :param y_name:
        :param x_description:
        :param y_description:
        :return:
        """

        # Add columns
        cls.column_info.append((x_name, float, str(x_unit), x_description))
        cls.column_info.append((y_name, float, str(y_unit), y_description))

        # Call the initialize function of the SmartTable table function
        return super(Curve, cls).initialize()

    # -----------------------------------------------------------------

    def add_point(self, x_value, y_value):

        """
        This function ...
        :param x_value:
        :param y_value:
        :return:
        """

        # Set values
        values = [x_value, y_value]

        # Add a row
        self.add_row(values)

        # Sort the table by the x values
        self.sort(self.colnames[0])

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the saveto function of the base class
        super(Curve, self).saveto(path)

    # -----------------------------------------------------------------

    @property
    def has_errors(self):

        """
        This function ...
        :return:
        """

        return "Error-" in self.colnames and "Error+" in self.colnames

# -----------------------------------------------------------------

class WavelengthCurve(Curve):

    """
    This function ...
    """

    @classmethod
    def initialize(cls, name, description, unit=None):

        """
        This class ...
        :param name:
        :param description:
        :param unit:
        :return:
        """

        # Call the initialize function of the base class
        return super(WavelengthCurve, cls).initialize("micron", unit, "Wavelength", name, "Wavelength", description)

    # -----------------------------------------------------------------

    @property
    def value_name(self):

        """
        This function ..
        :return:
        """

        return self.colnames[-1]

    # -----------------------------------------------------------------

    def value_for_wavelength(self, wavelength):

        """
        This function ...
        :return:
        """

        interpolated = interpolate.interp1d(self.wavelengths(unit="micron", asarray=True), self.values(asarray=True), kind='linear')
        return interpolated(wavelength.to("micron").value)

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self["Wavelength"], unit=unit)
        else: return tables.column_as_list(self["Wavelength"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def values(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self[self.value_name], unit=unit)
        else: return tables.column_as_list(self[self.value_name], unit=unit, add_unit=add_unit)

# -----------------------------------------------------------------

class FilterCurve(WavelengthCurve):

    """
    This function ...
    """

    column_info = ["Observatory", "Instrument", "Band"] # first three columns

    # -----------------------------------------------------------------

    def add_point(self, fltr, value):

        """
        This function ...
        :param fltr:
        :param value:
        :return:
        """

        values = [fltr.observatory, fltr.instrument, fltr.band, fltr.pivotwavelength(), value]
        self.add_row(values)

    # -----------------------------------------------------------------

    def value_for_band(self, instrument, band, unit=None, add_unit=True):

        """
        This function ...
        :param instrument:
        :param band:
        :param unit:
        :param add_unit:
        :return:
        """

        has_unit = self[self.value_name].unit is not None
        has_mask = hasattr(self[self.value_name], "mask")

        # If unit has to be converted, check whether the original unit is specified
        if not has_unit and unit is not None: raise ValueError("Cannot determine the unit of quantity so values cannot be converted to " + str(unit))

        # Loop over all the entries in the table
        for i in range(len(self)):

            instrument_entry = self["Instrument"][i]
            band_entry = self["Band"][i]

            if not (instrument_entry == instrument and band_entry == band): continue

            # The column has a unit, we can convert if necessary
            if has_unit:

                # Add the unit initially to be able to convert
                value = self[self.value_name][i] * self[self.value_name].unit

                # If a target unit is specified, convert
                if unit is not None: value = value.to(unit).value * Unit(unit)

                # Strip unit if requested
                if not add_unit: value = value.value

            # No unit for the column
            else: value = self[self.value_name][i]

            # Return the value / quantity ...
            return value

        # If no match is found, return None
        return None

    # -----------------------------------------------------------------

    def value_for_filter(self, fltr, unit=None, add_unit=True):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.value_for_band(fltr.instrument, fltr.band, unit, add_unit)

    # -----------------------------------------------------------------

    def instruments(self):

        """
        This function ...
        :return:
        """

        return tables.column_as_list(self["Instrument"])

    # -----------------------------------------------------------------

    def bands(self):

        """
        This function ...
        :return:
        """

        return tables.column_as_list(self["Band"])

    # -----------------------------------------------------------------

    def filter_names(self):

        """
        This function ...
        :return:
        """

        # Initialize
        names = []

        # Loop over all entries
        for i in range(len(self)):

            # Get the instrument and band
            instrument = self["Instrument"][i]
            band = self["Band"][i]

            # Add the filter name
            names.append(instrument + " " + band)

        # Return the list of filter names
        return names

    # -----------------------------------------------------------------

    def filters(self):

        """
        This function ...
        :return:
        """

        # Initialize
        filters = []

        # Loop over all entries
        for i in range(len(self)):

            # Get the instrument and band
            instrument = self["Instrument"][i]
            band = self["Band"][i]

            # Create the filter
            fltr = Filter.from_instrument_and_band(instrument, band)

            # Add the filter to the list
            filters.append(fltr)

        # Return the list of filters
        return filters

# -----------------------------------------------------------------
