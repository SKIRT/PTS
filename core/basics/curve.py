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

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        if "x_unit" in kwargs: from_astropy = False
        else: from_astropy = True

        if not from_astropy:

            # Get properties
            x_unit = kwargs.pop("x_unit", None)
            y_unit = kwargs.pop("y_unit", None)
            x_name = kwargs.pop("x_name", "x")
            y_name = kwargs.pop("y_name", "y")
            x_description = kwargs.pop("x_description", "x values")
            y_description = kwargs.pop("y_description", "y values")

        # Call the constructor of the base class
        super(Curve, self).__init__(*args, **kwargs)

        if not from_astropy:

            # Set the column info
            self.column_info.append((x_name, float, str(x_unit), x_description))
            self.column_info.append((y_name, float, str(y_unit), y_description))

            # Set x name and y name
            self.x_name = x_name
            self.y_name = y_name

    # -----------------------------------------------------------------

    @property
    def x_unit(self):

        """
        This function ...
        :return:
        """

        return self[self.x_name].unit

    # -----------------------------------------------------------------

    @property
    def y_unit(self):

        """
        This function ...
        :return:
        """

        return self[self.y_name].unit

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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        if "y_name" in kwargs: from_astropy = False
        else: from_astropy = True

        if not from_astropy:

            # Get properties
            name = kwargs.pop("y_name")
            description = kwargs.pop("y_description")
            unit = kwargs.pop("y_unit", None)

            x_unit = "micron"
            y_unit = unit
            x_name = "Wavelength"
            y_name = name
            x_description = "Wavelength"
            y_description = description

            kwargs["x_unit"] = x_unit
            kwargs["y_unit"] = y_unit
            kwargs["x_name"] = x_name
            kwargs["y_name"] = y_name
            kwargs["x_description"] = x_description
            kwargs["y_description"] = y_description

        # Call the constructor of the base class
        super(WavelengthCurve, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def value_name(self):

        """
        This function ..
        :return:
        """

        return self.colnames[-1]

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        return self.y_unit

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):

        """
        This function ...
        :return:
        """

        return self.x_unit

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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FilterCurve, self).__init__(*args, **kwargs)

        # Set column info
        self.column_info.insert(0, ("Band", str, None, "band"))
        self.column_info.insert(0, ("Instrument", str, None, "instrument"))
        self.column_info.insert(0, ("Observatory", str, None, "observatory"))

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
