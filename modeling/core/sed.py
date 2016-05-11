#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.sed Contains the SED class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit, spectral
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools import tables
from ...core.basics.errorbar import ErrorBar

# -----------------------------------------------------------------

class IntrinsicSED(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        # Attributes
        self.table = None

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

    def luminosities(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Luminosity"], unit=unit)
        else: return tables.column_as_list(self.table["Luminosity"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, skiprows=0):
        """
        This function ...
        :param path:
        :param skiprows:
        :return:
        """

        # Create a new SED
        sed = cls()

        wavelength_column, luminosity_column = np.loadtxt(path, dtype=float, unpack=True, skiprows=skiprows)
        sed.table = tables.new([wavelength_column, luminosity_column], ["Wavelength", "Luminosity"])
        sed.table["Wavelength"].unit = Unit("micron")
        sed.table["Luminosity"].unit = Unit("W/micron")

        # Return the SED
        return sed

# -----------------------------------------------------------------

class ObservedSED(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Attributes
        self.table = Table(names=["Observatory", "Instrument", "Band", "Wavelength", "Flux", "Error-", "Error+"], dtype=('S10', 'S10', 'S10', 'f8', 'f8', 'f8', 'f8'))
        self.table["Wavelength"].unit = Unit("micron")
        self.table["Flux"].unit = Unit("Jy")
        self.table["Error-"].unit = Unit("Jy")
        self.table["Error+"].unit = Unit("Jy")

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new observed SED
        sed = cls()


        #names = ["Observatory", "Instrument", "Band", "Wavelength", "Flux", "Error-", "Error+"]
        #observatory_column, instrument_column, band_column, wavelength_column, flux_column, error_min_column, error_plus_column = np.loadtxt(path, unpack=True, dtype=str)
        #wavelength_column = wavelength_column.astype(float)
        #flux_column = flux_column.astype(float)
        #error_min_column = error_min_column.astype(float)
        #error_plus_column = error_plus_column.astype(float)
        #sed.table = tables.new([observatory_column, instrument_column, band_column, wavelength_column, flux_column, error_min_column, error_plus_column], names)
        #sed.table["Wavelength"].unit = "micron"
        #sed.table["Flux"].unit = "Jy"
        #sed.table["Error-"].unit = "Jy"
        #sed.table["Error+"].unit = "Jy"

        # New
        sed.table = tables.from_file(path, format="ascii.ecsv")

        # Return the observed SED
        return sed

    # -----------------------------------------------------------------

    def add_entry(self, filter, flux, error):

        """
        This function ...
        :param filter:
        :param flux:
        :param error:
        :return:
        """

        self.table.add_row([filter.observatory, filter.instrument, filter.band, filter.pivotwavelength(), flux, error.lower, error.upper])

    # -----------------------------------------------------------------

    def instruments(self):

        """
        This function ...
        :return:
        """

        return tables.column_as_list(self.table["Instrument"])

    # -----------------------------------------------------------------

    def bands(self):

        """
        This function ...
        :return:
        """

        return tables.column_as_list(self.table["Band"])

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

    def fluxes(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Flux"], unit=unit)
        else: return tables.column_as_list(self.table["Flux"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def errors(self, unit=None, add_unit=True):

        """
        This function ...
        :param unit:
        :param add_unit:
        :return:
        """

        return tables.columns_as_objects([self.table["Error-"], self.table["Error+"]], ErrorBar, unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def flux_for_filter(self, filter, unit=None, add_unit=True):

        """
        This function ...
        :param filter:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.flux_for_band(filter.instrument, filter.band, unit, add_unit)

    # -----------------------------------------------------------------

    def flux_for_band(self, instrument, band, unit=None, add_unit=True):

        """
        This function ...
        :param instrument:
        :param band:
        :param unit:
        :param add_unit:
        :return:
        """

        has_unit = self.table["Flux"].unit is not None
        has_mask = hasattr(self.table["Flux"], "mask")

        # If unit has to be converted, check whether the original unit is specified
        if not has_unit and unit is not None: raise ValueError(
            "Cannot determine the unit of the flux column so values cannot be converted to " + str(unit))

        # Loop over all the entries in the table
        for i in range(len(self.table)):

            instrument_entry = self.table["Instrument"][i]
            band_entry = self.table["Band"][i]

            if not (instrument_entry == instrument and band_entry == band): continue

            if has_unit:

                # Add the unit initially to be able to convert
                flux = self.table["Flux"][i] * self.table["Flux"].unit

                # If a target unit is specified, convert
                if unit is not None: flux = flux.to(unit).value * Unit(unit)

                if not add_unit: flux = flux.value

            else: flux = self.table["Flux"][i]

            return flux

        # If no match is found, return None
        return None

    # -----------------------------------------------------------------

    def error_for_filter(self, filter, unit=None, add_unit=True):

        """
        This function ...
        :param filter:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.error_for_band(filter.instrument, filter.band, unit, add_unit)

    # -----------------------------------------------------------------

    def error_for_band(self, instrument, band, unit=None, add_unit=True):

        """
        This function ...
        :param instrument:
        :param band:
        :param unit:
        :param add_unit:
        :return:
        """

        has_unit = self.table["Error-"].unit is not None and self.table["Error+"].unit is not None
        has_mask = hasattr(self.table["Error-"], "mask")
        assert has_mask == hasattr(self.table["Error+"], "mask")

        # If unit has to be converted, check whether the original unit is specified
        if not has_unit and unit is not None: raise ValueError(
            "Cannot determine the unit of the error columns so values cannot be converted to " + str(unit))

        # Loop over all the entries in the table
        for i in range(len(self.table)):

            instrument_entry = self.table["Instrument"][i]
            band_entry = self.table["Band"][i]

            if not (instrument_entry == instrument and band_entry == band): continue

            if has_unit:

                # Add the unit initially to be able to convert
                error_min = self.table["Error-"][i] * self.table["Error-"].unit
                error_plus = self.table["Error+"][i] * self.table["Error+"].unit

                # If a target unit is specified, convert
                if unit is not None:

                    error_min = error_min.to(unit).value * Unit(unit)
                    error_plus = error_plus.to(unit).value * Unit(unit)

                if not add_unit:

                    error_min = error_min.value
                    error_plus = error_plus.value

                error = ErrorBar(error_min, error_plus)

            else: error = ErrorBar(self.table["Error-"][i], self.table["Error+"][i])

            return error

        # If no match is found, return None
        return None

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Sort the table by wavelength
        self.table.sort("Wavelength")

        # Write the observed SED
        tables.write(self.table, path, format="ascii.ecsv")

# -----------------------------------------------------------------

class SED(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Attributes
        self.table = None

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

    def fluxes(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Flux"], unit=unit)
        else: return tables.column_as_list(self.table["Flux"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def errors(self, unit=None, add_unit=True):

        """
        This function ...
        :param unit:
        :param add_unit:
        :return:
        """

        return tables.columns_as_objects([self.table["Error-"], self.table["Error+"]], ErrorBar, unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    @property
    def has_errors(self):

        """
        This function ...
        :return:
        """

        return "Error-" in self.table.colnames and "Error+" in self.table.colnames

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new SED
        sed = cls()

        # New
        sed.table = tables.from_file(path, format="ascii.ecsv")

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt(cls, path, skiprows=0, contribution="total"):

        """
        This function ...
        :param path:
        :param skiprows:
        :param contribution:
        :return:
        """

        # Create a new SED
        sed = cls()

        # SEDInstrument:
        # column 1: lambda (micron)
        # column 2: total flux; lambda*F_lambda (W/m2)

        # From FullInstrument:
        # column 1: lambda (micron)
        # column 2: total flux; lambda*F_lambda (W/m2)
        # column 3: direct stellar flux; lambda*F_lambda (W/m2)
        # column 4: scattered stellar flux; lambda*F_lambda (W/m2)
        # column 5: total dust emission flux; lambda*F_lambda (W/m2)
        # column 6: dust emission scattered flux; lambda*F_lambda (W/m2)
        # column 7: transparent flux; lambda*F_lambda (W/m2)

        from ..preparation import unitconversion

        # Open the SED table
        #sed.table = tables.from_file(path, format="ascii.no_header") # sometimes doesn't work ?? why ??
        #sed.table.rename_column("col1", "Wavelength")
        #sed.table.rename_column("col2", "Flux")

        if contribution == "total": columns = (0,1)
        elif contribution == "direct": columns = (0,2)
        elif contribution == "scattered": columns = (0,3)
        elif contribution == "dust": columns = (0,4)
        elif contribution == "dustscattered": columns = (0,5)
        elif contribution == "transparent": columns = (0,6)
        else: raise ValueError("Wrong value for 'contribution': should be 'total', 'direct', 'scattered', 'dust', 'dustscattered' or 'transparent'")

        wavelength_column, flux_column = np.loadtxt(path, dtype=float, unpack=True, skiprows=skiprows, usecols=columns)

        sed.table = tables.new([wavelength_column, flux_column], ["Wavelength", "Flux"])
        sed.table["Wavelength"].unit = Unit("micron")

        jansky_column = []

        for i in range(len(sed.table)):

            # Get the flux density in W / m2 and the wavelength in micron
            neutral_fluxdensity = sed.table["Flux"][i] * Unit("W/m2")
            wavelength = sed.table["Wavelength"][i] * Unit("micron")

            # Convert to Jansky (2 methods give same result)
            #jansky_ = unitconversion.neutral_fluxdensity_to_jansky(neutral_fluxdensity, wavelength)
            jansky = (neutral_fluxdensity / wavelength.to("Hz", equivalencies=spectral())).to("Jy").value

            # Add the fluxdensity in Jansky to the new column
            jansky_column.append(jansky)

        # Add the flux column in Jansky
        sed.table.remove_column("Flux")
        sed.table["Flux"] = jansky_column
        sed.table["Flux"].unit = "Jy"

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path
        :return:
        """

        # Write the SED table to file
        tables.write(self.table, path, format="ascii.ecsv")

# -----------------------------------------------------------------
