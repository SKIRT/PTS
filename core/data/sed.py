#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.data.sed Contains the SED and ObservedSED class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import spectral

# Import the relevant PTS classes and modules
from ..basics.curve import WavelengthCurve, FilterCurve
from ..basics.unit import PhotometricUnit
from ..tools import tables
from ..basics.filter import Filter
from ...magic.tools.colours import calculate_colour
from ...core.basics.errorbar import ErrorBar
from ..basics.unit import parse_unit as u

# -----------------------------------------------------------------

class SED(WavelengthCurve):
    
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        if "photometry_unit" in kwargs: from_astropy = False
        else: from_astropy = True

        if not from_astropy:

            # Call the initialize function of the base class
            unit = kwargs.pop("photometry_unit")
            density = kwargs.pop("density")
            unit = PhotometricUnit(unit, density=density)

            kwargs["y_name"] = "Photometry"
            kwargs["y_description"] = "Photometric points"
            kwargs["y_unit"] = unit

        # Call the constructor of the base class
        super(SED, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_text_file(cls, path, wavelength_unit, photometry_unit, skiprows=None, density=False):

        """
        This function ...
        :param path:
        :param wavelength_unit:
        :param photometry_unit:
        :param skiprows:
        :param density:
        :return:
        """

        # Load the data
        wavelength_column, luminosity_column = np.loadtxt(path, dtype=float, unpack=True, skiprows=skiprows)

        # Create the SED
        return cls.from_arrays(wavelength_column, luminosity_column, wavelength_unit, photometry_unit, density=density)

    # -----------------------------------------------------------------

    def photometry(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        return self.values(unit, asarray, add_unit)

    # -----------------------------------------------------------------

    def photometry_at(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return self.value_for_wavelength(wavelength)

    # -----------------------------------------------------------------

    @classmethod
    def from_arrays(cls, wavelengths, photometry, wavelength_unit, photometry_unit, density=False):

        """
        This function ...
        :param wavelengths:
        :param photometry:
        :param wavelength_unit:
        :param photometry_unit:
        :param density:
        :return:
        """

        # Set unit
        photometry_unit = PhotometricUnit(photometry_unit, density=density)

        # Create new SED
        sed = cls(photometry_unit=photometry_unit, density=density)

        # Parse units
        wavelength_unit = u(wavelength_unit)
        photometry_unit = u(photometry_unit)

        # Add the entries
        for index in range(len(wavelengths)):

            # Get wavelength and measurement
            wavelength = wavelengths[index] * wavelength_unit
            phot = photometry[index] * photometry_unit

            # Add
            sed.add_point(wavelength, phot)

        # Return the sed
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

        # Open the SED table
        # sed.table = tables.from_file(path, format="ascii.no_header") # sometimes doesn't work ?? why ??
        # sed.table.rename_column("col1", "Wavelength")
        # sed.table.rename_column("col2", "Flux")

        if contribution == "total": columns = (0, 1)
        elif contribution == "direct": columns = (0, 2)
        elif contribution == "scattered": columns = (0, 3)
        elif contribution == "dust": columns = (0, 4)
        elif contribution == "dustscattered": columns = (0, 5)
        elif contribution == "transparent": columns = (0, 6)
        else: raise ValueError("Wrong value for 'contribution': should be 'total', 'direct', 'scattered', 'dust', 'dustscattered' or 'transparent'")

        wavelength_column, flux_column = np.loadtxt(path, dtype=float, unpack=True, skiprows=skiprows, usecols=columns)

        #sed.table = tables.new([wavelength_column, flux_column], ["Wavelength", "Flux"])
        #sed.table["Wavelength"].unit = Unit("micron")

        jansky_column = []

        for i in range(len(wavelength_column)):

            # Get the flux density in W / m2 and the wavelength in micron
            neutral_fluxdensity = flux_column[i] * PhotometricUnit("W/m2")
            wavelength = wavelength_column[i] * u("micron")

            # Convert to Jansky (2 methods give same result)
            # jansky_ = unitconversion.neutral_fluxdensity_to_jansky(neutral_fluxdensity, wavelength)
            jansky = (neutral_fluxdensity / wavelength.to("Hz", equivalencies=spectral())).to("Jy").value

            # Add the fluxdensity in Jansky to the new column
            jansky_column.append(jansky)

        # Add the flux column in Jansky
        #sed.table.remove_column("Flux")
        #sed.table["Flux"] = jansky_column
        #sed.table["Flux"].unit = "Jy"

        # Create a new SED
        sed = cls("Jy")

        # Add the entries
        for index in range(len(wavelength_column)):

            # Get values
            wavelength = wavelength_column[index] * u("micron")
            flux = jansky_column * PhotometricUnit("Jy")

            # Add point
            sed.add_point(wavelength, flux)

        # Return the SED
        return sed

# -----------------------------------------------------------------

class ObservedSED(FilterCurve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        if "photometry_unit" not in kwargs: from_astropy = True
        else: from_astropy = False

        if not from_astropy:

            # Get properties
            unit = kwargs.pop("photometry_unit")
            density = kwargs.pop("density", False)
            unit = PhotometricUnit(unit, density=density)

            # Set kwargs
            kwargs["y_name"] = "Photometry"
            kwargs["y_description"] = "Photometric points"
            kwargs["y_unit"] = unit

        # Call the constructor of the base class
        super(ObservedSED, self).__init__(*args, **kwargs)

        if not from_astropy:

            # Add columns
            self.column_info.append(("Error-", float, unit, "Lower bound error"))
            self.column_info.append(("Error+", float, unit, "Upper bound error"))

    # -----------------------------------------------------------------

    def add_point(self, fltr, photometry, error=None):

        """
        This function ...
        :param fltr:
        :param photometry:
        :param error:
        :return:
        """

        if error is not None:
            if not isinstance(error, ErrorBar): error = ErrorBar(error)
            values = [fltr.observatory, fltr.instrument, fltr.band, fltr.pivotwavelength(), photometry, error.lower, error.upper]
        else: values = [fltr.observatory, fltr.instrument, fltr.band, fltr.pivotwavelength(), photometry, None, None]
        self.add_row(values)

    # -----------------------------------------------------------------

    def photometry(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        return self.values(unit, asarray, add_unit)

    # -----------------------------------------------------------------

    def photometry_at(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return self.value_for_wavelength(wavelength)

    # -----------------------------------------------------------------

    def photometry_for_band(self, instrument, band, unit=None, add_unit=True):

        """
        This function ...
        :param instrument:
        :param band:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.value_for_band(instrument, band, unit, add_unit)

    # -----------------------------------------------------------------

    def photometry_for_filter(self, fltr, unit=None, add_unit=True):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.value_for_filter(fltr, unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def errors(self, unit=None, add_unit=True):

        """
        This function ...
        :param unit:
        :param add_unit:
        :return:
        """

        return tables.columns_as_objects([self["Error-"], self["Error+"]], ErrorBar, unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def errors_min(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self["Error-"], unit=unit)
        else: return tables.column_as_list(self["Error-"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def errors_max(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self["Error+"], unit=unit)
        else: return tables.column_as_list(self["Error+"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def error_for_filter(self, fltr, unit=None, add_unit=True):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :return:
        """

        return self.error_for_band(fltr.instrument, fltr.band, unit, add_unit)

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

        has_unit = self["Error-"].unit is not None and self["Error+"].unit is not None
        has_mask = hasattr(self["Error-"], "mask")
        assert has_mask == hasattr(self["Error+"], "mask")

        # If unit has to be converted, check whether the original unit is specified
        if not has_unit and unit is not None: raise ValueError(
            "Cannot determine the unit of the error columns so values cannot be converted to " + str(unit))

        # Loop over all the entries in the table
        for i in range(len(self)):

            instrument_entry = self["Instrument"][i]
            band_entry = self["Band"][i]

            if not (instrument_entry == instrument and band_entry == band): continue

            if has_unit:

                # Add the unit initially to be able to convert
                error_min = self["Error-"][i] * self["Error-"].unit
                error_plus = self["Error+"][i] * self["Error+"].unit

                # If a target unit is specified, convert
                if unit is not None:

                    error_min = error_min.to(unit).value * u(unit)
                    error_plus = error_plus.to(unit).value * u(unit)

                if not add_unit:

                    error_min = error_min.value
                    error_plus = error_plus.value

                error = ErrorBar(error_min, error_plus)

            else: error = ErrorBar(self["Error-"][i], self["Error+"][i])

            return error

        # If no match is found, return None
        return None

    # -----------------------------------------------------------------

    def colour(self, filter_a, filter_b):

        """
        This function ...
        :param filter_a:
        :param filter_b:
        :return:
        """

        return calculate_colour(self.photometry_for_filter(filter_a), self.photometry_for_filter(filter_b))

    # -----------------------------------------------------------------

    @classmethod
    def from_caapr(cls, path):

        """
        This function ...
        :return:
        """

        # Create a new observed SED
        sed = cls()

        # Load the table
        caapr_table = tables.from_file(path, format="csv")

        fluxes = dict()
        errors = dict()

        # Loop over the columns of the table
        for colname in caapr_table.colnames:

            if colname == "name": continue
            if "ERR" in colname:
                instrument_band = colname.split("_ERR")[0]
                error = abs(caapr_table[colname][0])
                if not np.isnan(error): errors[instrument_band] = error
            else:
                instrument_band = colname
                flux = caapr_table[colname][0]
                if not np.isnan(flux): fluxes[instrument_band] = flux

        filter_column = []
        observatory_column = []
        instrument_column = []
        band_column = []
        wavelength_column = []
        flux_column = []
        fluxerrmin_column = []
        fluxerrmax_column = []

        for instrument_band in fluxes:

            if not instrument_band in errors: raise ValueError("No error for " + instrument_band)

            flux = fluxes[instrument_band]
            error = errors[instrument_band]

            instrument = instrument_band.split("_")[0]
            band = instrument_band.split("_")[1]

            # Create filter
            fltr = Filter(instrument + " " + band)

            # Get filter properties
            observatory = fltr.observatory
            instrument = fltr.instrument
            band = fltr.band
            wavelength = fltr.pivotwavelength()

            # Add entry to the columns
            filter_column.append(fltr)
            observatory_column.append(observatory)
            instrument_column.append(instrument)
            band_column.append(band)
            wavelength_column.append(wavelength)
            flux_column.append(flux)
            fluxerrmin_column.append(-error)
            fluxerrmax_column.append(error)

        # Create the SED table
        #data = [observatory_column, instrument_column, band_column, wavelength_column, flux_column, fluxerrmin_column, fluxerrmax_column]
        #names = ["Observatory", "Instrument", "Band", "Wavelength", "Flux", "Error-", "Error+"]
        #sed.table = tables.new(data, names)
        #sed.table["Wavelength"].unit = "micron"
        #sed.table["Flux"].unit = "Jy"
        #sed.table["Error-"].unit = "Jy"
        #sed.table["Error+"].unit = "Jy"

        # Return the SED
        #return sed

        # Initialize SED
        sed = cls("Jy")

        # Add entries
        for index in range(len(filter_column)):
            sed.add_point(filter_column[index], flux_column[index], ErrorBar(fluxerrmin_column[index], fluxerrmax_column[index]))

        # Return the sed
        return sed

# -----------------------------------------------------------------
