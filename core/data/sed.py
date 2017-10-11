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

# Import the relevant PTS classes and modules
from ..basics.curve import WavelengthCurve, FilterCurve
from ..units.unit import PhotometricUnit
from ..tools import tables
from ..filter.broad import BroadBandFilter
from ...magic.tools.colours import calculate_colour
from ...core.basics.errorbar import ErrorBar
from ..units.parsing import parse_unit as u
from ..filter.filter import parse_filter
from ..tools import arrays
from ..simulation import textfile

# -----------------------------------------------------------------

def is_from_skirt(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..tools import filesystem as fs

    # Read the two first lines
    first, second = fs.get_first_lines(path, nlines=2)

    # Check first
    if not first.startswith("# column 1: lambda"): return False

    # Check second
    if not second.startswith("# column 2:"): return False

    # Checks passed for SKIRT
    return True

# -----------------------------------------------------------------

def load_sed(path):

    """
    This function ...
    :param path:
    :return:
    """

    if is_from_skirt(path): return SED.from_skirt(path)
    else:

        table = tables.from_file(path)

        # Observed SED or SED?
        if "Observatory" in table.colnames and "Instrument" in table.colnames: return ObservedSED.from_file(path)
        else: return SED.from_file(path)

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
            density = kwargs.pop("density", False)
            brightness = kwargs.pop("brightness", False)
            unit = PhotometricUnit(unit, density=density, brightness=brightness)

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

    def get_photometry(self, index, add_unit=True, unit=None, density=False, brightness=False, conversion_info=None):

        """
        This function ...
        :param index:
        :param add_unit:
        :param unit:
        :param density:
        :param brightness:
        :param conversion_info:
        :return:
        """

        return self.photometry_for_index(index, add_unit=add_unit, unit=unit, density=density, brightness=brightness, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def photometry_for_index(self, index, add_unit=True, unit=None, density=False, brightness=False, conversion_info=None):

        """
        This function ...
        :param index:
        :param add_unit:
        :param unit:
        :param density:
        :param brightness:
        :param conversion_info:
        :return:
        """

        return self.value_for_index(index, add_unit=add_unit, unit=unit, density=density, brightness=brightness, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def photometry(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        return self.values(unit, asarray, add_unit, conversion_info=conversion_info, density=density, brightness=brightness,
                           min_wavelength=min_wavelength, max_wavelength=max_wavelength)

    # -----------------------------------------------------------------

    def photometry_at(self, wavelength, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param wavelength:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        return self.value_for_wavelength(wavelength, unit=unit, add_unit=add_unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    @classmethod
    def from_arrays(cls, wavelengths, photometry, wavelength_unit, photometry_unit, density=False, brightness=False):

        """
        This function ...
        :param wavelengths:
        :param photometry:
        :param wavelength_unit:
        :param photometry_unit:
        :param density:
        :param brightness:
        :return:
        """

        # Set unit
        photometry_unit = PhotometricUnit(photometry_unit, density=density, brightness=brightness)

        # Create new SED
        sed = cls(photometry_unit=photometry_unit, density=density, brightness=brightness)

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
    def from_skirt(cls, path, skiprows=0, contribution="total", unit=None, remote=None):

        """
        This function ...
        :param path:
        :param skiprows:
        :param contribution:
        :param unit: define the unit for the photometry for the SED
        :param remote:
        :return:
        """

        # SEDInstrument:
        # column 1: lambda (micron)
        # column 2: total flux; lambda*F_lambda (W/m2)
        # of:
        # column 1: lambda (micron)
        # column 2: total flux; F_nu (Jy)

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

        # Keep track of the units of the different columns
        units = textfile.get_units(path, remote=remote)

        # Define index of different columns
        contributions_index = dict()
        contributions_index["total"] = 1
        contributions_index["direct"] = 2
        contributions_index["scattered"] = 3
        contributions_index["dust"] = 4
        contributions_index["dustscattered"] = 5
        contributions_index["transparent"] = 6

        # Load the column data
        if contribution not in contributions_index: raise ValueError("Wrong value for 'contribution': should be 'total', 'direct', 'scattered', 'dust', 'dustscattered' or 'transparent'")
        columns = (0, contributions_index[contribution])

        # LOAD THE DATA
        if remote is not None:
            lines = remote.get_lines(path, add_sep=True)
            wavelength_column, photometry_column = np.loadtxt(lines, dtype=float, unpack=True, skiprows=skiprows, usecols=columns)
        else: wavelength_column, photometry_column = np.loadtxt(path, dtype=float, unpack=True, skiprows=skiprows, usecols=columns)

        # Get column units
        wavelength_unit = units[0]
        photometry_unit = units[contributions_index[contribution]]
        if unit is None: unit = photometry_unit

        # Create a new SED
        sed = cls(photometry_unit=unit)

        # Add the entries
        for index in range(len(wavelength_column)):

            # Get values
            wavelength = wavelength_column[index] * wavelength_unit
            photometry = photometry_column[index] * photometry_unit

            # Add point
            sed.add_point(wavelength, photometry)

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def convert_to(self, wavelength_unit=None, photometry_unit=None):

        """
        This function ...
        :param wavelength_unit:
        :param photometry_unit:
        :return:
        """

        # If wavelength unit has to be converted
        if wavelength_unit is not None:

            wavelength_unit = u(wavelength_unit)
            self["Wavelength"] = self.wavelengths(asarray=True, unit=wavelength_unit)
            self["Wavelength"].unit = wavelength_unit

        # If photometry unit has to be converted
        if photometry_unit is not None:

            photometry_unit = PhotometricUnit(photometry_unit)
            self["Photometry"] = self.photometry(asarray=True, unit=photometry_unit)
            self["Photometry"].unit = photometry_unit

            # Set whether this column is a spectral density
            if photometry_unit.density:
                if "density" not in self.meta: self.meta["density"] = []
                self.meta["density"].append("Photometry")

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
            self.add_column_info("Error-", float, unit, "Lower bound error")
            self.add_column_info("Error+", float, unit, "Upper bound error")

    # -----------------------------------------------------------------

    def add_point(self, fltr, photometry, error=None):

        """
        This function ...
        :param fltr:
        :param photometry:
        :param error:
        :return:
        """

        # Create fltr
        fltr = parse_filter(fltr)

        if error is not None:
            if not isinstance(error, ErrorBar): error = ErrorBar(error)
            values = [fltr.observatory, fltr.instrument, fltr.band, fltr.pivot, photometry, error.lower, error.upper]
        else: values = [fltr.observatory, fltr.instrument, fltr.band, fltr.pivot, photometry, None, None]
        self.add_row(values)

    # -----------------------------------------------------------------

    def photometry(self, unit=None, asarray=False, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        return self.values(unit, asarray, add_unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def photometry_at(self, wavelength, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param wavelength:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        return self.value_for_wavelength(wavelength, unit=unit, add_unit=add_unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def get_photometry(self, index, add_unit=True, unit=None, density=False, brightness=False, conversion_info=None):

        """
        This function ...
        :param index:
        :param add_unit:
        :param unit:
        :param density:
        :param brightness:
        :param conversion_info:
        :return:
        """

        return self.photometry_for_index(index, add_unit=add_unit, unit=unit, density=density, brightness=brightness, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def photometry_for_index(self, index, add_unit=True, unit=None, density=False, brightness=False, conversion_info=None):

        """
        This function ...
        :param index:
        :param add_unit:
        :param unit:
        :param density:
        :param brightness:
        :param conversion_info:
        :return:
        """

        return self.value_for_index(index, add_unit=add_unit, unit=unit, density=density, brightness=brightness, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def photometry_for_band(self, instrument, band, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param instrument:
        :param band:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        return self.value_for_band(instrument, band, unit, add_unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def photometry_for_filter(self, fltr, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        return self.value_for_filter(fltr, unit=unit, add_unit=add_unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def errors(self, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        column_units = [self.column_unit("Error-"), self.column_unit("Error+")]
        return tables.columns_as_objects([self["Error-"], self["Error+"]], ErrorBar, unit=unit, add_unit=add_unit, column_units=column_units, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def errors_min(self, unit=None, asarray=False, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        if asarray: return arrays.plain_array(self["Error-"], unit=unit, array_unit=self.column_unit("Error-"), density=density, brightness=brightness)
        else: return arrays.array_as_list(self["Error-"], unit=unit, add_unit=add_unit, array_unit=self.column_unit("Error-"), density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def errors_max(self, unit=None, asarray=False, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        if asarray: return arrays.plain_array(self["Error+"], unit=unit, array_unit=self.column_unit("Error+"), density=density, brightness=brightness)
        else: return arrays.array_as_list(self["Error+"], unit=unit, add_unit=add_unit, array_unit=self.column_unit("Error+"), density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def error_for_filter(self, fltr, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        return self.error_for_band(fltr.instrument, fltr.band, unit=unit, add_unit=add_unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def error_for_band(self, instrument, band, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param instrument:
        :param band:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        has_unit = self["Error-"].unit is not None and self["Error+"].unit is not None
        has_mask = hasattr(self["Error-"], "mask")
        assert has_mask == hasattr(self["Error+"], "mask")

        # If unit has to be converted, check whether the original unit is specified
        if not has_unit and unit is not None: raise ValueError("Cannot determine the unit of the error columns so values cannot be converted to " + str(unit))

        # Parse the unit
        if unit is not None: unit = u(unit, density=density, brightness=brightness)

        # Loop over all the entries in the table
        for i in range(len(self)):

            # Get instrument and band
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
            fltr = BroadBandFilter(instrument + " " + band)

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
        sed = cls(photometry_unit="Jy")

        # Add entries
        for index in range(len(filter_column)):
            sed.add_point(filter_column[index], flux_column[index], ErrorBar(fluxerrmin_column[index], fluxerrmax_column[index]))

        # Return the sed
        return sed

    # -----------------------------------------------------------------

    def convert_to(self, wavelength_unit=None, photometry_unit=None):

        """
        This function ...
        :param wavelength_unit:
        :param photometry_unit:
        :return:
        """

        # If wavelength unit has to be converted
        if wavelength_unit is not None:

            wavelength_unit = u(wavelength_unit)
            self["Wavelength"] = self.wavelengths(asarray=True, unit=wavelength_unit)
            self["Wavelength"].unit = wavelength_unit

        # If photometry unit has to be converted
        if photometry_unit is not None:

            photometry_unit = PhotometricUnit(photometry_unit)

            self["Photometry"] = self.photometry(asarray=True, unit=photometry_unit)
            self["Photometry"].unit = photometry_unit

            self["Error-"] = self.errors_min(asarray=True, unit=photometry_unit)
            self["Error-"].unit = photometry_unit

            self["Error+"] = self.errors_max(asarray=True, unit=photometry_unit)
            self["Error+"].unit = photometry_unit

            # Set whether this column is a spectral density
            if photometry_unit.density:
                if "density" not in self.meta: self.meta["density"] = []
                self.meta["density"].append("Photometry")

# -----------------------------------------------------------------
