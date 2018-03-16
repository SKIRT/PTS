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
from astropy.io.ascii.core import InconsistentTableError

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

class NotRealColumn(Exception):

    """
    This class ...
    """

    def __init__(self, message, column_name):

        """
        Thisf unction ...
        :param message:
        :param column_name:
        """

        # Call the base class constructor with the parameters it needs
        super(NotRealColumn, self).__init__(message)

        # The column name
        self.column_name = column_name

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

def load_sed(path, wavelength_unit=None, photometry_unit=None):

    """
    This function ...
    :param path:
    :param wavelength_unit:
    :param photometry_unit:
    :return:
    """

    if is_from_skirt(path): return SED.from_skirt(path)
    else:

        # Try PTS (or at least ECSV) format
        try:

            # Try loading the table as ECSV format
            table = tables.from_file(path)

            # Observed SED or SED?
            if "Observatory" in table.colnames and "Instrument" in table.colnames: return ObservedSED.from_file(path)
            else: return SED.from_file(path)

        # Not readable as ECSV table, try reading as data file (wavelength and photometry columns)
        except InconsistentTableError:

            # Try loading as regular SED
            try: return SED.from_text_file(path, wavelength_unit=wavelength_unit, photometry_unit=photometry_unit)

            # ValueError: could not convert string to float: try as observed SED with string column
            except NotRealColumn: return ObservedSED.from_text_file(path, wavelength_unit=wavelength_unit, photometry_unit=photometry_unit)

# -----------------------------------------------------------------

def load_multiple_seds(path, wavelength_unit=None, photometry_unit=None, as_dict=False):

    """
    This function ...
    :param path:
    :param wavelength_unit:
    :param photometry_unit:
    :param as_dict:
    :return:
    """

    from pts.core.tools import filesystem as fs

    seds = []
    names = []

    # From SKIRT
    if is_from_skirt(path):

        ncols = get_ncolumns(path)
        if ncols == 2:
            sed = SED.from_skirt(path)
            seds.append(sed)
            names.append('total')
        else:
            contributions = ['total', 'direct', 'scattered', 'dust', 'dustscattered', 'transparent']
            for contribution in contributions:
                sed = SED.from_skirt(path, contribution=contribution)
                seds.append(sed)
                names.append(contribution)

    # Not from SKIRT
    else:

        # Try PTS (or at least ECSV) format
        try:

            # Try loading the table as ECSV format
            table = tables.from_file(path)

            # Observed SED or SED?
            if "Observatory" in table.colnames and "Instrument" in table.colnames: sed = ObservedSED.from_file(path)
            else: sed = SED.from_file(path)

            seds.append(sed)
            names.append(fs.strip_extension(fs.name(path)))

        # Not readable as ECSV table, try reading as data file (wavelength and photometry columns)
        except InconsistentTableError:

            try:

                ncols = get_ncolumns(path)
                nphotometry_cols = ncols - 1
                column_names = fs.get_header_labels(path)

                for index in range(nphotometry_cols):

                    name = column_names[index+1]
                    sed = SED.from_text_file(path, wavelength_unit=wavelength_unit, photometry_unit=photometry_unit, photometry_column=index+1)
                    seds.append(sed)
                    names.append(name)

            # There is a string column
            except NotRealColumn:

                #print(path)
                sed = ObservedSED.from_text_file(path, wavelength_unit=wavelength_unit, photometry_unit=photometry_unit)
                seds.append(sed)
                names.append(fs.strip_extension(fs.name(path)))

    # Return the SEDs
    if as_dict: return dict(zip(names, seds))
    else: return seds

# -----------------------------------------------------------------

def get_ncolumns(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..tools import filesystem as fs

    # Load last 2 lines
    lines = fs.tail(path, 2)

    # Interpret
    return np.loadtxt(lines, dtype='str').shape[1]

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
    def from_text_file(cls, path, wavelength_unit=None, photometry_unit=None, density=False, wavelength_column=0, photometry_column=1):

        """
        This function ...
        :param path:
        :param wavelength_unit:
        :param photometry_unit:
        :param density:
        :param wavelength_column:
        :param photometry_column:
        :return:
        """

        # OLD IMPLEMENTATION WITH NUMPY
        # # Check whether units are passed
        # if wavelength_unit is None: raise ValueError("Wavelength unit has to be specified")
        # if photometry_unit is None: raise ValueError("Photometry unit has to be specified")
        #
        # # Determine the columns to use
        # columns = (0, photometry_column)
        #
        # # Load the data
        # wavelength_column, photometry_column = np.loadtxt(path, dtype=float, unpack=True, skiprows=skiprows, usecols=columns)
        #
        # # Create the SED
        # return cls.from_arrays(wavelength_column, photometry_column, wavelength_unit, photometry_unit, density=density)

        from ..tools import types

        # Load as table
        table = tables.from_file(path, format="ascii")

        # Find wavelength column
        wavelength_column_name = table.colnames[wavelength_column]
        if not types.is_real_column(table[wavelength_column_name]): raise NotRealColumn("Wavelength column a re not real numbers", wavelength_column_name)

        # Find photometry column
        photometry_column_name = table.colnames[photometry_column]
        if not types.is_real_column(table[photometry_column_name]): raise NotRealColumn("Photometry column are not real numbers", photometry_column_name)

        # Get wavelength unit
        if table[wavelength_column_name].unit is not None:
            wavelength_unit_table = table[wavelength_column_name].unit
            if wavelength_unit is not None:
                if wavelength_unit != wavelength_unit_table: raise ValueError("Wavelength unit does not correspond that found in file")
            else: wavelength_unit = wavelength_unit_table
        else:
            if wavelength_unit is None: raise ValueError("Wavelength unit must be specified")

        # Get photometry unit
        if table[photometry_column_name].unit is not None:
            photometry_unit_table = PhotometricUnit(table[photometry_column_name].unit)
            if photometry_unit is not None:
                photometry_unit = PhotometricUnit(photometry_unit, density=density)
                if photometry_unit != photometry_unit_table: raise ValueError("Photometric unit does not correspond that found in file")
            else: photometry_unit = photometry_unit_table
        else:
            if photometry_unit is None: raise ValueError("Photometry unit must be specified")
            photometry_unit = PhotometricUnit(photometry_unit, density=density)

        # Create sed
        wavelengths = table[wavelength_column_name]
        photometry = table[photometry_column_name]
        sed = cls.from_arrays(wavelengths, photometry, wavelength_unit=wavelength_unit, photometry_unit=photometry_unit)

        # Return the SED
        return sed

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

    def photometry(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False,
                   min_wavelength=None, max_wavelength=None):

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

    def photometry_at(self, wavelength, unit=None, add_unit=True, density=False, brightness=False, interpolate=True, conversion_info=None):

        """
        This function ...
        :param wavelength:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :param interpolate:
        :param conversion_info:
        :return:
        """

        return self.value_for_wavelength(wavelength, unit=unit, add_unit=add_unit, density=density, brightness=brightness, interpolate=interpolate, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def normalize_at_wavelength(self, wavelength, value=1.):

        """
        This function ...
        :param wavelength:
        :param value:
        :return:
        """

        value_for_wavelength = self.photometry_at(wavelength, unit=self.unit, add_unit=False)
        self[self.y_name] /= value_for_wavelength * value
        self[self.y_name].unit = None

    # -----------------------------------------------------------------

    def normalized_photometry(self, value=1.0, method="integral", asarray=False):

        """
        This function ...
        :param value:
        :param method:
        :param asarray:
        :return:
        """

        # Max method
        if method == "max":

            # Get array
            photometry = self.photometry(self.unit, asarray=True)
            max_value = np.max(photometry)
            factor = value / max_value

            # Return the normalized data
            if asarray: return photometry * factor
            else: return list(photometry * factor)

        # Integral method
        elif method == "integral": raise NotImplementedError("Not implemented yet")

        # Invalid method
        else: raise ValueError("Invalid option for 'method'")

    # -----------------------------------------------------------------

    def normalized_photometry_at_wavelength(self, wavelength, value=1.):

        """
        This function ...
        :param wavelength:
        :param value:
        :return:
        """

        value_for_wavelength = self.photometry_at(wavelength, unit=self.unit, add_unit=False)

        # Get array
        photometry = self.photometry(self.unit, asarray=True)
        factor = value / value_for_wavelength

        # Return the normalized data
        return photometry * factor

    # -----------------------------------------------------------------

    @classmethod
    def from_arrays(cls, wavelengths, photometry, wavelength_unit, photometry_unit, density=False, brightness=False, distance=None):

        """
        This function ...
        :param wavelengths:
        :param photometry:
        :param wavelength_unit:
        :param photometry_unit:
        :param density:
        :param brightness:
        :param distance:
        :return:
        """

        # Set unit
        photometry_unit = PhotometricUnit(photometry_unit, density=density, brightness=brightness)

        # Create new SED
        sed = cls(photometry_unit=photometry_unit, density=density, brightness=brightness, distance=distance)

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
    def from_skirt(cls, path, skiprows=0, contribution="total", unit=None, remote=None, distance=None):

        """
        This function ...
        :param path:
        :param skiprows:
        :param contribution:
        :param unit: define the unit for the photometry for the SED
        :param remote:
        :param distance:
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
        sed = cls(photometry_unit=unit, distance=distance)

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

    @classmethod
    def from_text_file(cls, path, photometry_unit=None, density=False, photometry_column=None, wavelength_unit="micron"):

        """
        Thisn function ...
        :param path:
        :param photometry_unit:
        :param density:
        :param photometry_column:
        :param wavelength_unit:
        :return:
        """

        from ..tools import types

        # Load as table
        table = tables.from_file(path, format="ascii")

        # Find string column for the filters
        filter_column_name = None
        for name in table.colnames:
            if types.is_string_column(table[name]):
                filter_column_name = name
                break
        if filter_column_name is None: raise ValueError("Cannot find filter column")

        # Find photometry column
        if photometry_column is not None: photometry_column_name = table.colnames[photometry_column]
        else:

            photometry_column_name = None
            first_real_column_name = None
            for name in table.colnames:
                if name == filter_column_name: continue
                if not types.is_real_column(table[name]): continue
                if first_real_column_name is None: first_real_column_name = name
                if "flux" in name.lower() or "lum" in name.lower():
                    photometry_column_name = name
                    break
            if photometry_column_name is None and first_real_column_name is not None: photometry_column_name = first_real_column_name
            if photometry_column_name is None: raise ValueError("Cannot find photometry column")

        # Find error column
        error_column_name = None
        for name in table.colnames:
            if name == filter_column_name: continue
            if name == photometry_column_name: continue
            if "err" in name.lower():
                error_column_name = name
                break

        # Get photometry unit
        if table[photometry_column_name].unit is not None:
            photometry_unit_table = PhotometricUnit(table[photometry_column_name].unit)
            if photometry_unit is not None:
                photometry_unit = PhotometricUnit(photometry_unit, density=density)
                if photometry_unit != photometry_unit_table: raise ValueError("Photometric unit does not correspond that found in file")
            else: photometry_unit = photometry_unit_table
        else:
            if photometry_unit is None: raise ValueError("Photometry unit must be specified")
            photometry_unit = PhotometricUnit(photometry_unit, density=density)

        # Get error unit
        if error_column_name is not None:
            if table[error_column_name].unit is not None:
                error_unit = PhotometricUnit(table[error_column_name].unit)
            else:
                if photometry_unit is None: raise ValueError("Photometry unit must be specified")
                error_unit = PhotometricUnit(photometry_unit, density=density)
        else: error_unit = None

        # Create sed
        sed = cls(wavelength_unit=wavelength_unit, photometry_unit=photometry_unit)
        sed._setup()

        # Add points
        nrows = len(table)
        for index in range(nrows):

            filter_name = table[filter_column_name][index]
            photometry = table[photometry_column_name][index]
            error = table[error_column_name][index] if error_column_name is not None else None
            if np.isnan(error): error = None

            # Add units
            photometry = photometry * photometry_unit
            if error is not None: error = error * error_unit

            # Parse filter
            fltr = parse_filter(filter_name)

            # Add point
            sed.add_point(fltr, photometry, error=error, sort=False)

        # Sort the SED on wavelength
        sed.sort("Wavelength")

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Multiply column with value
        self[self.value_name] *= value
        return self

    # -----------------------------------------------------------------

    def add_point(self, fltr, photometry, error=None, conversion_info=None, sort=True):

        """
        This function ...
        :param fltr:
        :param photometry:
        :param error:
        :param conversion_info:
        :param sort:
        :return:
        """

        # Create fltr
        fltr = parse_filter(fltr)

        if error is not None:
            if not isinstance(error, ErrorBar): error = ErrorBar(error)
            values = [fltr.observatory, fltr.instrument, fltr.band, fltr.pivot, photometry, error.lower, error.upper]
        else: values = [fltr.observatory, fltr.instrument, fltr.band, fltr.pivot, photometry, None, None]

        # Set conversion info
        if conversion_info is None:
            conversion_info_photometry = dict()
            conversion_info_photometry["wavelength"] = fltr.pivot
            conversion_info = {self.value_name: conversion_info_photometry}

        # Add the row
        self.add_row(values, conversion_info=conversion_info)

        # Sort
        if sort: self.sort("Wavelength")

    # -----------------------------------------------------------------

    def photometry(self, unit=None, asarray=False, add_unit=True, density=False, brightness=False, conversion_info=None):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param density:
        :param brightness:
        :param conversion_info:
        :return:
        """

        return self.values(unit, asarray, add_unit, density=density, brightness=brightness, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def photometry_at(self, wavelength, unit=None, add_unit=True, density=False, brightness=False, conversion_info=None):

        """
        This function ...
        :param wavelength:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :param conversion_info:
        :return:
        """

        return self.value_for_wavelength(wavelength, unit=unit, add_unit=add_unit, density=density, brightness=brightness, conversion_info=conversion_info)

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

    def get_lower_error(self, index, add_unit=True, unit=None, density=False, brightness=False, conversion_info=None):

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

        # Get the value
        value = self.get_value("Error-", index, add_unit=True)
        if value is None: return None

        # Convert unit if necessary
        if unit is not None:

            # Parse the unit
            unit = u(unit, density=density, brightness=brightness)

            # Create conversion info
            if conversion_info is None: conversion_info = dict()
            conversion_info["wavelength"] = self.wavelength_for_index(index)
            if self.distance is not None: conversion_info["distance"] = self.distance

            # Create converted value
            value = value.to(unit, **conversion_info)

        # Remove unit if requested
        if not add_unit: value = value.value

        # Return the value
        return value

    # -----------------------------------------------------------------

    def get_upper_error(self, index, add_unit=True, unit=None, density=False, brightness=False, conversion_info=None):

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

        # Get the value
        value = self.get_value("Error+", index, add_unit=True)
        if value is None: return None

        # Convert unit if necessary
        if unit is not None:

            # Parse the unit
            unit = u(unit, density=density, brightness=brightness)

            # Create conversion info
            if conversion_info is None: conversion_info = dict()
            conversion_info["wavelength"] = self.wavelength_for_index(index)
            if self.distance is not None: conversion_info["distance"] = self.distance

            # Create converted value
            value = value.to(unit, **conversion_info)

        # Remove unit if requested
        if not add_unit: value = value.value

        # Return the value
        return value

    # -----------------------------------------------------------------

    def set_lower_error(self, index, value):

        """
        This function ...
        :param index:
        :param value:
        :return:
        """

        self.set_value("Error-", index, value)

    # -----------------------------------------------------------------

    def set_upper_error(self, index, value):

        """
        This function ...
        :param index:
        :param value:
        :return:
        """

        self.set_value("Error+", index, value)

    # -----------------------------------------------------------------

    def add_relative_error(self, relative_error):

        """
        This function ...
        :param relative_error:
        :return:
        """

        # Loop over the points in the SED
        for index in range(len(self)):

            # Calculate the additional error
            #wavelength = self.get_wavelength(index)
            value = self.get_photometry(index)
            additional_error = value * relative_error
            unit = value.unit
            additional_error_value = additional_error.to(unit).value

            # Set the lower error
            lower = self.get_lower_error(index)
            if lower is not None:

                # Add the additional error in quadrature
                lower_value = lower.to(unit).value
                lower = -np.sqrt(lower_value ** 2 + additional_error_value ** 2) * unit

                # Set the new error
                self.set_lower_error(index, lower)

            # Set the upper error
            upper = self.get_upper_error(index)
            if upper is not None:

                # Add the additional error in quadrature
                upper_value = upper.to(unit).value
                upper = np.sqrt(upper_value**2 + additional_error_value**2) * unit

                # Set the new error bar
                self.set_upper_error(index, upper)

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

    def errors(self, unit=None, add_unit=True, density=False, brightness=False, conversion_info=None):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :param conversion_info:
        :return:
        """

        column_units = [self.column_unit("Error-"), self.column_unit("Error+")]
        return tables.columns_as_objects([self["Error-"], self["Error+"]], ErrorBar, unit=unit, add_unit=add_unit, column_units=column_units, density=density, brightness=brightness, conversion_info=conversion_info)

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
            instrument_entry = self.get_value("Instrument", i)
            band_entry = self.get_value("Band", i)

            if not (instrument_entry == instrument and band_entry == band): continue

            if has_unit:

                # Add the unit initially to be able to convert
                #error_min = self["Error-"][i] * self["Error-"].unit
                #error_plus = self["Error+"][i] * self["Error+"].unit
                error_min = self.get_value("Error-", i)
                error_plus = self.get_value("Error+", i)

                # If a target unit is specified, convert
                if unit is not None:

                    error_min = error_min.to(unit).value * u(unit)
                    error_plus = error_plus.to(unit).value * u(unit)

                if not add_unit:

                    error_min = error_min.value
                    error_plus = error_plus.value

                error = ErrorBar(error_min, error_plus)

            else: error = ErrorBar(self.get_value("Error-", i), self.get_value("Error+", i))

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
