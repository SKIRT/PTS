#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.table Contains the SkirtTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..basics.table import SmartTable, initialize_table
from ..units.unit import parse_unit
from ..tools.utils import memoize_method, lazyproperty
from ..tools import numbers
from .wavelengthgrid import WavelengthGrid
from ..basics.log import log
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class TruncatedSKIRTTableError(Exception):

    """
    This class ...
    """

    def __init__(self, message, path=None):

        """
        Thisf unction ...
        :param message:
        :param path:
        """

        # Call the base class constructor with the parameters it needs
        super(TruncatedSKIRTTableError, self).__init__(message)

        # The FITS file path
        self.path = path

# -----------------------------------------------------------------

def is_valid(path):

    """
    This function ...
    :param path:
    :return:
    """

    columns = np.loadtxt(path, unpack=True, ndmin=2)

    #print(columns)
    #print(len(columns))

    number_of_columns = len(columns)

    # THIS WAS BEFORE I DISCOVERED THE NDMIN PARAMETER
    # # Try to interpret the number of rows
    # try: number_of_rows = len(columns[0])
    # except TypeError:  # object of type 'numpy.float64' has no len()
    #     #raise TruncatedSKIRTTableError("The file only contains one line", path=path)
    #     return False

    number_of_rows = len(columns[0])

    #print("Number of columns: " + str(number_of_columns))
    #print("Number of rows: " + str(number_of_rows))

    # ONLY ONE ROW: NOT NORMAL
    if number_of_rows == 1:
        #raise TruncatedSKIRTTableError("The file only contains one line", path=path)
        return False

    # Get first ncolumn+1 rows of the file = the header
    header = fs.get_first_lines(path, number_of_columns + 1)

    # Loop over the lines
    has_title = False
    for index, line in enumerate(header):

        # Set flags
        first = index == 0
        last = index == number_of_columns

        if last and not has_title: break # one row less

        # Should start as comment
        if not line.startswith("#"): return False

        # Should state column index, no not always!
        if not line.startswith("# column "):
            if not first: return False
            else: has_title = True

    # Every check passed
    return True

# -----------------------------------------------------------------

class SkirtTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SkirtTable, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, expected_nrows=None, method="numpy", set_masks=True):

        """
        This function ...
        :param path:
        :param expected_nrows:
        :param method:
        :param set_masks:
        :return:
        """

        # Get the data
        data, names, units = get_skirt_data(path, expected_nrows=expected_nrows, method=method)
        #print(names, units)

        # Debugging
        log.debug("Creating the table ...")

        # Create a new table from the data
        table = cls(data=data, names=names, masked=True, copy=False)
        #print(table.colnames)

        # SET THE DATA
        # Set mask for each column from None values
        # Are there even None values? Data are Numpy arrays!
        #for column_index in range(len(names)): table[names[column_index]].mask = [value is None for value in data[column_index]]
        if set_masks:
            # NEW, not even necessary?
            # -> introduced as flag 'set_masks' to be able to turn this off
            for column_index in range(len(names)): table[names[column_index]].mask = np.isnan(data[column_index])

        # Set the column units
        for column_name in units: table[column_name].unit = parse_unit(units[column_name])

        # Initialize
        initialize_table(table, "SKIRT")

        # Return the table
        return table

    # -----------------------------------------------------------------

    def get_array(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        return np.asarray(self[column_name])

    # -----------------------------------------------------------------

    @memoize_method
    def get_mean(self, column_name, weights=None, weights_column_name=None):

        """
        This function ...
        :param column_name:
        :param weights:
        :param weights_column_name:
        :return:
        """

        # Check
        if weights is not None and weights_column_name is not None: raise ValueError("Cannot specify both weights as weights_column_name")

        # Weights from column
        if weights_column_name is not None: weights = self.get_array(weights_column_name)

        # Get the values
        values = self.get_array(column_name)

        # Weighed?
        if weights is not None: return numbers.weighed_arithmetic_mean(values, weights)

        # Not weighed
        else: return np.mean(values)

    # -----------------------------------------------------------------

    @memoize_method
    def get_median(self, column_name):

        """
        This function ...
        :param column_name:
        :return:
        """

        # Get the values
        values = self.get_array(column_name)

        # Return the median
        return np.nanmedian(values)

    # -----------------------------------------------------------------

    @memoize_method
    def get_stddev(self, column_name, weights=None, weights_column_name=None):

        """
        This function ...
        :param column_name:
        :param weights:
        :param weights_column_name:
        :return:
        """

        # Check
        if weights is not None and weights_column_name is not None: raise ValueError("Cannot specify both weights as weights_column_name")

        # Weights from column
        if weights_column_name is not None: weights = self.get_array(weights_column_name)

        # Get the values
        values = self.get_array(column_name)

        # Weighted?
        if weights is not None: return numbers.weighed_standard_deviation(values, weights)

        # Not weighed
        else: return np.nanstd(values)

# -----------------------------------------------------------------

class AbsorptionSpectraTable(SkirtTable):
    
    """
    This class ...
    """

    @classmethod
    def from_file(cls, path, ncells=None):

        """
        This function ...
        :param path:
        :param ncells: if the number of dust cells is known, pass this to check with the number of rows in the file (to check whether it was not truncated)
        :return:
        """

        # Get the data
        data, full_names, units = get_skirt_data(path, expected_nrows=ncells)

        # Set new column names
        nwavelengths = len(full_names)
        names = ["Luminosity" + str(index) for index in range(nwavelengths)]
        new_names = dict(zip(full_names, names)) # mapping from old names to new names

        # Create a new table from the data
        table = cls(data=data, names=names, masked=True)

        # Parse the wavelengths
        from ..units.parsing import parse_quantity
        wavelengths = [parse_quantity(name.split("lambda = ")[1]) for name in full_names]
        table.meta["wavelengths"] = wavelengths

        # SET THE DATA
        # Set mask for each column from None values
        for column_index in range(len(names)): table[names[column_index]].mask = [value is None for value in data[column_index]]

        # Set the column units
        for column_name in units: table[new_names[column_name]].unit = parse_unit(units[column_name])

        # Initialize
        initialize_table(table, "Absorption spectra")

        # Return the table
        return table

    # -----------------------------------------------------------------

    @property
    def ncells(self):
        return self.nrows

    # -----------------------------------------------------------------

    @property
    def wavelengths(self):
        return self.meta["wavelengths"]

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths(self):
        return len(self.wavelengths)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):
        return WavelengthGrid.from_wavelengths(self.wavelengths)

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):
        return self.wavelengths[0]

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):
        return self.wavelengths[-1]

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):
        return self.wavelength_grid.unit

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_array(self):
        return np.asarray(self.wavelength_grid.table["Wavelength"])

    # -----------------------------------------------------------------

    @property
    def _column_units(self):
        return [self.get_column_unit(name) for name in self.colnames]

    # -----------------------------------------------------------------

    @lazyproperty
    def luminosity_unit(self):
        from ..tools import sequences
        return sequences.get_all_equal_value(self._column_units)

    # -----------------------------------------------------------------

    def get_luminosities(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_row(index, as_list=True)

    # -----------------------------------------------------------------

    def get_luminosity_array(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return np.array(self.get_row(index, as_list=True, unit=self.luminosity_unit, add_unit=False))

    # -----------------------------------------------------------------

    def get_sed(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get the values
        luminosities = self.get_luminosity_array(index)

        from ..data.sed import SED
        return SED.from_arrays(self.wavelength_array, luminosities, self.wavelength_unit, self.luminosity_unit)

    # -----------------------------------------------------------------

    def get_luminosities_for_wavelength_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        colname = "Luminosity" + str(index)
        return self.get_column_values(colname, unit=self.luminosity_unit)

    # -----------------------------------------------------------------

    def get_luminosity_array_for_wavelength_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        colname = "Luminosity" + str(index)
        return self.get_column_array(colname, unit=self.luminosity_unit)

    # -----------------------------------------------------------------

    def get_average_luminosity_for_wavelength_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return np.nanmean(self.get_luminosity_array_for_wavelength_index(index))

    # -----------------------------------------------------------------

    @lazyproperty
    def average_sed(self):

        """
        This function ...
        :return:
        """

        # Calculate average luminosities
        luminosities = [self.get_average_luminosity_for_wavelength_index(index) for index in range(self.nwavelengths)]

        # Create and return
        from ..data.sed import SED
        return SED.from_arrays(self.wavelength_array, luminosities, self.wavelength_unit, self.luminosity_unit)

# -----------------------------------------------------------------

def get_skirt_data(path, expected_nrows=None, method="numpy"):

    """
    This function ...
    :param path:
    :param expected_nrows:
    :param method:
    :return:
    """

    # Debugging
    log.debug("Loading SKIRT table from '" + path + "' using " + method.capitalize() + " ...")

    # Get the columns
    columns = fs.get_columns(path, method=method)

    # Get number of columns and number of rows
    number_of_columns = len(columns)
    number_of_rows = len(columns[0])

    # THIS WAS BEFORE I DISCOVERED THE NDMIN PARAMETER
    # try: number_of_rows = len(columns[0])
    # except TypeError: # object of type 'numpy.float64' has no len()
    #     raise TruncatedSKIRTTableError("The file only contains one line", path=path)

    if number_of_rows == 1: raise TruncatedSKIRTTableError("The file only contains one line", path=path)

    # Check expected number of rows
    if expected_nrows is not None and number_of_rows != expected_nrows:
        raise IOError("Expected " + str(expected_nrows) + " rows but only found " + str(number_of_rows))

    data = []
    names = []
    units = dict()

    # Debugging
    log.debug("Reading header ...")

    # Get names and units
    with open(path) as table_file:

        for i in range(number_of_columns):

            line = table_file.next()
            column_number = i + 1

            if "column " + str(column_number) not in line:

                if i != 0:
                    raise IOError("Column name and unit of column " + str(i + 1) + " not found in file header")
                else:

                    # Determine offset
                    while "column" not in line: line = table_file.next()

            name_and_unit = line.split(": ")[1].split("\n")[0]
            if "(" in name_and_unit and ")" in name_and_unit:
                name = name_and_unit.split(" (")[0].capitalize() + name_and_unit.split(" (")[1].split(")")[1]
                unit = name_and_unit.split(" (")[1].split(")")[0]
                if "dimensionless" in unit: unit = None
            else:
                name = name_and_unit.capitalize()
                unit = None

            if ", i.e." in name: name = name.split(", i.e.")[0]

            # print(name)

            data.append(columns[i])
            names.append(name)

            if unit is not None: units[name] = unit

    # Return
    return data, names, units

# -----------------------------------------------------------------
