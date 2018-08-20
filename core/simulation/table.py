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

    from ..tools import filesystem as fs

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
    def from_file(cls, path, expected_nrows=None, method="numpy"):

        """
        This function ...
        :param path:
        :param expected_nrows:
        :param method:
        :return:
        """

        # Get the data
        data, names, units = get_skirt_data(path, expected_nrows=expected_nrows, method=method)

        # Debugging
        log.debug("Creating the table ...")

        # Create a new table from the data
        table = cls(data=data, names=names, masked=True)

        # SET THE DATA
        # Set mask for each column from None values
        for column_index in range(len(names)):
            table[names[column_index]].mask = [value is None for value in data[column_index]]

        # Set the column units
        for column_name in units:
            table[column_name].unit = parse_unit(units[column_name])

        # Initialize
        initialize_table(table)

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
        initialize_table(table)

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

# FOR ISRF:

# def load_isrf(self):
#
#     """
#     This function ...
#     :return:
#     """
#
#     # M81_ds_isrf.dat
#
#     # Mean field intensities for all dust cells with nonzero absorption
#     # column 1: dust cell index
#     # column 2: x coordinate of cell center (pc)
#     # column 3: y coordinate of cell center (pc)
#     # column 4: z coordinate of cell center (pc)
#     # column 5: J_lambda (W/m3/sr) for lambda = 0.1 micron
#     # column 6: J_lambda (W/m3/sr) for lambda = 0.121153 micron
#     # column 7: J_lambda (W/m3/sr) for lambda = 0.14678 micron
#     # column 8: J_lambda (W/m3/sr) for lambda = 0.177828 micron
#     # column 9: J_lambda (W/m3/sr) for lambda = 0.215443 micron
#     # column 10: J_lambda (W/m3/sr) for lambda = 0.261016 micron
#     # column 11: J_lambda (W/m3/sr) for lambda = 0.316228 micron
#     # column 12: J_lambda (W/m3/sr) for lambda = 0.383119 micron
#     # column 13: J_lambda (W/m3/sr) for lambda = 0.464159 micron
#     # column 14: J_lambda (W/m3/sr) for lambda = 0.562341 micron
#     # column 15: J_lambda (W/m3/sr) for lambda = 0.681292 micron
#     # column 16: J_lambda (W/m3/sr) for lambda = 0.825404 micron
#     # column 17: J_lambda (W/m3/sr) for lambda = 1 micron
#     # column 18: J_lambda (W/m3/sr) for lambda = 1.21153 micron
#     # column 19: J_lambda (W/m3/sr) for lambda = 1.4678 micron
#     # column 20: J_lambda (W/m3/sr) for lambda = 1.77828 micron
#     # column 21: J_lambda (W/m3/sr) for lambda = 2.15443 micron
#     # column 22: J_lambda (W/m3/sr) for lambda = 2.61016 micron
#     # column 23: J_lambda (W/m3/sr) for lambda = 3.16228 micron
#     # column 24: J_lambda (W/m3/sr) for lambda = 3.83119 micron
#     # column 25: J_lambda (W/m3/sr) for lambda = 4.64159 micron
#     # column 26: J_lambda (W/m3/sr) for lambda = 5.62341 micron
#     # column 27: J_lambda (W/m3/sr) for lambda = 6.81292 micron
#     # column 28: J_lambda (W/m3/sr) for lambda = 8.25404 micron
#     # column 29: J_lambda (W/m3/sr) for lambda = 10 micron
#
#     # Determine the indices of the columns that have to be imported
#     column_indices = [0,1,2,3]
#     for i in range(4, 4 + self.nwavelengths): column_indices.append(i)
#
#     # Loop over the different contributions
#     for contribution in contributions:
#
#         # Determine the path to the output directory of the simulation
#         output_path = self.analysis_run.heating_output_path_for_contribution(contribution)
#
#         # Determine the path to the ISFR data file
#         isrf_path = fs.join(output_path, self.galaxy_name + "_ds_isrf.dat")
#
#         # Load the ISRF file
#         columns = np.loadtxt(isrf_path, usecols=column_indices, unpack=True)
#
#         # columns[0]: dust cell index
#         # columns[1]: x coordinate of cell center
#         # columns[2]: y coordinate of cell center
#         # columns[3]: z coordinate of cell center
#         # columns[4 -> 4 + (nwavelengths - 1)]: J_lambda
#
#         # Integrate over the J_lambda values to get the total bolometric absorbed luminosity per cell
#         luminosities = self.integrate_over_wavelengths(columns[4:4+self.nwavelengths])
#
#         #IDtot, x, y, z, Ltot = np.loadtxt(totISRFfile, usecols=(0, 1, 2, 3, 4,), unpack=True)
#         #IDold, Lold = np.loadtxt(oldISRFfile, usecols=(0, 4,), unpack=True)
#         #IDyng, Lyng = np.loadtxt(yngISRFfile, usecols=(0, 4,), unpack=True)
#         #IDnew, Lnew = np.loadtxt(newISRFfile, usecols=(0, 4,), unpack=True)
#
#         # write out
#         #writeLabsTot('Labs_tot.dat', IDtot, x, y, z, Ltot)
#         #writeLabs('Labs_old.dat', IDold, Lold)
#         #writeLabs('Labs_yng.dat', IDyng, Lyng)
#         #writeLabs('Labs_new.dat', IDnew, Lnew)
#
# # -----------------------------------------------------------------
#
# def integrate_over_wavelengths(self, jlambdas):
#
#     """
#     This function ...
#     :param jlambdas:
#     :return:
#     """
#
#     lum_cells = []
#
#     # L = sum_lambda (j_lambda * d_lambda)
#
#     # Get the wavelength deltas
#     deltas = self.wavelength_grid.deltas(asarray=True, unit="m") # deltas in meter
#
#     # Loop over all dust cells
#     for i in range(jlambdas[0].size):
#
#         # Put the Jlambda values for this cell in an array
#         jlambdas_cell = np.array([jlambdas[k][i] for k in range(len(jlambdas))])
#
#         # Calculate the luminosity
#         #lum = np.sum(jlambdas_cell * MjySr_to_LsunMicron * deltas)
#         lum = np.sum(jlambdas_cell * deltas) # if deltas are in meter, this is in W/m3/sr * m -> W/m2/sr
#         #lum = np.sum(jlambdas_cell * deltas * ...)
#
#         # Add the luminosity to the list
#         lum_cells.append(lum)
#
#     # Return the list of luminosities
#     return lum_cells

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

    # Debugging
    log.debug("Reading data ...")

    # Get the column data
    if method == "numpy": columns = np.loadtxt(path, unpack=True, ndmin=2)
    elif method == "pandas":
        import pandas as pd
        df = pd.read_csv(path, sep=" ", comment="#", header=None)
        ncolumns = len(df.columns)
        columns = [df[index].values for index in range(ncolumns)]
    else: raise ValueError("Invalid method: must be 'numpy' or 'pandas'")

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
