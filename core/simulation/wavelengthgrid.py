#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.wavelengthgrid Contains the WavelengthGrid class.

# -----------------------------------------------------------------

# Import standard modules
import numpy as np
import copy

# Import astronomical modules
from astropy.table import Table
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.tools import tables

# -----------------------------------------------------------------

class WavelengthGrid(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Attributes
        #self.table = Table(names=["Wavelength", "Delta"], dtype=('f8', 'f8'))
        #self.table = Table(names=["Wavelength"], dtype=['f8'])
        #self.table["Wavelength"].unit = Unit("micron")
        #self.table["Delta"].unit = Unit("micron")

        self.table = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Load the table
        grid.table = tables.from_file(path, format="ascii.ecsv")

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @property
    def has_deltas(self):

        """
        This function ...
        :return:
        """

        return "Delta" in self.table.colnames

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.table)

    # -----------------------------------------------------------------

    @classmethod
    def from_wavelengths(cls, wavelengths, unit="micron"):

        """
        This function ...
        :param wavelengths:
        :param unit:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Add the wavelengths
        grid.table = Table()
        grid.table["Wavelength"] = wavelengths
        grid.table["Wavelength"].unit = unit

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt_output(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Load the table
        table = tables.from_file(path, format="ascii")

        # Set the column names and units
        table.rename_column("col1", "Wavelength")
        table.rename_column("col2", "Delta")
        table["Wavelength"].unit = "micron"
        table["Delta"].unit = "micron"

        # Set the table
        grid.table = table

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt_input(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new class instance
        grid = cls()

        wavelengths = np.loadtxt(path, unpack=True, skiprows=1)

        # Create the table
        table = Table()
        table["Wavelength"] = wavelengths
        table["Wavelength"].unit = "micron"

        # Set the table
        grid.table = table

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.table)

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):

        """
        This function ...
        :return:
        """

        return np.min(self.table["Wavelength"])

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):

        """
        This function ...
        :return:
        """

        return np.max(self.table["Wavelength"])

    # -----------------------------------------------------------------

    def __getitem__(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.table["Wavelength"][index] * self.table["Wavelength"].unit

    # -----------------------------------------------------------------

    def add_point(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        self.table.add_row([wavelength])

    # -----------------------------------------------------------------

    def closest_wavelength(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return self.table["Wavelength"][self.closest_wavelength_index(wavelength)]

    # -----------------------------------------------------------------

    def closest_wavelength_index(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return tables.find_closest_index(self.table, wavelength, column_name="Wavelength")

    # -----------------------------------------------------------------

    def closest_wavelength_above_index(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return tables.find_closest_above_index(self.table, wavelength, column_name="Wavelength")

    # -----------------------------------------------------------------

    def closest_wavelength_below_index(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return tables.find_closest_below_index(self.table, wavelength, column_name="Wavelength")

    # -----------------------------------------------------------------

    def wavelength_indices(self, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        if min_wavelength is None and max_wavelength is None: return list(range(len(self.table)))
        elif min_wavelength is None: return list(range(self.closest_wavelength_below_index(max_wavelength)+1))
        elif max_wavelength is None: return list(range(self.closest_wavelength_above_index(min_wavelength), len(self.table)))
        else: return list(range(self.closest_wavelength_above_index(min_wavelength), self.closest_wavelength_below_index(max_wavelength)+1))

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

    def deltas(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Delta"], unit=unit)
        else: return tables.column_as_list(self.table["Delta"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def to_skirt_input(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        table_copy = copy.deepcopy(self.table)

        # Write the wavelength grid
        table_copy.rename_column("Wavelength", str(len(self.table)))  # Trick to have the number of wavelengths in the first line (required for SKIRT)
        tables.write(table_copy, path, format="ascii")

    # -----------------------------------------------------------------

    def to_skirt_output(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Write the table
        tables.write(self.table, path, format="ascii.ecsv")

# -----------------------------------------------------------------
