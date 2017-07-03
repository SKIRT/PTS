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
from ...core.tools import arrays

# -----------------------------------------------------------------

class WavelengthGrid(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        self.table = None

        # The path
        self.path = None

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

        # Set the path
        grid.path = path

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
    def from_wavelengths(cls, wavelengths, unit=None):

        """
        This function ...
        :param wavelengths:
        :param unit:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Check if wavelengths are values or quantities
        for index in range(len(wavelengths)):

            if unit is None:
                if hasattr(wavelengths[index], "unit"):
                    unit = wavelengths[index].unit
                    wavelengths[index] = wavelengths[index].value
                else: pass
            else:
                if hasattr(wavelengths[index], "unit"):
                    wavelengths[index] = wavelengths[index].to(unit).value
                else: pass

        # Add the wavelengths
        grid.table = Table()
        grid.table["Wavelength"] = wavelengths
        grid.table["Wavelength"].unit = unit

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @classmethod
    def logarithmic(cls, wrange, npoints):

        """
        This function ...
        :param wrange:
        :param npoints:
        :return:
        """

        # Verify the grid parameters
        if npoints < 2: raise ValueError("the number of points in the grid should be at least 2")
        if wrange.min <= 0: raise ValueError("the shortest wavelength should be positive")

        # Calculate log of boundaries
        logmin = np.log10(float(wrange.min))
        logmax = np.log10(float(wrange.max))

        # Calculate the grid points
        base_grid = np.logspace(logmin, logmax, num=npoints, endpoint=True, base=10)

        # Create the wavelength grid
        wavelength_grid = cls.from_wavelengths(base_grid)

        # Return the wavelength grid
        return wavelength_grid

    # -----------------------------------------------------------------

    @classmethod
    def nested_log(cls, wrange, npoints, wrange_zoom, npoints_zoom):

        """
        This function ...
        :param wrange:
        :param npoints:
        :param wrange_zoom:
        :param npoints_zoom:
        :return:
        """

        # Verify the grid parameters
        if npoints < 2: raise ValueError("the number of points in the low-resolution grid should be at least 2")
        if npoints_zoom < 2: raise ValueError("the number of points in the high-resolution subgrid should be at least 2")
        if wrange.min <= 0: raise ValueError("the shortest wavelength should be positive")
        if (wrange_zoom.min <= wrange.min
            or wrange_zoom.max <= wrange_zoom.min
            or wrange.max <= wrange_zoom.max):
            raise ValueError("the high-resolution subgrid should be properly nested in the low-resolution grid")

        logmin = np.log10(float(wrange.min))
        logmax = np.log10(float(wrange.max))
        logmin_zoom = np.log10(float(wrange_zoom.min))
        logmax_zoom = np.log10(float(wrange_zoom.max))

        # Build the high- and low-resolution grids independently
        base_grid = np.logspace(logmin, logmax, num=npoints, endpoint=True, base=10)
        zoom_grid = np.logspace(logmin_zoom, logmax_zoom, num=npoints_zoom, endpoint=True, base=10)

        # Merge the two grids
        total_grid = []

        # Add the wavelengths of the low-resolution grid before the first wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength < wrange_zoom.min: total_grid.append(wavelength)

        # Add the wavelengths of the high-resolution grid
        for wavelength in zoom_grid: total_grid.append(wavelength)

        # Add the wavelengths of the low-resolution grid after the last wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength > wrange_zoom.max: total_grid.append(wavelength)

        # Create the wavelength grid
        wavelength_grid = cls.from_wavelengths(total_grid)

        # Return the wavelength grid
        return wavelength_grid

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

        return arrays.find_closest_index(self.table["Wavelength"], wavelength, array_unit=self.table["Wavelength"].unit)

    # -----------------------------------------------------------------

    def closest_wavelength_above_index(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return arrays.find_closest_above_index(self.table["Wavelength"], wavelength, array_unit=self.table["Wavelength"].unit)

    # -----------------------------------------------------------------

    def closest_wavelength_below_index(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return arrays.find_closest_below_index(self.table["Wavelength"], wavelength, array_unit=self.table["Wavelength"].unit)

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

        if asarray: return arrays.plain_array(self.table["Wavelength"], unit=unit, array_unit=self.table["Wavelength"].unit)
        else: return arrays.array_as_list(self.table["Wavelength"], unit=unit, add_unit=add_unit, array_unit=self.table["Wavelength"].unit)

    # -----------------------------------------------------------------

    def deltas(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return arrays.plain_array(self.table["Delta"], unit=unit, array_unit=self.table["Delta"].unit)
        else: return arrays.array_as_list(self.table["Delta"], unit=unit, add_unit=add_unit, array_unit=self.table["Delta"].unit)

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

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether the path is valid
        if self.path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Write the table
        tables.write(self.table, path, format="ascii.ecsv")

        # Update the path
        self.path = path

# -----------------------------------------------------------------
