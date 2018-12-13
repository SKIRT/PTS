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
import math
import copy

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools import tables
from ...core.tools import arrays
from ...core.tools.stringify import stringify_list_fancy
from ...core.tools import strings

# -----------------------------------------------------------------

def is_from_skirt_input(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..tools import filesystem as fs

    # Get the first two lines
    first, second = fs.get_first_lines(path, 2)

    # Check whether the first line is an integer
    if not strings.is_integer(first): return False

    # Check whether the value of the first line (= the number of wavelengths) is greater than the second value (= the actual first wavelength)
    # OK yeah this is a bit of a hack ...
    if float(second) > float(first): return False

    # Checks passed for SKIRT
    return True

# -----------------------------------------------------------------

def is_ecsv(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..tools import filesystem as fs
    return "ECSV" in fs.get_first_line(path)

# -----------------------------------------------------------------

def load_wavelength_grid(path, **kwargs):

    """
    This function ...
    :param path:
    :param kwargs:
    :return:
    """

    # From SKIRT?
    if is_from_skirt_input(path): return WavelengthGrid.from_skirt_input(path, **kwargs)

    # TODO: From SKIRT output?

    # ECSV?
    elif is_ecsv(path): return WavelengthGrid.from_file(path, **kwargs)

    # Other?
    else: return WavelengthGrid.from_text_file(path, kwargs.pop("unit"), **kwargs)

# -----------------------------------------------------------------

class WavelengthGrid(object):

    """
    This class ...
    """

    def __init__(self, table=None):

        """
        The constructor ...
        """

        # The wavelength grid table
        self.table = table

        # The path
        self.path = None

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return WavelengthGrid(table=self.table.copy())

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, format="ascii.ecsv"):

        """
        This function ...
        :param path:
        :param format:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Load the table
        grid.table = tables.from_file(path, format=format)

        # Set the path
        grid.path = path

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @classmethod
    def from_text_file(cls, path, unit, column=0, skiprows=0):

        """
        This function ...
        :param path:
        :param unit:
        :param column:
        :param skiprows:
        :return:
        """

        # Determine the columns to use
        columns = (column)

        # Load the data
        wavelength_column = np.loadtxt(path, dtype=float, unpack=True, skiprows=skiprows, usecols=columns)

        # Create the SED
        return cls.from_wavelengths(wavelength_column, unit=unit)

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
    def from_sed(cls, sed):

        """
        Thisf unction ...
        :param sed:
        :return:
        """

        # Create from the wavelengths of the SED
        return cls.from_wavelengths(sed.wavelengths())

    # -----------------------------------------------------------------

    @classmethod
    def from_wavelengths(cls, wavelengths, unit=None, sort=False):

        """
        This function ...
        :param wavelengths:
        :param unit:
        :param sort:
        :return:
        """

        # Make a copy of the list of wavelengths
        wavelengths = wavelengths[:]

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

        # Sort?
        if sort: grid.sort()

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
    def from_skirt_output(cls, path, ignore_deltas=False):

        """
        This function ...
        :param path:
        :param ignore_deltas:
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
        #table["Delta"].unit = "micron"

        # Get delta values
        deltas_skirt = np.array(table["Delta"])

        # Remove delta column
        table.remove_column("Delta")

        # Set the table
        grid.table = table

        # Compare deltas
        deltas = grid.deltas(unit="micron", asarray=True)

        # Get absolute and relative differences
        #absdiff = abs(deltas_skirt - deltas)
        #reldiff = absdiff / deltas
        #print("max abs:", np.max(absdiff))
        #print("max rel:", np.max(reldiff))
        if not ignore_deltas and not np.all(np.isclose(deltas_skirt, deltas, atol=1e-5, rtol=1e-6)): raise IOError("The delta wavelengths are not consistent with the wavelength points. Call with ignore_deltas=True to ignore this and to use re-computed delta values.")

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt_input(cls, path, remote=None):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        # Create a new class instance
        grid = cls()

        if remote is not None:
            lines = remote.get_lines(path, add_sep=True)
            wavelengths = np.loadtxt(lines, unpack=True, skiprows=1)
        else: wavelengths = np.loadtxt(path, unpack=True, skiprows=1)

        # Create the table
        table = Table()
        table["Wavelength"] = wavelengths
        table["Wavelength"].unit = "micron"

        # Set the table
        grid.table = table

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        #return str(self.table)
        return "[" + stringify_list_fancy(self.wavelengths(unit="micron", add_unit=False))[1] + "] micron"

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.table)

    # -----------------------------------------------------------------

    def covers(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        scalar = wavelength.to("micron").value
        return self.min_wavelength < scalar < self.max_wavelength

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):
        return np.min(self.table["Wavelength"])

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):
        return np.max(self.table["Wavelength"])

    # -----------------------------------------------------------------

    @property
    def range(self):
        from ..basics.range import QuantityRange
        return QuantityRange(self.min_wavelength, self.max_wavelength, unit=self.unit)

    # -----------------------------------------------------------------

    @property
    def mean_wavelength(self):
        from ..tools import numbers
        return numbers.arithmetic_mean(*self.wavelengths(unit=self.unit, add_unit=False)) * self.unit

    # -----------------------------------------------------------------

    @property
    def geometric_mean_wavelength(self):
        from ..tools import numbers
        return numbers.geometric_mean(*self.wavelengths(unit=self.unit, add_unit=False)) * self.unit

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_indices):

        """
        This function ...
        :param index_or_indices:
        :return:
        """

        from ..tools import types

        # Single index
        if types.is_integer_type(index_or_indices): return self.table["Wavelength"][index_or_indices] * self.table["Wavelength"].unit

        # Multiple indices
        elif types.is_integer_sequence(index_or_indices) or types.is_integer_array(index_or_indices) or types.is_integer_tuple(index_or_indices):

            # Get the wavelengths and create a new grid
            wavelengths = self.wavelengths(asarray=True, unit=self.unit)[list(index_or_indices)]
            return WavelengthGrid.from_wavelengths(wavelengths, unit=self.unit)

        # Invalid
        else: raise ValueError("Invalid index or indices: '" + str(index_or_indices) + "'")

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

        return self.table["Wavelength"][self.closest_wavelength_index(wavelength)] * self.table["Wavelength"].unit

    # -----------------------------------------------------------------

    def closest_wavelength_index(self, wavelength, return_wavelength=False):

        """
        This function ...
        :param wavelength:
        :param return_wavelength:
        :return:
        """

        index = arrays.find_closest_index(self.table["Wavelength"], wavelength, array_unit=self.table["Wavelength"].unit)
        if return_wavelength: return index, self.table["Wavelength"][index] * self.table["Wavelength"].unit
        else: return index

    # -----------------------------------------------------------------

    def closest_wavelength_above_index(self, wavelength, inclusive=False):

        """
        This function ...
        :param wavelength:
        :param inclusive:
        :return:
        """

        return arrays.find_closest_above_index(self.table["Wavelength"], wavelength, array_unit=self.table["Wavelength"].unit, inclusive=inclusive)

    # -----------------------------------------------------------------

    def closest_wavelength_below_index(self, wavelength, inclusive=False):

        """
        This function ...
        :param wavelength:
        :param inclusive:
        :return:
        """

        return arrays.find_closest_below_index(self.table["Wavelength"], wavelength, array_unit=self.table["Wavelength"].unit, inclusive=inclusive)

    # -----------------------------------------------------------------

    def wavelength_indices(self, min_wavelength=None, max_wavelength=None, wavelength_range=None, include_min=False, include_max=False):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param wavelength_range:
        :param include_min:
        :param include_max:
        :return:
        """

        # Set min and max from range
        if wavelength_range is not None:
            if min_wavelength is not None: raise ValueError("Minimum wavelength cannot be specified when wavelength range is specified")
            if max_wavelength is not None: raise ValueError("Maximum wavelength cannot be specified when wavelength range is specified")
            min_wavelength = wavelength_range.min
            max_wavelength = wavelength_range.max

        # Get indices
        if min_wavelength is None and max_wavelength is None: return list(range(len(self.table)))
        elif min_wavelength is None:
            last_index = self.closest_wavelength_below_index(max_wavelength, inclusive=include_max)
            return list(range(last_index+1))
        elif max_wavelength is None:
            first_index = self.closest_wavelength_above_index(min_wavelength, inclusive=include_min)
            return list(range(first_index, len(self.table)))
        else:
            first_index = self.closest_wavelength_above_index(min_wavelength, inclusive=include_min)
            last_index = self.closest_wavelength_below_index(max_wavelength, inclusive=include_max)
            return list(range(first_index, last_index+1))

    # -----------------------------------------------------------------

    def wavelengths_mask(self, min_wavelength=None, max_wavelength=None, inclusive=True, wavelength_range=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param inclusive:
        :param wavelength_range:
        :return:
        """

        # Set min and max from range
        if wavelength_range is not None:
            if min_wavelength is not None: raise ValueError("Minimum wavelength cannot be specified when wavelength range is specified")
            if max_wavelength is not None: raise ValueError("Maximum wavelength cannot be specified when wavelength range is specified")
            min_wavelength = wavelength_range.min
            max_wavelength = wavelength_range.max

        # Initialize mask
        mask = np.zeros(len(self), dtype=bool)

        # Loop over the wavelengths, check them
        for index, wavelength in enumerate(self.wavelengths()):

            if min_wavelength is not None:
                if inclusive and wavelength < min_wavelength: mask[index] = True
                if not inclusive and wavelength <= min_wavelength: mask[index] = True
            if max_wavelength is not None:
                if inclusive and wavelength > max_wavelength: mask[index] = True
                if not inclusive and wavelength >= max_wavelength: mask[index] = True

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def get_subgrid(self, min_wavelength=None, max_wavelength=None, wavelength_range=None, return_indices=False,
                    include_min=False, include_max=False):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param wavelength_range:
        :param return_indices:
        :param include_min:
        :param include_max:
        :return:
        """

        # Get the indices
        indices = self.wavelength_indices(min_wavelength=min_wavelength, max_wavelength=max_wavelength,
                                          wavelength_range=wavelength_range, include_min=include_min, include_max=include_max)

        # Return the subgrid
        if return_indices: return self[indices], indices
        else: return self[indices]

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True, min_wavelength=None, max_wavelength=None, inclusive=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param min_wavelength:
        :param max_wavelength:
        :param inclusive:
        :return:
        """

        # Create mask
        if min_wavelength is not None or max_wavelength is not None: mask = self.wavelengths_mask(min_wavelength=min_wavelength, max_wavelength=max_wavelength, inclusive=inclusive)
        else: mask = None

        if asarray: return arrays.plain_array(self.table["Wavelength"], unit=unit, array_unit=self.table["Wavelength"].unit, mask=mask)
        else: return arrays.array_as_list(self.table["Wavelength"], unit=unit, add_unit=add_unit, array_unit=self.table["Wavelength"].unit, mask=mask)

    # -----------------------------------------------------------------

    @property
    def unit(self):
        return self.table["Wavelength"].unit

    # -----------------------------------------------------------------

    def minima(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        # Initialize a list to contain the values
        values = []

        # Check
        if asarray and unit is None: raise ValueError("You cannot know which units the values are going to be if you don't specify the target unit")
        if not add_unit and unit is None: raise ValueError("You cannot know which units the values are going to be if you don't specifiy the target unit and you put add_unit to False")

        # Loop over the wavelengths
        for index in range(self.nwavelengths):

            # Get the value for the minimum
            if index == 0: value = self.table["Wavelength"][0] * self.unit
            else: value = math.sqrt(self.table["Wavelength"][index - 1] * self.table["Wavelength"][index]) * self.unit

            # Convert to specified unit
            if unit is not None: value = value.to(unit)

            # Add unit?
            if not add_unit or asarray: value = value.value

            # Add
            values.append(value)

        # Return
        if asarray: return np.array(values)
        else: return values

    # -----------------------------------------------------------------

    def maxima(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        # Initialize a list to contain the values
        values = []

        # Check
        if asarray and unit is None: raise ValueError("You cannot know which units the values are going to be if you don't specify the target unit")
        if not add_unit and unit is None: raise ValueError("You cannot know which units the values are going to be if you don't specifiy the target unit and you put add_unit to False")

        # Loop over the wavelengths
        for index in range(self.nwavelengths):

            # Get the value for the maximum
            if index == self.nwavelengths - 1: value = self.table["Wavelength"][-1] * self.unit
            else: value = math.sqrt(self.table["Wavelength"][index] * self.table["Wavelength"][index + 1]) * self.unit

            # Convert to specified unit
            if unit is not None: value = value.to(unit)

            # Add unit?
            if not add_unit or asarray: value = value.value

            # Add
            values.append(value)

        # Return
        if asarray: return np.array(values)
        else: return values

    # -----------------------------------------------------------------

    def deltas(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        values = []
        minima = self.minima()
        maxima = self.maxima()

        # Check
        if asarray and unit is None: raise ValueError("You cannot know which units the values are going to be if you don't specify the target unit")
        if not add_unit and unit is None: raise ValueError("You cannot know which units the values are going to be if you don't specifiy the target unit and you put add_unit to False")

        # Loop over the wavelengths
        for index in range(self.nwavelengths):

            # Calcualte delta
            value = maxima[index] - minima[index]

            # Convert to specified unit
            if unit is not None: value = value.to(unit)

            # Add unit?
            if not add_unit or asarray: value = value.value

            # Add
            values.append(value)

        # Return
        if asarray: return np.array(values)
        else: return values

    # -----------------------------------------------------------------

    def logminima(self, unit, asarray=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :return:
        """

        minima = self.minima(unit=unit, asarray=True)
        logminima = np.log10(minima)
        if asarray: return logminima
        else: return list(logminima)

    # -----------------------------------------------------------------

    def logmaxima(self, unit, asarray=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :return:
        """

        maxima = self.maxima(unit=unit, asarray=True)
        logminima = np.log10(maxima)
        if asarray: return logminima
        else: return list(logminima)

    # -----------------------------------------------------------------

    def logdeltas(self, asarray=False):

        """
        Thisf unction ...
        :param asarray:
        :return:
        """

        values = []
        # Unit do not matter as long as they're the same
        logminima = self.logminima(unit=self.unit)
        logmaxima = self.logmaxima(unit=self.unit)

        # Loop over the wavelengths
        for index in range(self.nwavelengths):

            # Calculate log delta
            value = logmaxima[index] - logminima[index]

            # Add value
            values.append(value)

        # Return
        if asarray: return np.array(values)
        else: return values

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

    def sort(self):

        """
        This function ...
        :return:
        """

        self.table.sort("Wavelength")

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

    def plot(self, min_y=0, max_y=1, y_value=0.5, unit=None, size=(10,2), color="b", marker="o", markersize=2):

        """
        This function ...
        :param min_y:
        :param max_y:
        :param y_value:
        :param unit:
        :param size:
        :param color:
        :param marker:
        :param markersize:
        :return:
        """

        import matplotlib.pyplot as plt

        plt.figure(figsize=size)

        if unit is None: unit = self.unit
        wavelengths = self.wavelengths(unit=unit, asarray=True)
        y_values = [y_value for _ in wavelengths]
        plt.scatter(wavelengths, y_values, c=color, marker=marker, s=markersize)

        # Labels
        plt.xlabel("Wavelength (" + str(unit) + ")")

        # Set limits and scaling
        min_wavelength = 0.9 * np.min(wavelengths)
        max_wavelength = 1.1 * np.max(wavelengths)
        axes = plt.gca()
        axes.set_xlim(min_wavelength, max_wavelength)
        axes.set_ylim(min_y, max_y)
        axes.set_xscale("log")

        from ..basics.plot import ScalarFormatterForceFormat

        # Set scalar formatting
        xaxis = axes.get_xaxis()
        xaxis.set_major_formatter(ScalarFormatterForceFormat(useOffset=True, useMathText=True))

        # Hide y axis
        yaxis = axes.get_yaxis()
        yaxis.set_visible(False)
        yaxis.set_ticks([])

        # Fix
        plt.tight_layout()

        # Show
        plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_deltas(self, unit=None, color="b"):

        """
        This function ...
        :param unit:
        :param color:
        :return:
        """

        import matplotlib.pyplot as plt

        plt.figure()

        # Get wavelengths and deltas
        if unit is None: unit = self.unit
        wavelengths = self.wavelengths(unit=unit, asarray=True)
        deltas = self.deltas(unit=unit, asarray=True)

        # Plot
        plt.plot(wavelengths, deltas, c=color)

        # Labels
        plt.ylabel("Delta wavelength (" + str(unit) + ")")
        plt.xlabel("Wavelength (" + str(unit) + ")")

        # Set limits and scaling
        min_wavelength = 0.9 * np.min(wavelengths)
        max_wavelength = 1.1 * np.max(wavelengths)
        axes = plt.gca()
        axes.set_xlim(min_wavelength, max_wavelength)
        #axes.set_ylim(min_y, max_y)
        axes.set_xscale("log")
        axes.set_yscale("log")

        from ..basics.plot import ScalarFormatterForceFormat

        # Set scalar formatting
        xaxis = axes.get_xaxis()
        xaxis.set_major_formatter(ScalarFormatterForceFormat(useOffset=True, useMathText=True))
        yaxis = axes.get_yaxis()
        yaxis.set_major_formatter(ScalarFormatterForceFormat(useOffset=True, useMathText=True))

        # Fix
        plt.tight_layout()

        # Show
        plt.show()
        plt.close()

    # -----------------------------------------------------------------

    def plot_logdeltas(self, unit=None, color="b"):

        """
        This function ...
        :param unit:
        :param color:
        :return:
        """

        import matplotlib.pyplot as plt

        plt.figure()

        # Get wavelengths and log-deltas
        if unit is None: unit = self.unit
        wavelengths = self.wavelengths(unit=unit, asarray=True)
        logdeltas = self.logdeltas(asarray=True)

        # Plot
        plt.plot(wavelengths, logdeltas, c=color)

        # Labels
        plt.ylabel("Log-delta wavelength (dex)")
        plt.xlabel("Wavelength (" + str(unit) + ")")

        # Set limits and scaling
        min_wavelength = 0.9 * np.min(wavelengths)
        max_wavelength = 1.1 * np.max(wavelengths)
        axes = plt.gca()
        axes.set_xlim(min_wavelength, max_wavelength)
        # axes.set_ylim(min_y, max_y)
        axes.set_xscale("log")
        #axes.set_yscale("log")

        from ..basics.plot import ScalarFormatterForceFormat

        # Set scalar formatting
        xaxis = axes.get_xaxis()
        xaxis.set_major_formatter(ScalarFormatterForceFormat(useOffset=True, useMathText=True))

        # Fix
        plt.tight_layout()

        # Show
        plt.show()
        plt.close()

# -----------------------------------------------------------------
