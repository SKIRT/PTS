#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.data Contains the Data3D class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.tools.utils import lazyproperty, ForbiddenOperation
from ...core.basics.table import SmartTable
from ...core.units.unit import PhotometricUnit, get_conversion_factor
from ...magic.core.frame import AllZeroError
from ...core.units.quantity import get_value_and_unit, add_with_units, subtract_with_units, multiply_with_units, divide_with_units
from ...core.basics.log import log
from ...core.tools import sequences, types
from ...core.tools.stringify import tostr, yes_or_no
from ...magic.core.frame import nan_value, inf_value, zero_value
from ...core.tools import formatting as fmt
from ...core.tools.numbers import weighed_arithmetic_mean_numpy, arithmetic_mean_numpy, weighed_median_numpy, median_numpy, weighed_standard_deviation_numpy, standard_deviation_numpy
from ...core.basics.range import RealRange, QuantityRange
from ...core.tools import filesystem as fs
from ...core.tools import strings
from ...core.units.parsing import parse_quantity as q
from ...core.units.parsing import parse_unit as u
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.data.sed import SED
from ...core.tools.parsing import real_list

# -----------------------------------------------------------------

x_coordinate_name = "x coordinates"
y_coordinate_name = "y coordinate"
z_coordinate_name = "z coordinate"
weight_name = "weight"

# -----------------------------------------------------------------

class Data3D(object):

    """
    This class represents 3D data THAT IS STATIC! (not like a table, which is modifiable)
    """

    def __init__(self, name, x, y, z, values, weights=None, length_unit=None, unit=None, description=None, distance=None,
                 wavelength=None, solid_angle=None):

        """
        The constructor ...
        :param x:
        :param y:
        :param z:
        :param values:
        :param weights:
        :param length_unit:
        :param unit:
        :param description:
        :param distance:
        :param wavelength:
        :param solid_angle:
        """

        # Set the name and description for the data
        self.name = name
        self.description = description

        # Set coordinates
        self.x = x
        self.y = y
        self.z = z

        # Set values
        self.values = values

        # Set weights
        self.weights = weights

        # Set units
        self.length_unit = u(length_unit) if length_unit is not None else None
        self.unit = u(unit) if unit is not None else None

        # Set conversion info
        self.distance = distance
        self.wavelength = wavelength
        self.solid_angle = solid_angle

        # Check sizes?
        self.check()

        # Path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, length_unit=None, unit=None, return_table=False):

        """
        This function ...
        :param path:
        :param length_unit:
        :param unit:
        :param return_table:
        :return:
        """

        # Debugging
        log.debug("Loading the 3D data in table format from '" + path + "' ...")

        # Read the table
        table = SmartTable.from_file(path, method="pandas")

        # Find the name of the variable
        standard_column_names = [x_coordinate_name.capitalize(), y_coordinate_name.capitalize(), z_coordinate_name.capitalize(), weight_name.capitalize()]
        column_name = sequences.get_single_other(table.column_names, standard_column_names, none="error", method="error")
        name = column_name.lower()

        # Get the units
        if length_unit is None: length_unit = table.get_column_unit(x_coordinate_name.capitalize())
        if unit is None: unit = table.get_column_unit(column_name)

        # Get the data
        x = table.get_column_array(x_coordinate_name.capitalize(), unit=length_unit, masked=False)
        y = table.get_column_array(y_coordinate_name.capitalize(), unit=length_unit, masked=False)
        z = table.get_column_array(z_coordinate_name.capitalize(), unit=length_unit, masked=False)
        values = table.get_column_array(column_name, unit=unit, masked=False)
        weights = table.get_column_array(weight_name.capitalize(), masked=False) if weight_name.capitalize() in table.colnames else None

        # Get the description
        description = table.meta["description"] if "description" in table.meta else None

        # Get the conversion info
        distance = table.meta["distance"] if "distance" in table.meta else None
        wavelength = table.meta["wavelength"] if "wavelength" in table.meta else None
        solid_angle = table.meta["solid_angle"] if "solid_angle" in table.meta else None

        # Create
        data = cls(name, x, y, z, values, weights=weights, length_unit=length_unit, unit=unit, description=description,
                   distance=distance, wavelength=wavelength, solid_angle=solid_angle)

        # Set the path
        data.path = path

        # Return the data
        if return_table: return data, table
        else: return data

    # -----------------------------------------------------------------

    @property
    def has_description(self):
        return self.description is not None

    # -----------------------------------------------------------------

    @property
    def has_length_unit(self):
        return self.length_unit is not None

    # -----------------------------------------------------------------

    @property
    def has_unit(self):
        return self.unit is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def radii(self):
        return np.sqrt(self.x ** 2 + self.y ** 2)

    # -----------------------------------------------------------------

    @property
    def nx(self):
        return len(self.x)

    # -----------------------------------------------------------------

    @property
    def ny(self):
        return len(self.y)

    # -----------------------------------------------------------------

    @property
    def nz(self):
        return len(self.z)

    # -----------------------------------------------------------------

    @property
    def nvalues(self):
        return len(self.values)

    # -----------------------------------------------------------------

    @property
    def has_weights(self):
        return self.weights is not None

    # -----------------------------------------------------------------

    @property
    def nweights(self):
        return len(self.weights) if self.has_weights else 0

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return self.nvalues

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        # Check types
        if not types.is_real_array(self.x): raise ValueError("x data must be an array of real values")
        if not types.is_real_array(self.y): raise ValueError("y data must be an array of real values")
        if not types.is_real_array(self.z): raise ValueError("z data must be an array of real values")
        if not types.is_real_array(self.values): raise ValueError("values must be an array of real values")
        if self.has_weights and not types.is_real_array(self.weights): raise ValueError("weights must be an array of real values")

        # Check
        sizes = [self.nx, self.ny, self.nz, self.nvalues]
        if not sequences.all_equal(sizes): raise ValueError("Sizes of data arrays are not equal")

        # Check weights
        if self.has_weights and self.nvalues != self.nweights: raise ValueError("Number of weights must be equal to number of values")

    # -----------------------------------------------------------------

    @lazyproperty
    def nans(self):
        return np.isnan(self.values)

    # -----------------------------------------------------------------

    @lazyproperty
    def infs(self):
        return np.isinf(self.values)

    # -----------------------------------------------------------------

    @lazyproperty
    def invalid(self):
        return self.nans + self.infs

    # -----------------------------------------------------------------

    @lazyproperty
    def ninvalid(self):
        return np.sum(self.invalid)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid(self):
        return np.logical_not(self.invalid)

    # -----------------------------------------------------------------

    @lazyproperty
    def nvalid(self):
        return np.sum(self.valid)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_x(self):
        return self.x[self.valid]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_y(self):
        return self.y[self.valid]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_z(self):
        return self.z[self.valid]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_values(self):
        return self.values[self.valid]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_weights(self):
        if not self.has_weights: return None
        return self.weights[self.valid]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_radii(self):
        return self.radii[self.valid]

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        return self.distance is not None

    # -----------------------------------------------------------------

    @property
    def has_wavelength(self):
        return self.wavelength is not None

    # -----------------------------------------------------------------

    @property
    def has_solid_angle(self):
        return self.solid_angle is not None

    # -----------------------------------------------------------------

    @property
    def conversion_info(self):

        """
        Thisfunction ...
        :return:
        """

        info = dict()
        if self.has_distance: info["distance"] = self.distance
        if self.has_wavelength: info["wavelength"] = self.wavelength
        if self.has_solid_angle: info["solid_angle"] = self.solid_angle
        return info

    # -----------------------------------------------------------------

    def sum(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        # Check
        unit = conversion_factor = None

        # Calculate
        result = np.nansum(self.values)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def mean(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        # Check
        unit = conversion_factor = None

        # Calculate
        result = weighed_arithmetic_mean_numpy(self.values, weights=self.weights) if self.has_weights else arithmetic_mean_numpy(self.values)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def median(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        # Check
        unit = conversion_factor = None

        # Calculate
        result = weighed_median_numpy(self.values, weights=self.weights) if self.has_weights else median_numpy(self.values)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def stddev(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        # Check
        unit = conversion_factor = None

        # Calculate
        result = weighed_standard_deviation_numpy(self.values, weights=self.weights) if self.has_weights else standard_deviation_numpy(self.values)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def __add__(self, other):

        """
        This fucntion ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values = add_with_units(self.values, self.unit, other_value, other_unit=other_unit, conversion_info=self.conversion_info)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=self.unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __iadd__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "addition")

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values = subtract_with_units(self.values, self.unit, other_value, other_unit=other_unit, conversion_info=self.conversion_info)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=self.unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __isub__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "subtraction")

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values, new_unit = multiply_with_units(self.values, self.unit, other_value, other_unit=other_unit)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=new_unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __imul__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "multiplication")

    # -----------------------------------------------------------------

    def __div__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values, new_unit = divide_with_units(self.values, self.unit, other_value, other_unit=other_unit)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=new_unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __idiv__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "division")

    # -----------------------------------------------------------------

    __truediv__ = __div__

    # -----------------------------------------------------------------

    __itruediv__ = __idiv__

    # -----------------------------------------------------------------

    def copy(self, copy_coordinates=False):

        """
        This function ...
        :param copy_coordinates:
        :return:
        """

        # Make copies of the data
        values = np.copy(self.values)
        weights = np.copy(self.weights) if self.has_weights else None
        x = np.copy(self.x) if copy_coordinates else self.x
        y = np.copy(self.y) if copy_coordinates else self.y
        z = np.copy(self.z) if copy_coordinates else self.z

        # Return the new copy
        return self.__class__(self.name, x, y, z, values, weights=weights, length_unit=self.length_unit, unit=self.unit,
                              description=self.description, distance=self.distance, wavelength=self.wavelength,
                              solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def normalized(self, to=1., return_sum=False, silent=False):

        """
        This function ...
        :param to:
        :param return_sum:
        :param silent:
        :return:
        """

        # Get the sum
        sum = self.sum()

        # Check whether the sum is nonnegative
        if sum < 0: raise RuntimeError("The sum of the data is negative")

        # Check if the sum is not zero
        if sum == 0: raise AllZeroError("The data cannot be normalized")

        # Calculate the conversion factor
        if hasattr(to, "unit"):  # quantity
            factor = to.value / sum
            unit = to.unit
        else:
            factor = to / sum
            unit = None

        # Debugging
        if not silent: log.debug("Multiplying the data with a factor of " + tostr(factor) + " to normalize to " + tostr(to) + " ...")

        # Create the normalized data and return
        new = self.converted_by_factor(factor, unit)
        if return_sum: return new, sum
        else: return new

    # -----------------------------------------------------------------

    def converted_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False,
                     brightness_strict=False, wavelength=None, solid_angle=None, silent=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param wavelength:
        :param solid_angle:
        :param silent:
        :return:
        """

        # Check the unit of the data is defined
        if not self.has_unit: raise ValueError("The unit of the data is not defined")

        # Parse "to unit": VERY IMPORTANT, BECAUSE DOING SELF.UNIT = TO_UNIT WILL OTHERWISE REPARSE AND WILL BE OFTEN INCORRECT!! (NO DENSITY OR BRIGHTNESS INFO)
        to_unit = u(to_unit, density=density, brightness=brightness, brightness_strict=brightness_strict, density_strict=density_strict)

        # Already in the correct unit
        if to_unit == self.unit:
            if not silent: log.debug("Data is already in the desired unit")
            return self.copy()

        # Debugging
        if not silent: log.debug("Converting the data from unit " + tostr(self.unit, add_physical_type=True) + " to unit " + tostr(to_unit, add_physical_type=True) + " ...")

        # Get the conversion factor
        factor = self._get_conversion_factor(to_unit, distance=distance, wavelength=wavelength, solid_angle=solid_angle, silent=silent)

        # Debugging
        if not silent: log.debug("Conversion factor: " + str(factor))

        # Return converted data
        return self.converted_by_factor(factor, to_unit)

    # -----------------------------------------------------------------

    def converted_by_factor(self, factor, new_unit, new_name=None, new_description=None):

        """
        This function ...
        :param factor:
        :param new_unit:
        :param new_name:
        :param new_description:
        :return:
        """

        # Create multiplicated data
        new = self * factor

        # Set the new unit
        new.unit = new_unit

        # Set new name?
        if new_name is not None: new.name = new_name
        if new_description is not None: new.description = new_description

        # Return the new data instance
        return new

    # -----------------------------------------------------------------

    @property
    def is_photometric(self):
        return self.has_unit and isinstance(self.unit, PhotometricUnit)

    # -----------------------------------------------------------------

    def _get_conversion_factor(self, to_unit, distance=None, wavelength=None, solid_angle=None, silent=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param wavelength:
        :param solid_angle:
        :param silent:
        :return:
        """

        # New: one central place to implement this
        return get_conversion_factor(self.unit, to_unit, parse=False, silent=silent, distance=distance,
                                     wavelength=wavelength, solid_angle=solid_angle, conversion_info=self.conversion_info)

        # # The data has a photometric unit
        # if self.is_photometric:
        #
        #     # Check that the target unit is also photometric
        #     if not isinstance(to_unit, PhotometricUnit): raise ValueError("Target unit is not photometric, while the data is")
        #
        #     # Set the conversion info
        #     #conversion_info = dict()
        #     conversion_info = self.conversion_info
        #     if distance is not None: conversion_info["distance"] = distance
        #     if wavelength is not None: conversion_info["wavelength"] = wavelength
        #     if solid_angle is not None: conversion_info["solid_angle"] = solid_angle
        #
        #     # Calculate the conversion factor
        #     factor = self.unit.conversion_factor(to_unit, silent=silent, **conversion_info)
        #
        # # The data does not have a photometric unit
        # else:
        #
        #     # Check whether target unit is also not photometric
        #     if isinstance(to_unit, PhotometricUnit): raise ValueError("Target unit is photometric, while the data is not")
        #
        #     # Calculate the conversion factor
        #     factor = self.unit.to(to_unit, silent=True)
        #
        # # Return
        # return factor

    # -----------------------------------------------------------------

    def where(self, value):
        return np.equal(self.values, value)

    # -----------------------------------------------------------------

    def where_not(self, value):
        return np.not_equal(self.values, value)

    # -----------------------------------------------------------------

    def where_smaller_than(self, value):
        return np.less(self.values, value)

    # -----------------------------------------------------------------

    def where_smaller_than_or_equal(self, value):
        return np.less_equal(self.values, value)

    # -----------------------------------------------------------------

    def where_greater_than(self, value):
        return np.greater(self.values, value)

    # -----------------------------------------------------------------

    def where_greater_than_or_equal(self, value):
        return np.greater_equal(self.values, value)

    # -----------------------------------------------------------------

    @lazyproperty
    def absolute_values(self):

        """
        This function ...
        :return:
        """

        return np.abs(self.values)

    # -----------------------------------------------------------------

    @property
    def nnans(self):
        return np.sum(self.nans)

    # -----------------------------------------------------------------

    @property
    def relative_nnans(self):
        return float(self.nnans) / self.nvalues

    # -----------------------------------------------------------------

    @property
    def has_nans(self):
        return np.any(self.nans)

    # -----------------------------------------------------------------

    @property
    def all_nans(self):
        return np.all(self.nans)

    # -----------------------------------------------------------------

    @property
    def ninfs(self):
        return np.sum(self.infs)

    # -----------------------------------------------------------------

    @property
    def relative_ninfs(self):
        return float(self.ninfs) / self.nvalues

    # -----------------------------------------------------------------

    @property
    def has_infs(self):
        return np.any(self.infs)

    # -----------------------------------------------------------------

    @property
    def all_infs(self):
        return np.all(self.infs)

    # -----------------------------------------------------------------

    @property
    def zeroes(self):
        return self.where(zero_value)

    # -----------------------------------------------------------------

    @property
    def nzeroes(self):
        return np.sum(self.zeroes)

    # -----------------------------------------------------------------

    @property
    def relative_nzeroes(self):
        return float(self.nzeroes) / self.nvalues

    # -----------------------------------------------------------------

    @property
    def has_zeroes(self):
        return np.any(self.zeroes)

    # -----------------------------------------------------------------

    @property
    def all_zeroes(self):
        return np.all(self.zeroes)

    # -----------------------------------------------------------------

    @property
    def min_x(self):
        return np.min(self.x)

    # -----------------------------------------------------------------

    @property
    def max_x(self):
        return np.max(self.x)

    # -----------------------------------------------------------------

    def x_range(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        if add_unit and self.has_length_unit: return QuantityRange(self.min_x, self.max_x, unit=self.length_unit)
        else: return RealRange(self.min_x, self.max_x)

    # -----------------------------------------------------------------

    @lazyproperty
    def x_span(self):
        return self.x_range().span

    # -----------------------------------------------------------------

    @lazyproperty
    def y_span(self):
        return self.y_range().span

    # -----------------------------------------------------------------

    @lazyproperty
    def z_span(self):
        return self.z_range().span

    # -----------------------------------------------------------------

    @property
    def min_y(self):
        return np.min(self.y)

    # -----------------------------------------------------------------

    @property
    def max_y(self):
        return np.max(self.y)

    # -----------------------------------------------------------------

    def y_range(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        if add_unit and self.has_length_unit: return QuantityRange(self.min_y, self.max_y, unit=self.length_unit)
        else: return RealRange(self.min_y, self.max_y)

    # -----------------------------------------------------------------

    @property
    def min_z(self):
        return np.min(self.z)

    # -----------------------------------------------------------------

    @property
    def max_z(self):
        return np.max(self.z)

    # -----------------------------------------------------------------

    def z_range(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        if add_unit and self.has_length_unit: return QuantityRange(self.min_z, self.max_z, unit=self.length_unit)
        else: return RealRange(self.min_z, self.max_z)

    # -----------------------------------------------------------------

    @lazyproperty
    def min_value(self):
        return np.nanmin(self.values)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_value(self):
        return np.nanmax(self.values)

    # -----------------------------------------------------------------

    def value_range(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        if add_unit and self.has_unit: return QuantityRange(self.min_value, self.max_value, unit=self.unit)
        else: return RealRange(self.min_value, self.max_value)

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

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Debugging
        log.debug("Saving the 3D data in table format to '" + path + "' ...")

        # Set the columns
        columns = [self.x, self.y, self.z, self.values]
        names = [x_coordinate_name, y_coordinate_name, z_coordinate_name, self.name]
        if self.has_weights:
            columns.append(self.weights)
            columns.append(weight_name)

        # Capitalize column names
        names = [name.capitalize() for name in names]

        # Set units
        units = [self.length_unit, self.length_unit, self.length_unit, self.unit]
        if self.has_weights: units.append(None)

        # Debugging
        log.debug("Creating table ...")

        # Create table
        table = SmartTable.from_columns(*columns, names=names, units=units, as_columns=True)

        # Debugging
        log.debug("Writing table ...")

        # Set the description
        if self.has_description: table.meta["description"] = self.description

        # Set the conversion info
        if self.has_distance: table.meta["distance"] = self.distance
        if self.has_wavelength: table.meta["wavelength"] = self.wavelength
        if self.has_solid_angle: table.meta["solid_angle"] = self.solid_angle

        # Save the table
        table.saveto(path)

        # Save the path
        if update_path: self.path = path

# -----------------------------------------------------------------

def show_data_properties(data):

    """
    This function ...
    :param data:
    :return:
    """

    print("")
    print(" - " + fmt.bold + "Name: " + fmt.reset_bold + data.name)
    if data.has_description: print(" - " + fmt.bold + "Description: " + fmt.reset_bold + data.description)
    if data.has_unit: print(" - " + fmt.bold + "Unit: " + fmt.reset_bold + tostr(data.unit))
    if data.has_length_unit: print(" - " + fmt.bold + "Length unit: " + fmt.reset_bold + tostr(data.length_unit))
    print(" - " + fmt.bold + "Number of points: " + fmt.reset_bold + str(data.nvalues))
    print(" - " + fmt.bold + "Weighed: " + fmt.reset_bold + yes_or_no(data.has_weights))
    print(" - " + fmt.bold + "Number of NaNs: " + fmt.reset_bold + str(data.nnans))
    print(" - " + fmt.bold + "Number of infs: " + fmt.reset_bold + str(data.ninfs))
    print(" - " + fmt.bold + "Number of zeroes: " + fmt.reset_bold + str(data.nzeroes))
    print(" - " + fmt.bold + "Number of invalid points: " + fmt.reset_bold + str(data.ninvalid))
    print(" - " + fmt.bold + "Number of valid points: " + fmt.reset_bold + str(data.nvalid))
    if data.has_distance: print(" - " + fmt.bold + "Distance: " + fmt.reset_bold + tostr(data.distance))
    if data.has_wavelength: print(" - " + fmt.bold + "Wavelength: " + fmt.reset_bold + tostr(data.wavelength))
    if data.has_solid_angle: print(" - " + fmt.bold + "Solid angle: " + fmt.reset_bold + tostr(data.solid_angle))
    print(" - " + fmt.bold + "Sum of all values: " + fmt.reset_bold + tostr(data.sum(add_unit=True)))
    print(" - " + fmt.bold + "Average value: " + fmt.reset_bold + tostr(data.mean(add_unit=True)))
    print(" - " + fmt.bold + "Median value: " + fmt.reset_bold + tostr(data.median(add_unit=True)))
    print(" - " + fmt.bold + "Standard deviation from the mean: " + fmt.reset_bold + tostr(data.stddev(add_unit=True)))
    print(" - " + fmt.bold + "Range of x coordinates: " + fmt.reset_bold + tostr(data.x_range(add_unit=True)))
    print(" - " + fmt.bold + "Range of y coordinates: " + fmt.reset_bold + tostr(data.y_range(add_unit=True)))
    print(" - " + fmt.bold + "Range of z coordinates: " + fmt.reset_bold + tostr(data.z_range(add_unit=True)))
    print(" - " + fmt.bold + "Range of values: " + fmt.reset_bold + tostr(data.value_range(add_unit=True)))
    print("")

# -----------------------------------------------------------------

class SpectralData3D(object):

    """
    This function ...
    """

    def __init__(self, name, x, y, z, wavelength_grid, filepath, length_unit=None, unit=None, description=None,
                 distance=None, solid_angle=None):

        """
        The constructor ...
        :param x:
        :param y:
        :param z:
        :param wavelength_grid:
        :param filepath: filepath of the data
        :param unit:
        :param description:
        :param distance:
        :param solid_angle:
        """

        # Set the name and description for the data
        self.name = name
        self.description = description

        # Set coordinates
        self.x = x
        self.y = y
        self.z = z

        # Set wavelength grid
        self.wavelength_grid = wavelength_grid

        # Set data filepath
        self.filepath = filepath

        # Set unit of the data and length unit
        self.length_unit = u(length_unit) if length_unit is not None else None
        self.unit = u(unit) if unit is not None else None

        # Set conversion info
        self.distance = distance
        self.solid_angle = solid_angle

        # Check sizes?
        self.check()

        # Initialize data structures
        self.spectral_data = [None for _ in range(self.ncoordinates)]
        self.spatial_data = [None for _ in range(self.nwavelengths)]

        # Path
        self.path = None

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        # Check types
        if not types.is_real_array(self.x): raise ValueError("x data must be an array of real values")
        if not types.is_real_array(self.y): raise ValueError("y data must be an array of real values")
        if not types.is_real_array(self.z): raise ValueError("z data must be an array of real values")

        # Check
        sizes = [self.nx, self.ny, self.nz, self.nrows]
        if not sequences.all_equal(sizes): raise ValueError("Sizes of data and file arrays are not equal")

    # -----------------------------------------------------------------

    @classmethod
    def from_table_file(cls, path, x, y, z, length_unit=None, name=None, description=None):

        """
        This function ...
        :param path:
        :param x:
        :param y:
        :param z:
        :param length_unit:
        :param name:
        :param description:
        :return:
        """

        # Load the header
        names = fs.get_column_names(path)

        # Get wavelengths
        wavelengths = []
        if strings.all_contains(names, "for lambda = "):

            # Loop over the lines
            for colname in names:
                wavelength_string = colname.split("for lambda = ")[1]
                wavelength = q(wavelength_string)
                wavelengths.append(wavelength)

        # Load the data
        #if name is None:
        common_name = strings.common_part(*names)
        if "for lambda =" in common_name: common_name = common_name.split("for lambda =")[0].strip()
        if common_name.endswith(")"):
            common_name, unit_string = strings.split_at_last(common_name[:-1], "(")
            common_name.strip()
            unit = u(unit_string)
        else: unit = None
        if name is None: name = common_name

        # Create wavelength grid
        wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths)

        # Return
        return cls(name, x, y, z, wavelength_grid, path, length_unit=length_unit, unit=unit, description=description)

    # -----------------------------------------------------------------

    @lazyproperty
    def header_lines(self):
        return fs.get_header_lines(self.filepath)

    # -----------------------------------------------------------------

    @property
    def nheader_lines(self):
        return len(self.header_lines)

    # -----------------------------------------------------------------

    @lazyproperty
    def nrows(self):
        return fs.get_nlines(self.filepath) - self.nheader_lines

    # -----------------------------------------------------------------

    @property
    def has_description(self):
        return self.description is not None

    # -----------------------------------------------------------------

    @property
    def has_length_unit(self):
        return self.length_unit is not None

    # -----------------------------------------------------------------

    @property
    def has_unit(self):
        return self.unit is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def radii(self):
        return np.sqrt(self.x ** 2 + self.y ** 2)

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        return self.distance is not None

    # -----------------------------------------------------------------

    @property
    def has_solid_angle(self):
        return self.solid_angle is not None

    # -----------------------------------------------------------------

    @property
    def conversion_info(self):

        """
        Thisfunction ...
        :return:
        """

        info = dict()
        if self.has_distance: info["distance"] = self.distance
        if self.has_solid_angle: info["solid_angle"] = self.solid_angle
        return info

    # -----------------------------------------------------------------

    @lazyproperty
    def nx(self):
        return len(self.x)

    # -----------------------------------------------------------------

    @lazyproperty
    def ny(self):
        return len(self.y)

    # -----------------------------------------------------------------

    @lazyproperty
    def nz(self):
        return len(self.z)

    # -----------------------------------------------------------------

    @property
    def ncoordinates(self):
        return self.nx

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths(self):
        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    def get_wavelength(self, index):
        return self.wavelength_grid[index]

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        return self.wavelength_grid.wavelengths(unit=unit, asarray=asarray, add_unit=add_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_array(self):
        return self.wavelengths(unit=self.wavelength_unit, asarray=True)

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):
        return self.wavelength_grid.unit

    # -----------------------------------------------------------------

    @property
    def wavelength_unit_string(self):
        return tostr(self.wavelength_unit)

    # -----------------------------------------------------------------

    def wavelengths_string(self, unit=None, add_unit=True, delimiter=", "):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param delimiter:
        :return:
        """

        # Set unit
        if unit is None: unit = self.wavelength_unit

        # Get scalar values
        wavelengths = self.wavelengths(unit=unit, add_unit=False)
        wavelength_strings = [tostr(wavelength, scientific=False) for wavelength in wavelengths]

        # Create string
        string = delimiter.join(wavelength_strings)

        # Add unit?
        if add_unit: string += " " + tostr(unit)

        # Return the string
        return string

    # -----------------------------------------------------------------

    def wavelength_deltas(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        return self.wavelength_grid.deltas(unit=unit, asarray=asarray, add_unit=add_unit)

    # -----------------------------------------------------------------

    def wavelength_indices(self, min_wavelength=None, max_wavelength=None, include_min=True, include_max=True):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param include_min:
        :param include_max:
        :return:
        """

        return self.wavelength_grid.wavelength_indices(min_wavelength, max_wavelength, include_min=include_min, include_max=include_max)

    # -----------------------------------------------------------------

    def get_indices(self, min_wavelength=None, max_wavelength=None, include_min=True, include_max=True):
        return self.wavelength_indices(min_wavelength=min_wavelength, max_wavelength=max_wavelength, include_min=include_min, include_max=include_max)

    # -----------------------------------------------------------------
    # SPECTRAL (FOR A COORDINATE)
    # -----------------------------------------------------------------

    def has_spectral_array(self, i):
        return self.spectral_data[i] is not None

    # -----------------------------------------------------------------

    def get_spectral_array(self, i):

        """
        This function ...
        :param i: 
        :return: 
        """

        # Read the data
        if self.has_spectral_array(i): self.spectral_data[i] = fs.get_row(self.filepath, i)

        # Return
        return self.spectral_data[i]

    # -----------------------------------------------------------------

    def get_sed(self, i):

        """
        This function ...
        :param i:
        :return:
        """

        # Get data
        array = self.get_spectral_array(i)

        # Create SED
        return SED.from_arrays(self.wavelength_array, array, self.wavelength_unit, self.unit)

    # -----------------------------------------------------------------
    # SPATIAL (FOR A WAVELENGTH)
    # -----------------------------------------------------------------

    def has_spatial_array(self, j):
        return self.spatial_data[j] is not None

    # -----------------------------------------------------------------

    def get_spatial_array(self, j):

        """
        Thisf ucntion ...
        :param j: 
        :return: 
        """

        return fs.get_column(self.filepath, j, float)

    # -----------------------------------------------------------------

    def get_data3d(self, j):

        """
        This function ...
        :param j:
        :return:
        """

        # Get data
        array = self.get_spatial_array(j)

        # Get the wavelength
        wavelength = self.get_wavelength(j)

        # Create Data3D
        return Data3D(self.name, self.x, self.y, self.z, array, length_unit=None, unit=self.unit, description=self.description,
                      wavelength=wavelength, distance=self.distance, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------`

    def get_index_for_wavelength(self, wavelength, return_wavelength=False):

        """
        This function ...
        :param wavelength:
        :param return_wavelength:
        :return:
        """

        return self.wavelength_grid.closest_wavelength_index(wavelength, return_wavelength=return_wavelength)

    # -----------------------------------------------------------------

    def get_data3d_for_wavelength(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        j = self.get_index_for_wavelength(wavelength)
        return self.get_data3d(j)

    # -----------------------------------------------------------------
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

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Debugging
        log.debug("Saving the spectral 3D data in table format to '" + path + "' ...")

        # Set the columns
        columns = [self.x, self.y, self.z]
        names = [x_coordinate_name, y_coordinate_name, z_coordinate_name]
        #if self.has_weights:
        #    columns.append(self.weights)
        #    columns.append(weight_name)

        # Capitalize column names
        names = [name.capitalize() for name in names]

        # Set units
        units = [self.length_unit, self.length_unit, self.length_unit]
        #if self.has_weights: units.append(None)

        # Debugging
        log.debug("Creating table ...")

        # Create table
        table = SmartTable.from_columns(*columns, names=names, units=units, as_columns=True)

        # Debugging
        log.debug("Writing table ...")

        # Set the name, description and unit
        table.meta["name"] = self.name
        if self.has_description: table.meta["description"] = self.description
        if self.has_unit: table.meta["unit"] = self.unit

        # Set the conversion info
        if self.has_distance: table.meta["distance"] = self.distance
        if self.has_solid_angle: table.meta["solid_angle"] = self.solid_angle

        # Set the data filepath
        table.meta["filepath"] = self.filepath

        # Set the wavelengths
        table.meta["wavelength_unit"] = self.wavelength_unit
        table.meta["wavelengths"] = self.wavelengths_string(add_unit=False)

        # Save the table
        table.saveto(path)

        # Update the path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Loading the spectral 3D data in table format from '" + path + "' ...")

        # Read the table
        table = SmartTable.from_file(path, method="pandas")

        # Get the units
        length_unit = table.get_column_unit(x_coordinate_name.capitalize())

        # Get the data
        x = table.get_column_array(x_coordinate_name.capitalize(), unit=length_unit)
        y = table.get_column_array(y_coordinate_name.capitalize(), unit=length_unit)
        z = table.get_column_array(z_coordinate_name.capitalize(), unit=length_unit)
        #values = table.get_column_array(column_name, unit=unit)
        #weights = table.get_column_array(weight_name.capitalize()) if weight_name.capitalize() in table.colnames else None

        # Get the name description, and unit
        name = table.meta["name"]
        description = table.meta["description"] if "description" in table.meta else None
        unit = table.meta["unit"] if "unit" in table.meta else None

        # Get the conversion info
        distance = table.meta["distance"] if "distance" in table.meta else None
        solid_angle = table.meta["solid_angle"] if "solid_angle" in table.meta else None

        # Get the filepath
        filepath = table.meta["filepath"]

        # Get the wavelengths
        wavelengths = real_list(table.meta["wavelengths"])
        wavelength_unit = table.meta["wavelength_unit"]
        wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths, unit=wavelength_unit)

        # Create
        # name, x, y, z, wavelength_grid, filepath, length_unit=None, unit=None, description=None, distance=None, solid_angle=None
        data = cls(name, x, y, z, wavelength_grid, filepath, length_unit=length_unit, unit=unit, description=description, distance=distance, solid_angle=solid_angle)

        # Set the path
        data.path = path

        # Return the data
        return data

# -----------------------------------------------------------------
