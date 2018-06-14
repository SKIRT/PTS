#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.curve Contains the Curve class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import copy
from scipy.interpolate import interp1d

# Import astronomical modules
from astropy.units import spectral

# Import the relevant PTS classes and modules
from .table import SmartTable
from ..filter.filter import parse_filter
from ..units.parsing import parse_unit as u
from ..tools import arrays
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter

# -----------------------------------------------------------------

class Curve(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        if kwargs.get("from_astropy", None) is None:
            if "x_unit" in kwargs: from_astropy = False
            else: from_astropy = True
        else: from_astropy = kwargs.pop("from_astropy")

        if not from_astropy:

            # Get properties
            x_unit = kwargs.pop("x_unit", None)
            y_unit = kwargs.pop("y_unit", None)
            x_name = kwargs.pop("x_name", "x")
            y_name = kwargs.pop("y_name", "y")
            x_description = kwargs.pop("x_description", "x values")
            y_description = kwargs.pop("y_description", "y values")

        # Call the constructor of the base class
        super(Curve, self).__init__(*args, **kwargs)

        if not from_astropy:

            # Set the column info
            self.add_column_info(x_name, float, x_unit, x_description)
            self.add_column_info(y_name, float, y_unit, y_description)

            # Set x name and y name
            self.x_name = x_name
            self.y_name = y_name

        #else: print(kwargs)

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        new = copy.deepcopy(self)

        if hasattr(self, "x_name"): new.x_name = self.x_name
        if hasattr(self, "y_name"): new.y_name = self.y_name

        return new

    # -----------------------------------------------------------------

    @classmethod
    def from_columns(cls, *columns, **kwargs):

        """
        This function ...
        :param columns:
        :param kwargs:
        :return:
        """

        kwargs["from_astropy"] = False

        # x and y name
        if kwargs.get("names", None) is not None:
            names = kwargs.get("names")
            x_name = names[0]
            y_name = names[1]
            kwargs["x_name"] = x_name
            kwargs["y_name"] = y_name

        # x and y unit
        if kwargs.get("units", None) is not None:
            units = kwargs.get("units")
            x_unit = units[0]
            y_unit = units[1]
            kwargs["x_unit"] = x_unit
            kwargs["y_unit"] = y_unit

        # x and y descriptions
        #if kwargs.get("descriptions", None) is not None:
        descriptions = kwargs.get("descriptions", ("no description", "no description"))
        x_description = descriptions[0]
        y_description = descriptions[1]
        kwargs["x_description"] = x_description
        kwargs["y_description"] = y_description

        # USe the base class implementation
        curve = super(Curve, cls).from_columns(*columns, **kwargs)

        # Set x name and y name
        curve.x_name = curve.column_names[0]
        curve.y_name = curve.column_names[1]

        # Return the curve
        return curve

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        # Load the curve
        curve = super(Curve, cls).from_file(path)

        # Set x name and y name
        curve.x_name = curve.column_info[0][0]
        curve.y_name = curve.column_info[1][0]

        # Return the curve
        return curve

    # -----------------------------------------------------------------

    @classmethod
    def from_data_file(cls, path, x_name="x", x_description="x values", x_unit=None, y_name="y",
                       y_description="y values", y_unit=None, x_column=0, y_column=1, skiprows=0, conversion_info=None):

        """
        This function ...
        :param path:
        :param x_name:
        :param x_description:
        :param x_unit:
        :param y_name:
        :param y_description:
        :param y_unit:
        :param x_column:
        :param y_column:
        :param skiprows:
        :param conversion_info:
        :return:
        """

        # Define columns
        columns = (x_column, y_column)

        kwargs = dict()

        kwargs["x_name"] = x_name
        kwargs["x_description"] = x_description
        kwargs["x_unit"] = x_unit

        kwargs["y_name"] = y_name
        kwargs["y_description"] = y_description
        kwargs["y_unit"] = y_unit

        # Create the curve
        curve = cls(**kwargs)
        x_values, y_values = np.loadtxt(path, unpack=True, usecols=columns, skiprows=skiprows)
        for x_value, y_value in zip(x_values, y_values): curve.add_point(x_value, y_value, conversion_info=conversion_info)

        # Return the curve
        return curve

    # -----------------------------------------------------------------

    @property
    def x_unit(self):

        """
        This function ...
        :return:
        """

        return self[self.x_name].unit

    # -----------------------------------------------------------------

    @property
    def y_unit(self):

        """
        This function ...
        :return:
        """

        return self[self.y_name].unit

    # -----------------------------------------------------------------

    @property
    def x_data(self):

        """
        This function ...
        :return:
        """

        return self[self.x_name]

    # -----------------------------------------------------------------

    @property
    def y_data(self):

        """
        This function ...
        :return:
        """

        return self[self.y_name]

    # -----------------------------------------------------------------

    def add_point(self, x_value, y_value, conversion_info=None, sort=True):

        """
        This function ...
        :param x_value:
        :param y_value:
        :param conversion_info:
        :param sort:
        :return:
        """

        # Set values
        values = [x_value, y_value]

        # Add a row
        self.add_row(values, conversion_info=conversion_info)

        # Sort the table by the x values
        if sort: self.sort(self.x_name)

    # -----------------------------------------------------------------

    @property
    def has_errors(self):

        """
        This function ...
        :return:
        """

        return "Error-" in self.colnames and "Error+" in self.colnames

    # -----------------------------------------------------------------

    def normalize(self, value=1.0, method="integral"):

        """
        This function ...
        :param value:
        :param method:
        :return:
        """

        if method == "max":

            max_value = np.max(self[self.y_name])
            factor = value / max_value
            self[self.y_name] *= factor
            self[self.y_name].unit = None

        elif method == "sum":

            sum_value = np.sum(self[self.y_name])
            factor = value / sum_value
            self[self.y_name] *= factor
            self[self.y_name].unit = None

        elif method == "integral": raise NotImplementedError("Not implemented yet")
        else: raise ValueError("Invalid option for 'method'")

    # -----------------------------------------------------------------

    def get_x(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :return:
        """

        # Create and return
        if asarray: return arrays.plain_array(self[self.x_name], unit=unit, array_unit=self.column_unit(self.x_name),
                                      conversion_info=conversion_info, density=density, brightness=brightness)
        else: return arrays.array_as_list(self[self.x_name], unit=unit, add_unit=add_unit,
                                        array_unit=self.column_unit(self.x_name), conversion_info=conversion_info,
                                        density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def get_y(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :return:
        """

        # Create and return
        if asarray: return arrays.plain_array(self[self.y_name], unit=unit, array_unit=self.column_unit(self.y_name),
                                      conversion_info=conversion_info, density=density, brightness=brightness)
        else: return arrays.array_as_list(self[self.y_name], unit=unit, add_unit=add_unit,
                                        array_unit=self.column_unit(self.y_name), conversion_info=conversion_info,
                                        density=density, brightness=brightness)

    # -----------------------------------------------------------------

    @property
    def npoints(self):

        """
        Thisf unction
        :return:
        """

        return len(self)

    # -----------------------------------------------------------------

    def get_indices(self, x_min=None, x_max=None, include_min=True, include_max=True):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param include_min:
        :param include_max:
        :return:
        """

        # No limits given
        if x_min is None and x_max is None: return list(range(self.npoints))

        # Get the values
        x_values = self.get_x()

        indices = []

        # Loop over the values
        for index, value in enumerate(x_values):

            # Checks
            if x_min is not None:
                if include_min:
                    if value < x_min: continue
                else:
                    if value <= x_min: continue
            if x_max is not None:
                if include_max:
                    if value > x_max: continue
                else:
                    if value >= x_max: continue

            # Add the index
            indices.append(index)

        # Return the indices
        return indices

    # -----------------------------------------------------------------

    def splice(self, x_min=None, x_max=None, include_min=True, include_max=True):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param include_min:
        :param include_max:
        :return:
        """

        # Get the indices
        indices = self.get_indices(x_min, x_max, include_min=include_min, include_max=include_max) # don't name arguments because of re-definition of function in WavelengthCurve class

        # Set the x and y unit
        x_unit = self.x_unit
        y_unit = self.y_unit

        # Get the values
        x_values = [self.get_value(self.x_name, index, unit=x_unit, add_unit=False) for index in indices]
        y_values = [self.get_value(self.y_name, index, unit=y_unit, add_unit=False) for index in indices]

        # Set the names and units
        names = (self.x_name, self.y_name)
        units = (x_unit, y_unit)

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def flatten_above(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        from ..units.parsing import parse_quantity

        # Check value with unit
        if self.x_unit is not None:
            if hasattr(value, "unit"): pass
            elif isinstance(value, basestring): value = parse_quantity(value)
            else: raise ValueError("Unit of the value is not defined")
        elif hasattr(value, "unit") or isinstance(value, basestring): raise ValueError("Unit of the value is defined, but column unit is not")

        # Get the indices of the values to be flattened
        #indices = self.get_indices(x_min=value, include_min=include)
        indices = self.get_indices(value, None, include_min=include) # works also with derived class implementation

        # Set flatten value with unit
        if self.y_unit is not None:
            if hasattr(flatten_value, "unit"): pass
            elif flatten_value == 0.: flatten_value = flatten_value * self.y_unit
            else: raise ValueError("Unit of the flatten value is not defined")
        elif hasattr(flatten_value, "unit") or isinstance(flatten_value, basestring): raise ValueError("Unit of the flatten value is defined, but column unit is not")
        #print(flatten_value)

        # Flatten values
        for index in indices: self.set_value(self.y_name, index, flatten_value)

    # -----------------------------------------------------------------

    def flatten_below(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        from ..units.parsing import parse_quantity

        # Check value with unit
        if self.x_unit is not None:
            if hasattr(value, "unit"): pass
            elif isinstance(value, basestring): value = parse_quantity(value)
            else: raise ValueError("Unit of the value is not defined")
        elif hasattr(value, "unit") or isinstance(value, basestring): raise ValueError("Unit of the value is defined, but column unit is not")

        # Get the indices of the values to be flattened
        #indices = self.get_indices(x_max=value, include_max=include)
        indices = self.get_indices(None, value, include_max=include) # works also with derived class implementation

        # Set flatten value with unit
        if self.y_unit is not None:
            if hasattr(flatten_value, "unit"): pass
            elif isinstance(flatten_value, basestring): flatten_value = parse_quantity(flatten_value)
            elif flatten_value == 0.: flatten_value = flatten_value * self.y_unit
            else: raise ValueError("Unit of the flatten value is not defined")
        elif hasattr(flatten_value, "unit") or isinstance(flatten_value, basestring): raise ValueError("Unit of the flatten value is defined, but column unit is not")
        #print(flatten_value)

        # Flatten values
        for index in indices: self.set_value(self.y_name, index, flatten_value)

    # -----------------------------------------------------------------

    def flattened_above(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        # Make copy
        new = self.copy()
        new.flatten_above(value, flatten_value=flatten_value, include=include)
        return new

    # -----------------------------------------------------------------

    def flattened_below(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :return:
        """

        # Make copy
        new = self.copy()
        new.flatten_below(value, flatten_value=flatten_value, include=include)
        return new

    # -----------------------------------------------------------------

    def __iadd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # HACKS
        #if hasattr(self, "x_name"): print("self", self.x_name)
        #if hasattr(other, "x_name"): print("other", other.x_name)
        #if hasattr(self, "y_name"): print("self", self.y_name)
        #if hasattr(other, "y_name"): print("other", other.y_name)

        #x_name = self.x_name if hasattr(self, "x_name") else other.x_name
        #y_name = self.y_name if hasattr(self, "y_name") else other.y_name
        #x_unit = self.x_unit if hasattr(self, "x_unit") and self.x_unit is not None else other.x_unit
        #y_unit = self.y_unit if hasattr(self, "y_unit") and self.y_unit is not None else other.y_unit
        #print(x_unit, y_unit)

        # Check whether names are the same
        if self.x_name != other.x_name: raise ValueError("x name must be the same")
        if self.y_name != other.y_name: raise ValueError("y name must be the same")
        y_name = self.y_name
        y_unit = self.y_unit

        # Check whether the lengths are the same
        if self.npoints != other.npoints: raise ValueError("number of data points must be the same")

        # Get the conversion factor
        conversion_factor = other.y_unit.conversion_factor(y_unit)

        # Get the values of other in the same unit
        #other_values = np.array([other.get_value(self.y_name, index).to(self.y_unit).value for index in range(self.npoints)])
        other_values = np.asarray(other[y_name]) * conversion_factor

        # Add the data
        self[y_name] += other_values

        # HACKS
        #if not hasattr(self, "x_name"): self.x_name = x_name
        #if not hasattr(self, "y_name"): self.y_name = y_name
        #if not hasattr(self, "x_unit") or self.x_unit is None: self[x_name].unit = x_unit
        #if not hasattr(self, "y_unit") or self.y_unit is None: self[y_name].unit = y_unit

        # Return
        return self

    # -----------------------------------------------------------------

    def __add__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Check whether x name is the same
        if self.x_name != other.x_name: raise ValueError("x name must be the same (" + str(self.x_name) + " and " + str(other.x_name) + ")")
        x_name = self.x_name

        # Get the units
        x_unit = self.x_unit
        y_unit = self.y_unit
        if x_unit is not None and other.x_unit is None: raise ValueError("unit of '" + other.x_name + "' is not defined")
        if y_unit is not None and other.y_unit is None: raise ValueError("unit of '" + other.y_name + "' is not defined")

        # Initialize a list for the x values and y values
        x_values = []
        y_values = []

        # Loop over the values of this curve and the other curve simultaneously
        i = 0
        j = 0
        while True:

            # Get the values
            x_a = self.get_value(self.x_name, i, unit=x_unit, add_unit=False)
            x_b = other.get_value(other.x_name, j, unit=x_unit, add_unit=False)

            # Value is the same: add
            if x_a == x_b:

                result = self.get_value(self.y_name, i, unit=y_unit, add_unit=False) + other.get_value(other.y_name, j, unit=y_unit, add_unit=False)

                x_values.append(x_a)
                y_values.append(result)

                # Increment
                i += 1
                j += 1

            # x of a is greater than x of b
            elif x_a > x_b: j += 1

            # x of b is greater than x of a
            else: i += 1

            # Check for termination
            if i >= self.npoints: break
            if j >= other.npoints: break

        # Set the new y name
        if self.y_name == other.y_name: y_name = self.y_name
        else: y_name = self.y_name + " + " + other.y_name

        # Set the names and units
        names = (x_name, y_name)
        units = (x_unit, y_unit)
        #print("names", names)
        #print("units", units)

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Check whether x name is the same
        if self.x_name != other.x_name: raise ValueError("x name must be the same (" + str(self.x_name) + " and " + str(other.x_name) + ")")
        x_name = self.x_name

        # Get the units
        x_unit = self.x_unit
        y_unit = self.y_unit
        if x_unit is not None and other.x_unit is None: raise ValueError("unit of '" + other.x_name + "' is not defined")
        if y_unit is not None and other.y_unit is None: raise ValueError("unit of '" + other.y_name + "' is not defined")

        # Initialize a list for the x values and y values
        x_values = []
        y_values = []

        # Loop over the values of this curve and the other curve simultaneously
        i = 0
        j = 0
        while True:

            # Get the values
            x_a = self.get_value(self.x_name, i, unit=x_unit, add_unit=False)
            x_b = other.get_value(other.x_name, j, unit=x_unit, add_unit=False)

            # Value is the same: add
            if x_a == x_b:

                result = self.get_value(self.y_name, i, unit=y_unit, add_unit=False) - other.get_value(other.y_name, j, unit=y_unit, add_unit=False)

                x_values.append(x_a)
                y_values.append(result)

                # Increment
                i += 1
                j += 1

            # x of a is greater than x of b
            elif x_a > x_b: j += 1

            # x of b is greater than x of a
            else: i += 1

            # Check for termination
            if i >= self.npoints: break
            if j >= other.npoints: break

        # Set the new y name
        if self.y_name == other.y_name: y_name = self.y_name
        else: y_name = self.y_name + " - " + other.y_name

        # Set the names and units
        names = (x_name, y_name)
        units = (x_unit, y_unit)
        #print("names", names)
        #print("units", units)

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Check whether x name is the same
        if self.x_name != other.x_name: raise ValueError("x name must be the same (" + str(self.x_name) + " and " + str(other.x_name) + ")")
        x_name = self.x_name

        # Initialize a list for the x values and y values
        x_values = []
        y_values = []

        # Loop over the values of this curve and the other curve simultaneously
        i = 0
        j = 0
        while True:

            # Get the values
            x_a = self.get_value(self.x_name, i)
            x_b = other.get_value(other.x_name, j)

            # Value is the same: add
            if x_a == x_b:

                result = self.get_value(self.y_name, i) * other.get_value(other.y_name, j)

                x_values.append(x_a)
                y_values.append(result)

                # Increment
                i += 1
                j += 1

            # x of a is greater than x of b
            elif x_a > x_b: j += 1

            # x of b is greater than x of a
            else: i += 1

            # Check for termination
            if i >= self.npoints: break
            if j >= other.npoints: break

        # Set the new y name
        if self.y_name == other.y_name: y_name = self.y_name + "1" + " * " + other.y_name + "2"
        else: y_name = self.y_name + " * " + other.y_name

        # Create new curve
        return Curve.from_columns(x_values, y_values, names=[x_name, y_name])

    # -----------------------------------------------------------------

    def __div__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Check whether x name is the same
        if self.x_name != other.x_name: raise ValueError("x name must be the same (" + str(self.x_name) + " and " + str(other.x_name) + ")")
        x_name = self.x_name

        # Initialize a list for the x values and y values
        x_values = []
        y_values = []

        # Loop over the values of this curve and the other curve simultaneously
        i = 0
        j = 0
        while True:

            # Get the values
            x_a = self.get_value(self.x_name, i)
            x_b = other.get_value(other.x_name, j)

            # Value is the same: add
            if x_a == x_b:

                result = self.get_value(self.y_name, i) / other.get_value(other.y_name, j)

                x_values.append(x_a)
                y_values.append(result)

                # Increment
                i += 1
                j += 1

            # x of a is greater than x of b
            elif x_a > x_b: j += 1

            # x of b is greater than x of a
            else: i += 1

            # Check for termination
            if i >= self.npoints: break
            if j >= other.npoints: break

        # Set the new y name
        if self.y_name == other.y_name: y_name = self.y_name + "1" + " / " + other.y_name + "2"
        else: y_name = self.y_name + " / " + other.y_name

        # Create new curve
        return Curve.from_columns(x_values, y_values, names=[x_name, y_name])

    # -----------------------------------------------------------------

    __truediv__ = __div__

# -----------------------------------------------------------------

class WavelengthCurve(Curve):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        if "y_name" in kwargs: from_astropy = False
        else: from_astropy = True

        if not from_astropy:

            # Get x properties
            wavelength_unit = kwargs.pop("wavelength_unit", "micron")
            if wavelength_unit is None: wavelength_unit = "micron"

            # Get y properties
            name = kwargs.pop("y_name")
            description = kwargs.pop("y_description")
            unit = kwargs.pop("y_unit", None)

            x_unit = wavelength_unit
            y_unit = unit
            x_name = "Wavelength"
            y_name = name
            x_description = "Wavelength"
            y_description = description

            kwargs["x_unit"] = x_unit
            kwargs["y_unit"] = y_unit
            kwargs["x_name"] = x_name
            kwargs["y_name"] = y_name
            kwargs["x_description"] = x_description
            kwargs["y_description"] = y_description

            # Set distance
            self.distance = kwargs.pop("distance", None)

        # From astropy call ...
        else: self.distance = None

        # Call the constructor of the base class
        super(WavelengthCurve, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def get_indices(self, min_wavelength=None, max_wavelength=None, include_min=True, include_max=True):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param include_min:
        :param include_max:
        :return:
        """

        # Call the implementation of the base class
        return super(WavelengthCurve, self).get_indices(x_min=min_wavelength, x_max=max_wavelength, include_min=include_min, include_max=include_max)

    # -----------------------------------------------------------------

    def splice(self, min_wavelength=None, max_wavelength=None, include_min=True, include_max=True):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param include_min:
        :param include_max:
        :return:
        """

        # Call the implementation of the base class
        return super(WavelengthCurve, self).splice(x_min=min_wavelength, x_max=max_wavelength, include_min=include_min, include_max=include_max)

    # -----------------------------------------------------------------

    def add_point(self, wavelength, value, conversion_info=None, sort=True):

        """
        This function ...
        :param wavelength:
        :param value:
        :param conversion_info:
        :param sort:
        :return:
        """

        # Set conversion info
        if conversion_info is None:
            conversion_info_value = dict()
            conversion_info_value["wavelength"] = wavelength
            if self.distance is not None: conversion_info_value["distance"] = self.distance
            conversion_info = {self.value_name: conversion_info_value}

        # Add the point, passing the conversion info
        super(WavelengthCurve, self).add_point(wavelength, value, conversion_info=conversion_info, sort=sort)

    # -----------------------------------------------------------------

    @property
    def value_name(self):

        """
        This function ..
        :return:
        """

        # Check if setup has been performed
        if len(self.colnames) == 0: self._setup()

        for index in reversed(range(len(self.colnames))):
            name = self.colnames[index]
            if name == "Error+" or name == "Error-": continue
            return self.colnames[index]

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        return self.y_unit

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):

        """
        This function ...
        :return:
        """

        return self.x_unit

    # -----------------------------------------------------------------

    def wavelength_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Wavelength", index)

    # -----------------------------------------------------------------

    def get_wavelength(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Wavelength", index)

    # -----------------------------------------------------------------

    def value_for_index(self, index, add_unit=True, unit=None, density=False, brightness=False, conversion_info=None):

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
        value = self.get_value(self.value_name, index, add_unit=True)

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

    def value_for_wavelength(self, wavelength, unit=None, add_unit=True, density=False, brightness=False,
                             interpolate=True, conversion_info=None):

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

        if interpolate:
            interpolated = interp1d(self.wavelengths(unit="micron", asarray=True), self.values(self.unit, asarray=True), kind='linear')
            value = interpolated(wavelength.to("micron").value) * self.unit
        else:
            from ..tools import sequences
            index = sequences.find_closest_index(self.wavelengths(unit="micron", add_unit=False), wavelength.to("micron").value)
            value = self[self.value_name][index] * self.unit

        # Convert unit if necessary
        if unit is not None:

            # Create conversion info
            if conversion_info is None: conversion_info = dict()
            conversion_info["wavelength"] = wavelength
            if self.distance is not None: conversion_info["distance"] = self.distance

            unit = u(unit, density=density, brightness=brightness)
            value = value.to(unit, **conversion_info)

        # Remove unit if requested
        if not add_unit: value = value.value

        return value

    # -----------------------------------------------------------------

    def wavelengths_mask(self, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :return: 
        """

        # Initialize mask
        mask = np.zeros(len(self), dtype=bool)

        # Loop over the wavelengths, check them
        for index, wavelength in enumerate(self.wavelengths()):

            if min_wavelength is not None and wavelength < min_wavelength: mask[index] = True
            if max_wavelength is not None and wavelength > max_wavelength: mask[index] = True

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):

        """
        This ufnction ...
        :return:
        """

        return self.get_value("Wavelength", 0)

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.get_value("Wavelength", -1)

    # -----------------------------------------------------------------

    def covers(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return self.min_wavelength < wavelength < self.max_wavelength

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        # Create mask
        if min_wavelength is not None or max_wavelength is not None: mask = self.wavelengths_mask(min_wavelength=min_wavelength, max_wavelength=max_wavelength)
        else: mask = None

        # Create and return
        if asarray: return arrays.plain_array(self["Wavelength"], unit=unit, array_unit=self.column_unit("Wavelength"), mask=mask)
        else: return arrays.array_as_list(self["Wavelength"], unit=unit, add_unit=add_unit, array_unit=self.column_unit("Wavelength"), mask=mask)

    # -----------------------------------------------------------------

    def wavelength_grid(self, unit=None, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param unit:
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid
        return WavelengthGrid.from_wavelengths(self.wavelengths(unit=unit, min_wavelength=min_wavelength, max_wavelength=max_wavelength))

    # -----------------------------------------------------------------

    def wavelength_deltas(self, unit=None, asarray=False, add_unit=True, min_wavelength=None, max_wavelength=None):

        """
        Thisn function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        return self.wavelength_grid(unit=unit, min_wavelength=min_wavelength, max_wavelength=max_wavelength).deltas(unit=unit, asarray=asarray, add_unit=add_unit)

    # -----------------------------------------------------------------

    def frequencies(self, unit="Hz", asarray=False, add_unit=True, min_wavelength=None, max_wavelength=None):

        """
        This function ...
        :param unit:
        :param asarray: 
        :param add_unit: 
        :param min_wavelength:
        :param max_wavelength:
        :return: 
        """

        # Create mask
        if min_wavelength is not None or max_wavelength is not None: mask = self.wavelengths_mask(min_wavelength=min_wavelength, max_wavelength=max_wavelength)
        else: mask = None

        # Create and return
        if asarray: return arrays.plain_array(self["Wavelength"], unit=unit, array_unit=self.column_unit("Wavelength"), equivalencies=spectral(), mask=mask)
        else: return arrays.array_as_list(self["Wavelength"], unit=unit, add_unit=add_unit, array_unit=self.column_unit("Wavelength"), equivalencies=spectral(), mask=mask)

    # -----------------------------------------------------------------

    def values(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False,
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

        # Create mask
        if min_wavelength is not None or max_wavelength is not None: mask = self.wavelengths_mask(min_wavelength=min_wavelength, max_wavelength=max_wavelength)
        else: mask = None

        # Create conversion info
        if conversion_info is None: conversion_info = dict()
        conversion_info["wavelengths"] = self.wavelengths()
        if self.distance is not None: conversion_info["distance"] = self.distance

        # Create and return
        if asarray: return arrays.plain_array(self[self.value_name], unit=unit, array_unit=self.column_unit(self.value_name),
                                              conversion_info=conversion_info, density=density, brightness=brightness,
                                              mask=mask)
        else: return arrays.array_as_list(self[self.value_name], unit=unit, add_unit=add_unit,
                                          array_unit=self.column_unit(self.value_name), conversion_info=conversion_info,
                                          density=density, brightness=brightness, mask=mask)

# -----------------------------------------------------------------

class FilterCurve(WavelengthCurve):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FilterCurve, self).__init__(*args, **kwargs)

        # Set column info
        self.column_info.insert(0, ("Band", str, None, "band"))
        self.column_info.insert(0, ("Instrument", str, None, "instrument"))
        self.column_info.insert(0, ("Observatory", str, None, "observatory"))

    # -----------------------------------------------------------------

    def add_point(self, fltr, value, conversion_info=None):

        """
        This function ...
        :param fltr:
        :param value:
        :param conversion_info:
        :return:
        """

        # Construct the
        values = [fltr.observatory, fltr.instrument, fltr.band, fltr.wavelength, value]

        # Create conversion info
        if conversion_info is None:
            conversion_info_value = dict()
            conversion_info_value["wavelength"] = fltr.wavelength
            conversion_info = {self.value_name: conversion_info_value}

        # Add the row
        self.add_row(values, conversion_info=conversion_info)

        # Sort the table by the x values
        self.sort(self.x_name)

    # -----------------------------------------------------------------

    def value_for_band(self, instrument, band, unit=None, add_unit=True, density=False, brightness=False):

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

        has_unit = self[self.value_name].unit is not None
        has_mask = hasattr(self[self.value_name], "mask")

        # If unit has to be converted, check whether the original unit is specified
        if not has_unit and unit is not None: raise ValueError("Cannot determine the unit of quantity so values cannot be converted to " + str(unit))

        # Parse the target unit
        if unit is not None: unit = u(unit, density=density, brightness=brightness)

        # Loop over all the entries in the table
        for i in range(len(self)):

            # Get instrument and band
            instrument_entry = self.get_value("Instrument", i)
            band_entry = self.get_value("Band", i)

            if not (instrument_entry == instrument and band_entry == band): continue

            # if the entry is masked, return None
            if has_mask and self.is_masked_value(self.value_name, i): return None

            # The column has a unit, we can convert if necessary
            if has_unit:

                # Add the unit initially to be able to convert
                #value = self[self.value_name][i] * self[self.value_name].unit
                value = self.get_value(self.value_name, i)

                # If a target unit is specified, convert
                if unit is not None: value = value.to(unit).value * u(unit)

                # Strip unit if requested
                if not add_unit: value = value.value

            # No unit for the column
            else: value = self.get_value(self.value_name, i)

            # Return the value / quantity ...
            return value

        # If no match is found, return None
        return None

    # -----------------------------------------------------------------

    def value_for_filter(self, fltr, unit=None, add_unit=True, density=False, brightness=False):

        """
        This function ...
        :param fltr:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :return:
        """

        return self.value_for_band(fltr.instrument, fltr.band, unit, add_unit, density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def instruments(self):

        """
        This function ...
        :return:
        """

        return arrays.array_as_list(self["Instrument"])

    # -----------------------------------------------------------------

    def bands(self):

        """
        This function ...
        :return:
        """

        return arrays.array_as_list(self["Band"])

    # -----------------------------------------------------------------

    def filter_names(self):

        """
        This function ...
        :return:
        """

        # Initialize
        names = []

        # Loop over all entries
        for i in range(len(self)):

            # Get the instrument and band
            instrument = str(self.get_value("Instrument", i))
            band = str(self.get_value("Band", i))

            # Add the filter name
            names.append(instrument + " " + band)

        # Return the list of filter names
        return names

    # -----------------------------------------------------------------

    def filters(self):

        """
        This function ...
        :return:
        """

        # Initialize
        filters = []

        # Loop over all entries
        for i in range(len(self)):

            # Get the filter
            fltr = self.get_filter(i)

            # Add the filter to the list
            filters.append(fltr)

        # Return the list of filters
        return filters

    # -----------------------------------------------------------------

    def has_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Loop over the rows, check instrument and band
        for index in range(len(self)):

            # Match found
            if self.get_value("Instrument", index) == fltr.instrument and self.get_value("Band", index) == fltr.band: return True

        # No match found
        return False

    # -----------------------------------------------------------------

    def index_for_filter(self, fltr, return_none=False):

        """
        This function ...
        :param fltr:
        :param return_none:
        :return:
        """

        # Loop over the rows, check instrument and band
        for index in range(len(self)):

            # Match?
            if self.get_value("Instrument", index) == fltr.instrument and self.get_value("Band", index) == fltr.band: return index

        # No match
        if return_none: return None
        else: raise ValueError("Filter not in the curve")

    # -----------------------------------------------------------------

    def only_broad_band(self):

        """
        This function ...
        :return:
        """

        # Make a copy of this curve
        new = self.copy()

        # Loop over the rows, remove the row if it does not correspond to a broad band filter
        is_broad = self.broad_band_filters()
        for index in reversed(range(len(self))):
            if not is_broad[index]: new.remove_row(index)

        # Return the new SED
        return new

    # -----------------------------------------------------------------

    def only_narrow_band(self):

        """
        This function ...
        :return:
        """

        # Make a copy of this curve
        new = self.copy()

        # Loop over the rows, remove the row if it does not correspond to a narrow band filter
        is_narrow = self.narrow_band_filters()
        for index in reversed(range(len(self))):
            if not is_narrow[index]: new.remove_row(index)

        # Return the new SED
        return new

    # -----------------------------------------------------------------

    def broad_band_filters(self):

        """
        This function ...
        :return:
        """

        return [isinstance(fltr, BroadBandFilter) for fltr in self.filters()]

    # -----------------------------------------------------------------

    def narrow_band_filters(self):

        """
        This function ...
        :return:
        """

        return [isinstance(fltr, NarrowBandFilter) for fltr in self.filters()]

    # -----------------------------------------------------------------

    def get_filter(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get instrument and band
        instrument = self.get_value("Instrument", index)
        band = self.get_value("Band", index)

        # Parse the filter
        return parse_filter(instrument + " " + band)

# -----------------------------------------------------------------
