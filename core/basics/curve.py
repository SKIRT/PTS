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
import warnings
from scipy.interpolate import interp1d

# Import astronomical modules
from astropy.units import spectral

# Import the relevant PTS classes and modules
from ..filter.filter import parse_filter
from ..units.parsing import parse_unit as u
from ..units.unit import get_common_unit
from ..tools import arrays
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter
from .relation import Relation
from .range import RealRange

# -----------------------------------------------------------------

class Curve(Relation):

    """
    This class ...
    """

    def add_point(self, x_value, y_value, conversion_info=None, sort=True):

        """
        This function ...
        :param x_value:
        :param y_value:
        :param conversion_info:
        :param sort: DEFAULT IS TRUE HERE
        :return:
        """

        # Call the implementation of the base class
        super(Curve, self).add_point(x_value, y_value, conversion_info=conversion_info, sort=sort)

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

    def __iadd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

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

        # Return
        return self

    # -----------------------------------------------------------------

    def __add__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        from ..tools import types
        from ..units.unit import get_converted_value

        # Real or integer value
        if types.is_real_or_integer(other):

            # Get real value
            other = float(other)

            # Set columns
            x_values = self.x_array
            y_values = self.y_array + other

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, self.y_unit)

        # Add quantity
        elif types.is_quantity(other):

            # Get value in correct unit
            other = get_converted_value(other, self.y_unit)

            # Set columns
            x_values = self.x_array
            y_values = self.y_array + other

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, self.y_unit)

        # Add other curve
        elif isinstance(other, Curve):

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
            #print("ADDING CURVES")
            while True:

                # Get the values
                x_a = self.get_value(self.x_name, i, unit=x_unit, add_unit=False)
                x_b = other.get_value(other.x_name, j, unit=x_unit, add_unit=False)

                #print(x_a, x_b, x_a==x_b)

                # Value is the same: add
                #if x_a == x_b:
                isclose = np.isclose(x_a, x_b)
                #if not isclose: print("NOT CLOSE:", x_a, x_b)
                if isclose:

                    # Try to get wavelength and distance for unit conversion
                    conversion_info = {}

                    # Set wavelength
                    if isinstance(self, WavelengthCurve): conversion_info["wavelength"] = self.get_wavelength(i)
                    elif isinstance(other, WavelengthCurve): conversion_info["wavelength"] = other.get_wavelength(j)

                    # Set distance
                    if hasattr(self, "distance") and self.distance is not None: conversion_info["distance"] = self.distance
                    elif hasattr(other, "distance") and other.distance is not None: conversion_info["distance"] = other.distance
                    if len(conversion_info) == 0: conversion_info = None

                    # Calculate the sum of the y values
                    result = self.get_value(self.y_name, i, unit=y_unit, add_unit=False, conversion_info=conversion_info) + other.get_value(other.y_name, j, unit=y_unit, add_unit=False, conversion_info=conversion_info)

                    # Add the x value and the new y value
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

        # Invalid type
        else: raise TypeError("Cannot add " + self.__class__.__name__ + " and a " + str(type(other).__name__) + " object")

        # Check
        nvalues = len(x_values)
        if nvalues == 0: warnings.warn("The resulting curve will have no points")

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def __radd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self.__add__(other)

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        from ..tools import types
        from ..units.unit import get_converted_value

        # Real or integer value
        if types.is_real_or_integer(other):

            # Get real value
            other = float(other)

            # Set columns
            x_values = self.x_array
            y_values = self.y_array - other

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, self.y_unit)

        # Add quantity
        elif types.is_quantity(other):

            # Get value in correct unit
            other = get_converted_value(other, self.y_unit)

            # Set columns
            x_values = self.x_array
            y_values = self.y_array - other

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, self.y_unit)

        # Add other curve
        elif isinstance(other, Curve):

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
            #print("SUBTRACTING CURVES")
            while True:

                # Get the values
                x_a = self.get_value(self.x_name, i, unit=x_unit, add_unit=False)
                x_b = other.get_value(other.x_name, j, unit=x_unit, add_unit=False)

                #print(x_a, x_b, x_a==x_b)

                # Value is the same: subtract
                #if x_a == x_b:
                isclose = np.isclose(x_a, x_b)
                #if not isclose: print("NOT CLOSE:", x_a, x_b)
                if isclose:

                    # Try to get wavelength and distance for unit conversion
                    conversion_info = {}

                    # Set wavelength
                    if isinstance(self, WavelengthCurve): conversion_info["wavelength"] = self.get_wavelength(i)
                    elif isinstance(other, WavelengthCurve): conversion_info["wavelength"] = other.get_wavelength(j)

                    # Set distance
                    if hasattr(self, "distance") and self.distance is not None: conversion_info["distance"] = self.distance
                    elif hasattr(other, "distance") and other.distance is not None: conversion_info["distance"] = other.distance
                    if len(conversion_info) == 0: conversion_info = None

                    # Calculate the difference of the y values
                    result = self.get_value(self.y_name, i, unit=y_unit, add_unit=False, conversion_info=conversion_info) - other.get_value(other.y_name, j, unit=y_unit, add_unit=False, conversion_info=conversion_info)

                    # Add the x value and the new y value
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

        # Invalid type
        else: raise TypeError("Cannot subtract " + self.__class__.__name__ + " and a " + str(type(other).__name__) + " object")

        # Check
        nvalues = len(x_values)
        if nvalues == 0: warnings.warn("The resulting curve will have no points")

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def __rsub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self.__sub__(other) * -1

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        from ..tools import types
        from ..units.quantity import is_scalar, get_scalar

        # Real or integer value
        if types.is_real_or_integer(other):

            # Get real value
            other = float(other)

            # Set columns
            x_values = self.x_array
            y_values = self.y_array * other

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, self.y_unit)

        # Add quantity
        elif types.is_quantity(other):

            # Get value in correct unit
            #other = get_converted_value(other, self.y_unit)
            value = other.value
            unit = other.unit
            new_unit = self.y_unit * unit

            # Set columns
            x_values = self.x_array
            y_values = self.y_array * value

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, new_unit)

        # Add other curve
        elif isinstance(other, Curve):

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

                    # Convert from dimensionless quantity to scalar if necessary
                    if is_scalar(result): result = get_scalar(result)

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

            # Set the new units
            x_unit = self.x_unit
            y_unit = self.y_unit * other.y_unit

            # Set the names and units
            names = (x_name, y_name)
            units = (x_unit, y_unit)

        # Invalid type
        else: raise TypeError("Cannot multiply " + self.__class__.__name__ + " and a " + str(type(other).__name__) + " object")

        # Check
        nvalues = len(x_values)
        if nvalues == 0: warnings.warn("The resulting curve will have no points")

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def __rmul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return self.__mul__(other)

    # -----------------------------------------------------------------

    def __div__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        from ..units.quantity import is_scalar, get_scalar
        from ..tools import types

        # Real or integer value
        if types.is_real_or_integer(other):

            # Get real value
            other = float(other)

            # Set columns
            x_values = self.x_array
            y_values = self.y_array / other

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, self.y_unit)

        # Add quantity
        elif types.is_quantity(other):

            # Get value in correct unit
            value = other.value
            unit = other.unit
            new_unit = self.y_unit / unit

            # Set columns
            x_values = self.x_array
            y_values = self.y_array / value

            # Set the names and units
            names = (self.x_name, self.y_name)
            units = (self.x_unit, new_unit)

        # Add other curve
        elif isinstance(other, Curve):

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

                # Value of x is the same
                if x_a == x_b:

                    result = self.get_value(self.y_name, i) / other.get_value(other.y_name, j)

                    # Convert from dimensionless quantity to scalar if necessary
                    if is_scalar(result): result = get_scalar(result)

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

            # Set units
            x_unit = self.x_unit
            y_unit = self.y_unit / other.y_unit

            # Set the names and units
            names = (x_name, y_name)
            units = (x_unit, y_unit)

        # Invalid type
        else: raise TypeError("Cannot multiply " + self.__class__.__name__ + " and a " + str(type(other).__name__) + " object")

        # Check
        nvalues = len(x_values)
        if nvalues == 0: warnings.warn("The resulting curve will have no points")

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    __truediv__ = __div__

    # -----------------------------------------------------------------

    def extrapolate_from(self, from_x, regression_from_x, xlog=False, ylog=False, degree=1, replace_nan=None):

        """
        This function ...
        :param from_x:
        :param regression_from_x:
        :param xlog:
        :param ylog:
        :param degree: 1 is default, means linear regression
        :param replace_nan:
        :return:
        """

        # Get the monotonic decreasing or increasing part (or fittable part in general)
        main_part = self.splice(regression_from_x, from_x)

        # Get x and y values for this part
        x_values = main_part.get_x(unit=main_part.x_unit, asarray=True)
        y_values = main_part.get_y(unit=main_part.y_unit, asarray=True)

        # To log?
        if xlog: x_values = np.log10(x_values)
        if ylog: y_values = np.log10(y_values)

        # Plot
        #from ...magic.tools import plotting
        #plotting.plot_xy(x_values, y_values)

        # Polynomial regression
        polyfun = np.poly1d(np.polyfit(x_values, y_values, degree))

        # Extrapolate to larger x
        larger_x, indices = self.get_x_splice(x_min=from_x, return_indices=True)
        if xlog: larger_x = np.log10(larger_x)
        larger_y = polyfun(larger_x) # is array

        # Correct offset
        original_first = self.get_value(self.y_name, indices[0], add_unit=False)
        if ylog: original_first = np.log10(original_first)
        longer_y_first = larger_y[0]
        diff = longer_y_first - original_first
        larger_y = larger_y - diff

        # Plot
        #plotting.plot_xy(larger_x, larger_y)

        # Replace with extrapolated data
        for index, y in zip(indices, larger_y):

            # Set
            if ylog: y = 10**y

            # Replace Nan?
            if replace_nan is not None and np.isnan(y): y = replace_nan

            # Set value
            self.y_data[index] = y

    # -----------------------------------------------------------------

    def extrapolated_from(self, from_x, regression_from_x, xlog=False, ylog=False, replace_nan=None):

        """
        This function ...
        :param from_x:
        :param regression_from_x:
        :param xlog:
        :param ylog:
        :param replace_nan:
        :return:
        """

        # Get copy
        sed = self.copy()

        # Extrapolate
        sed.extrapolate_from(from_x, regression_from_x, xlog=xlog, ylog=ylog, replace_nan=replace_nan)

        # Return
        return sed

    # -----------------------------------------------------------------

    @property
    def x_values(self):
        return self.get_x(unit=self.x_unit, asarray=True)

    # -----------------------------------------------------------------

    @property
    def y_values(self):
        return self.get_y(unit=self.y_unit, asarray=True)

    # -----------------------------------------------------------------

    def get_x_mask(self, lower=None, upper=None, close=None, rtol=1.e-5, atol=1.e-8):

        """
        This function returns a mask of the rows where x meets certain conditions
        :param lower:
        :param upper:
        :param close:
        :param rtol
        :param atol:
        :return:
        """

        # Initialize mask
        mask = np.ones_like(self.x_array, dtype=bool)

        # Lower limit
        if lower is not None:

            if hasattr(lower, "unit"): lower = lower.to(self.x_unit).value
            mask *= self.x_array > lower

        # Upper limit
        if upper is not None:

            if hasattr(upper, "unit"): upper = upper.to(self.x_unit).value
            mask *= self.x_array < upper

        # Close
        if close is not None:

            if hasattr(close, "unit"): close = close.to(self.x_unit).value
            mask *= np.isclose(self.x_array, close, rtol=rtol, atol=atol)

        # Return
        return mask

    # -----------------------------------------------------------------

    def get_y_mask(self, lower=None, upper=None, close=None, rtol=1.e-5, atol=1.e-8):

        """
        This function returns a mask of the rows where y meets certain conditions
        :param lower:
        :param upper:
        :param close:
        :param rtol:
        :param atol:
        :return:
        """

        # Initialize mask
        mask = np.ones_like(self.y_array, dtype=bool)

        # Lower limit
        if lower is not None:

            if hasattr(lower, "unit"): lower = lower.to(self.y_unit).value
            mask *= self.y_array > lower

        # Upper limit
        if upper is not None:

            if hasattr(upper, "unit"): upper = upper.to(self.y_unit).value
            mask *= self.y_array < upper

        # Close
        if close is not None:

            if hasattr(close, "unit"): close = close.to(self.y_unit).value
            mask *= np.isclose(self.y_array, close, rtol=rtol, atol=atol)

        # Return
        return mask

    # -----------------------------------------------------------------

    def get_zero_mask(self):
        return self.y_values == 0

    # -----------------------------------------------------------------

    def get_zero_indices(self):
        return np.where(self.get_zero_mask())[0]

    # -----------------------------------------------------------------

    def get_nonzero_mask(self):
        return self.y_values != 0

    # -----------------------------------------------------------------

    def get_nonzero_indices(self):
        return np.nonzero(self.y_values)

    # -----------------------------------------------------------------

    def get_positive_mask(self):
        return self.y_values > 0

    # -----------------------------------------------------------------

    def get_positive_indices(self):
        return np.where(self.get_positive_mask())[0]

    # -----------------------------------------------------------------

    def get_negative_mask(self):
        return self.y_values < 0

    # -----------------------------------------------------------------

    def get_negative_indices(self):
        return np.where(self.get_negative_mask())[0]

    # -----------------------------------------------------------------

    def get_nonnegative_mask(self):
        return self.y_values >= 0

    # -----------------------------------------------------------------

    def get_nonnegative_indices(self):
        return np.where(self.get_nonnegative_mask())[0]

    # -----------------------------------------------------------------

    def get_nonpositive_mask(self):
        return self.y_values <= 0

    # -----------------------------------------------------------------

    def get_nonpositive_indices(self):
        return np.where(self.get_nonpositive_mask())[0]

    # -----------------------------------------------------------------

    def stripped_zeroes(self):

        """
        This function ...
        :return:
        """

        # Get indices
        indices = self.get_nonzero_indices()
        first = self.get_value(self.x_name, indices[0])
        last = self.get_value(self.x_name, indices[-1])

        # Return the new curve
        return self.splice(first, last, include_min=True, include_max=True)

    # -----------------------------------------------------------------

    def stripped_negatives_and_zeroes(self):

        """
        Thisf unction ...
        :return:
        """

        # Get indices
        indices = self.get_positive_indices()
        first = self.get_value(self.x_name, indices[0])
        last = self.get_value(self.x_name, indices[-1])

        # Return the new curve
        return self.splice(first, last, include_min=True, include_max=True)

    # -----------------------------------------------------------------

    def set_negatives_to_zero(self):

        """
        This function ...
        :return:
        """

        # Get indices
        indices = self.get_negative_indices()

        # Replace
        for index in indices: self[self.y_name][index] = 0.0

    # -----------------------------------------------------------------

    def replace_zeros_by_lowest(self, factor=1):

        """
        This function ...
        :param factor:
        :return:
        """

        # Get indices of zeroes
        indices = self.get_zero_indices()

        # Get minimum value (except for zero)
        min_value = self.get_min_y_value(self.y_unit, add_unit=False, ignore_zero=True)
        value = min_value * factor

        # Replace
        for index in indices: self[self.y_name][index] = value

    # -----------------------------------------------------------------

    def extended_to_right(self, to_x, logscale=False, value=0., points=10):

        """
        This function ...
        :param to_x:
        :param logscale:
        :param value:
        :param points: CAN BE INTEGER OR LIST OF VALUES TO PICK FROM
        :return:
        """

        # Convert to_x
        if self.x_unit is not None: to_x = to_x.to(self.x_unit).value
        elif hasattr(to_x, "unit"): raise RuntimeError("Unexpected value: has unit")

        # Get values, as lists
        x_values = self.get_x(unit=self.x_unit, add_unit=False)
        y_values = self.get_y(unit=self.y_unit, add_unit=False)

        # Get units and names
        names = [self.x_name, self.y_name]
        units = [self.x_unit, self.y_unit]

        # Add values if necessary
        if to_x > x_values[-1]:

            # Set range of new values
            new_x_range = RealRange(x_values[-1], to_x, inclusive=False)

            # Get new x points
            if isinstance(points, int): new_x_values = generate_values(new_x_range, points, logscale=logscale)
            else:
                if self.x_unit is not None: points = [point.to(self.x_unit).value for point in points]
                new_x_values = generate_values(new_x_range, pick_from=points)
            npoints = len(new_x_values)

            # Add new values
            x_values.extend(new_x_values)
            y_values.extend([value] * npoints)

            #for new_x_value in new_x_values:
                #print("Adding " + str(new_x_value) + ", " + str(value) + " ...")
                #x_values.append(new_x_value)
                #y_values.append(value)

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def extended_to_left(self, to_x, logscale=False, value=0., points=10):

        """
        This function ...
        :param to_x:
        :param logscale:
        :param value:
        :param points:
        :return:
        """

        # Convert to_x
        if self.x_unit is not None: to_x = to_x.to(self.x_unit).value
        elif hasattr(to_x, "unit"): raise RuntimeError("Unexpected value: has unit")

        # Get values, as lists
        x_values = self.get_x(unit=self.x_unit, add_unit=False)
        y_values = self.get_y(unit=self.y_unit, add_unit=False)

        # Get units and names
        names = [self.x_name, self.y_name]
        units = [self.x_unit, self.y_unit]

        # Add values if necessary
        if to_x < x_values[0]:

            # Set range of new values
            new_x_range = RealRange(to_x, x_values[0], inclusive=False)

            # Get new x points
            if isinstance(points, int): new_x_values = generate_values(new_x_range, points, logscale=logscale)
            else:
                if self.x_unit is not None: points = [point.to(self.x_unit).value for point in points]
                new_x_values = generate_values(new_x_range, pick_from=points)
            npoints = len(new_x_values)

            # Add new values
            x_values = new_x_values + x_values
            y_values = [value] * npoints + y_values

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

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

    @classmethod
    def from_wavelengths_and_values(cls, name, wavelengths, values, wavelength_unit="micron", value_unit=None, description=None):

        """
        This function ...
        :param name:
        :param wavelengths:
        :param values:
        :param wavelength_unit:
        :param value_unit:
        :param description:
        :return:
        """

        # Determine the units
        if value_unit is None: value_unit = get_common_unit(values)

        # Create the curve
        curve = cls(y_name=name, y_unit=value_unit, y_description=description, wavelength_unit=wavelength_unit)

        # Add the points
        for wavelength, value in zip(wavelengths, values):
            curve.add_point(wavelength, value)

        # Return the curve
        return curve

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

        # X name?
        if self.y_name in self.colnames: return self.y_name

        # Search
        for index in reversed(range(len(self.colnames))):
            name = self.colnames[index]
            if name == "Error+" or name == "Error-": continue
            return self.colnames[index]

    # -----------------------------------------------------------------

    @property
    def wavelength_name(self):
        return self.x_name

    # -----------------------------------------------------------------

    @property
    def unit(self):
        return self.y_unit

    # -----------------------------------------------------------------

    @property
    def has_unit(self):
        return self.unit is not None

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):
        return self.x_unit

    # -----------------------------------------------------------------

    @property
    def has_wavelength_unit(self):
        return self.wavelength_unit is not None

    # -----------------------------------------------------------------

    def wavelength_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value(self.wavelength_name, index)

    # -----------------------------------------------------------------

    def get_wavelength(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value(self.wavelength_name, index)

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
                             interpolate=True, conversion_info=None, distance=None):

        """
        This function ...
        :param wavelength:
        :param unit:
        :param add_unit:
        :param density:
        :param brightness:
        :param interpolate:
        :param conversion_info:
        :param distance:
        :return:
        """

        if interpolate:
            interpolated = interp1d(self.wavelengths(unit="micron", asarray=True), self.values(self.unit, asarray=True), kind='linear')
            value = interpolated(wavelength.to("micron").value)
            if self.has_unit: value = value * self.unit
        else:
            from ..tools import sequences
            index = sequences.find_closest_index(self.wavelengths(unit="micron", add_unit=False), wavelength.to("micron").value)
            value = self[self.value_name][index]
            if self.has_unit: value = value * self.unit

        # Convert unit if necessary
        if unit is not None:

            # Create conversion info
            if conversion_info is None: conversion_info = dict()
            conversion_info["wavelength"] = wavelength
            if distance is not None: conversion_info["distance"] = distance
            elif self.distance is not None: conversion_info["distance"] = self.distance

            unit = u(unit, density=density, brightness=brightness)
            value = value.to(unit, **conversion_info)

        # Remove unit if requested
        if not add_unit and self.has_unit: value = value.value

        # Return
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
        return self.get_value(self.wavelength_name, 0)

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):
        return self.get_value(self.wavelength_name, -1)

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
        if asarray: return arrays.plain_array(self[self.wavelength_name], unit=unit, array_unit=self.column_unit(self.wavelength_name), mask=mask)
        else: return arrays.array_as_list(self[self.wavelength_name], unit=unit, add_unit=add_unit, array_unit=self.column_unit(self.wavelength_name), mask=mask)

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
        if asarray: return arrays.plain_array(self[self.wavelength_name], unit=unit, array_unit=self.column_unit(self.wavelength_name), equivalencies=spectral(), mask=mask)
        else: return arrays.array_as_list(self[self.wavelength_name], unit=unit, add_unit=add_unit, array_unit=self.column_unit(self.wavelength_name), equivalencies=spectral(), mask=mask)

    # -----------------------------------------------------------------

    def values(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False,
               min_wavelength=None, max_wavelength=None, distance=None):

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
        :param distance:
        :return:
        """

        #print("VALUE NAME", self.value_name)

        # Create mask
        if min_wavelength is not None or max_wavelength is not None: mask = self.wavelengths_mask(min_wavelength=min_wavelength, max_wavelength=max_wavelength)
        else: mask = None

        # Create conversion info
        if conversion_info is None: conversion_info = dict()
        conversion_info["wavelengths"] = self.wavelengths()
        if distance is not None: conversion_info["distance"] = distance
        elif self.distance is not None: conversion_info["distance"] = self.distance

        #print(self.colnames)
        #print(self.column_info)

        #print(self.value_name)
        #print(unit, self.column_unit(self.value_name))

        # Create and return
        if asarray: return arrays.plain_array(self[self.value_name], unit=unit, array_unit=self.column_unit(self.value_name),
                                              conversion_info=conversion_info, density=density, brightness=brightness,
                                              mask=mask)
        else: return arrays.array_as_list(self[self.value_name], unit=unit, add_unit=add_unit,
                                          array_unit=self.column_unit(self.value_name), conversion_info=conversion_info,
                                          density=density, brightness=brightness, mask=mask)

    # -----------------------------------------------------------------

    def get_positive_wavelength_mask(self, lower=None, upper=None):

        """
        This function ...
        :param lower:
        :param upper:
        :return:
        """

        # Get mask
        return self.get_positive_mask() * self.get_x_mask(lower=lower, upper=upper)

    # -----------------------------------------------------------------

    def get_positive_wavelength_indices(self, lower=None, upper=None):

        """
        This function ...
        :param lower:
        :param upper:
        :return:
        """

        return np.where(self.get_positive_wavelength_mask(lower=lower, upper=upper))[0]

    # -----------------------------------------------------------------

    def get_min_positive_wavelength_index(self, lower=None, upper=None):

        """
        This function returns the minimum wavelength for which the photometry is positive, with additional (optional) conditions
        :return:
        """

        # Return the min index
        return min(self.get_positive_wavelength_indices(lower=lower, upper=upper))

    # -----------------------------------------------------------------

    def get_min_positive_wavelength(self, lower=None, upper=None):

        """
        This function ...
        :param lower:
        :param upper:
        :return:
        """

        # Get the index
        index = self.get_min_positive_wavelength_index(lower=lower, upper=upper)

        # Return the wavelength
        return self.get_wavelength(index)

    # -----------------------------------------------------------------

    def get_max_positive_wavelength_index(self, lower=None, upper=None):

        """
        This function returns the maximum wavelength for which the photometry is positive, with additional (optional) conditions
        :return:
        """

        # Return the max index
        return max(self.get_positive_wavelength_indices(lower=lower, upper=upper))

    # -----------------------------------------------------------------

    def get_max_positive_wavelength(self, lower=None, upper=None):
        
        """
        Thisf unction ...
        :param lower: 
        :param upper: 
        :return: 
        """
        
        # Get the index
        index = self.get_max_positive_wavelength_index(lower=lower, upper=upper)
        
        # Return the wavelength
        return self.get_wavelength(index)

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

        # Get options
        extra_columns = kwargs.pop("extra_columns", [])

        # Call the constructor of the base class
        super(FilterCurve, self).__init__(*args, **kwargs)

        # Set column info
        self.column_info.insert(0, ("Band", str, None, "band"))
        self.column_info.insert(0, ("Instrument", str, None, "instrument"))
        self.column_info.insert(0, ("Observatory", str, None, "observatory"))

        # Extra columns
        for coldef in extra_columns: self.column_info.append(coldef)

    # -----------------------------------------------------------------

    def add_point(self, fltr, value, extra=None, conversion_info=None):

        """
        This function ...
        :param fltr:
        :param value:
        :param extra:
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

        # Add extra?
        if extra is not None: values.extend(extra)

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
        return arrays.array_as_list(self["Instrument"])

    # -----------------------------------------------------------------

    def bands(self):
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

    def has_broad_band(self):
        return self.nbroad_band_filters > 0

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

    def has_narrow_band(self):
        return self.nnarrow_band_filters > 0

    # -----------------------------------------------------------------

    def broad_band_filters(self):
        return [isinstance(fltr, BroadBandFilter) for fltr in self.filters()]

    # -----------------------------------------------------------------

    @property
    def nbroad_band_filters(self):
        return len(self.broad_band_filters())

    # -----------------------------------------------------------------

    def narrow_band_filters(self):
        return [isinstance(fltr, NarrowBandFilter) for fltr in self.filters()]

    # -----------------------------------------------------------------

    @property
    def nnarrow_band_filters(self):
        return len(self.narrow_band_filters())

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

    def for_filters(self, filters):

        """
        This function ...
        :param filters:
        :return:
        """

        # Make a copy of this SED
        new = self.copy()

        # Loop over the rows, remove the row if it does not correspond to a broad band filter
        for index in reversed(range(len(self))):

            # Remove if filter not in list
            fltr = self.get_filter(index)
            if fltr not in filters: new.remove_row(index)

        # Return the new SED
        return new

# -----------------------------------------------------------------

def generate_values(value_range, npoints=None, logscale=False, pick_from=None):

    """
    This function ...
    :param value_range:
    :param npoints:
    :param logscale:
    :param pick_from:
    :return:
    """

    # Pick from list, within range
    if pick_from is not None: return value_range.values_in_range(pick_from)
    else:
        if npoints is None: raise ValueError("Number of points must be defined")
        if logscale: return value_range.log(npoints)
        else: return value_range.linear(npoints)

# -----------------------------------------------------------------
