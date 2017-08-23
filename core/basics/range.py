#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.range Contains the IntegerRange, RealRange and QuantityRange classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import bisect
import numpy as np

# Import astronomical modules
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from ..units.parsing import parse_unit as u
from ..units.stringify import represent_unit as ru
from ..units.parsing import parse_quantity
from ..tools import types

# -----------------------------------------------------------------

# TODO: what to do with the combination of inclusive=False and invert=True ??
# Define inclusive for minimum value and maximum value seperately??

# -----------------------------------------------------------------

def range_around(value, rel_min, rel_max):

    """
    This function ...
    :param value:
    :param rel_min:
    :param rel_max:
    :return:
    """

    if types.is_integer_type(value): return IntegerRange.around(value, rel_min, rel_max)
    elif types.is_real_type(value): return RealRange.around(value, rel_min, rel_max)
    elif types.is_quantity(value): return QuantityRange.around(value, rel_min, rel_max)
    else: raise ValueError("Value has unknown type '" + str(type(value)) + "'")

# -----------------------------------------------------------------

class Range(object):
    
    """
    This class ...
    """
    
    def __init__(self, min_value, max_value, inclusive=True, invert=False, rearrange=False):
        
        """
        The constructor ...
        :param min_value:
        :param max_value:
        :param inclusive:
        """

        self._min = min_value
        self._max = max_value

        # Check whether min and max are in the right order
        if self._max < self._min:
            if rearrange:
                oldmin = self._min
                self._min = self._max
                self._max = oldmin
            else: raise ValueError("Minimum value should be lower than maximum value")

        self.inclusive = inclusive
        self.invert = invert

    # -----------------------------------------------------------------

    @classmethod
    def infinitesimal(cls, value):

        """
        This function ...
        :param value:
        :return:
        """

        return cls(value, value)

    # -----------------------------------------------------------------

    @classmethod
    def around(cls, value, rel_min, rel_max, inclusive=True, invert=False):

        """
        This function ...
        :param value:
        :param rel_min:
        :param rel_max:
        :param inclusive:
        :param invert:
        :return:
        """

        return cls(value * rel_min, value * rel_max, inclusive, invert)

    # -----------------------------------------------------------------

    @classmethod
    def around_magnitudes(cls, value, min_magnitude, max_magnitude, inclusive=True, invert=False):

        """
        This function ...
        :param value:
        :param min_magnitude:
        :param max_magnitude:
        :param inclusive:
        :param invert:
        :return:
        """

        # Determine rel min and max
        # -2 becomes 0.01
        # 2 becomes 100
        rel_min = 10**float(min_magnitude)
        rel_max = 10**float(max_magnitude)

        # Create and return
        return cls.around(value, rel_min, rel_max, inclusive=inclusive, invert=invert)

    # -----------------------------------------------------------------

    @classmethod
    def around_magnitude(cls, value, magnitude, inclusive=True, invert=False):

        """
        Symmetrical around -magnitude and +magnitude
        :param value:
        :param magnitude:
        :param inclusive:
        :param invert:
        :return:
        """

        return cls.around_magnitudes(value, -magnitude, magnitude, inclusive=inclusive, invert=invert)

    # -----------------------------------------------------------------

    @property
    def min(self):

        """
        This function ...
        :return:
        """

        return self._min

    # -----------------------------------------------------------------

    @property
    def max(self):

        """
        This function ...
        :return:
        """

        return self._max

    # -----------------------------------------------------------------

    @property
    def log_min(self):

        """
        This function ...
        :return:
        """

        return np.log10(self._min)

    # -----------------------------------------------------------------

    @property
    def log_max(self):

        """
        This function ...
        :return:
        """

        return np.log10(self._max)

    # -----------------------------------------------------------------

    @property
    def mean(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.min + self.max)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return self.mean

    # -----------------------------------------------------------------

    def as_tuple(self):

        """
        This function ...
        :return:
        """

        return (self.min, self.max)

    # -----------------------------------------------------------------

    def linear(self, npoints, as_list=False, fancy=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :param fancy:
        """

        # If only one point, return the arithmetic mean
        if npoints == 1:
            if as_list: return [0.5 * (self._min + self._max)]
            else: return np.array([0.5 * (self._min + self._max)])

        if fancy:
            span = self.max - self.min
            maxnpoints = npoints
            ticksize = best_tick(span, maxnpoints)
            fancy_min = round_to_1(self.min)
            values = np.arange(fancy_min, self.max, step=ticksize)
        else: values = np.linspace(self._min, self._max, npoints, endpoint=self.inclusive)
        if self.invert: values = np.flipud(values)

        values = sorted(list(set(values)))

        if as_list: return list(values)
        else: return values

    # -----------------------------------------------------------------

    def log(self, npoints, as_list=False, fancy=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :param fancy:
        :return:
        """

        # If only one point, return the geometric mean
        if npoints == 1:
            if as_list: return [np.sqrt(self._min * self._max)]
            else: return np.array([np.sqrt(self._min * self._max)])

        if fancy:

            values = []
            new_npoints = npoints

            while len(values) < 0.8 * npoints:

                #new_npoints = int(new_npoints * 1.2)
                #values = [round_to_1(self.min)] * new_npoints

                if len(values) > 0:
                    factor = npoints / len(values)
                    new_npoints = int(new_npoints * factor)
                values = [round_to_1(self.min)] * new_npoints

                #print(new_npoints)
                ticksize = best_tick_log(self.max / self.min, new_npoints)

                for i in range(1, new_npoints):
                    new_value = values[i-1] * ticksize
                    if new_value > self.max: break
                    else: values[i] = new_value
                values = np.array(values)
                values = sorted(list(set(values)))

                #print(values)

                # DOESNT WORK, JUST BREAK AFTER THE FIRST
                break

        else: values = np.logspace(self.log_min, self.log_max, npoints, endpoint=self.inclusive)
        if self.invert: values = np.flipud(values)

        if as_list: return list(values)
        else: return values

    # -----------------------------------------------------------------

    def sqrt(self, npoints, as_list=False, fancy=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :param fancy:
        :return:
        """

        if fancy: raise NotImplementedError("Not implemented")

        width = self._max - self._min
        normalized = np.linspace(0.0, 1.0, npoints, endpoint=self.inclusive)
        values = self._min + normalized * width
        if self.invert: values = np.flipud(values)

        values = sorted(list(set(values)))

        if as_list: return list(values)
        else: return values

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return str(self.min) + " > " + str(self.max)

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return repr(self.min) + " > " + repr(self.max)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from ..tools import types

        if isinstance(value, Quantity): return QuantityRange(self.min * value, self.max * value)
        elif types.is_real_type(value): return RealRange(self.min * value, self.max * value)
        elif types.is_integer_type(value) and isinstance(self, IntegerRange):
            return IntegerRange(self.min * value, self.max * value)
        elif types.is_integer_type(value): return RealRange(self.min * value, self.max * value)
        else: raise ValueError("Value must be Quantity, real or integer value")

    # -----------------------------------------------------------------

    def __rmul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__mul__(value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from ..tools import types

        if isinstance(value, Quantity): return QuantityRange(self.min / value, self.max / value)
        elif types.is_real_type(value) or types.is_integer_type(value):
            return RealRange(self.min / value, self.max / value)
        else: raise ValueError("Value must be Quantity, real or integer value")

    # -----------------------------------------------------------------

    def __add__(self, value):

        """
        This fucntion ...
        :param value: 
        :return: 
        """

        return self.__class__(self.min + value, self.max + value)

    # -----------------------------------------------------------------

    def __iadd__(self, value):

        """
        This function ...
        :param value: 
        :return: 
        """

        self._min += value
        self._max += value
        return self

    # -----------------------------------------------------------------

    def __sub__(self, value):

        """
        This function ...
        :param value: 
        :return: 
        """

        return self.__class__(self.min - value, self.max - value)

    # -----------------------------------------------------------------

    def __isub__(self, value):

        """
        This fucntion ...
        :param value: 
        :return: 
        """

        self._min -= value
        self._max -= value
        return self

    # -----------------------------------------------------------------

    def adjust(self, value_or_range):

        """
        This function ...
        :param value_or_range:
        :return:
        """

        # Range
        if isinstance(value_or_range, Range):

            if value_or_range.min < self.min: self._min = value_or_range.min
            if value_or_range.max > self.max: self._max = value_or_range.max

        # Value
        else:

            # Adjust minimal and maximal speedup
            if value_or_range < self.min: self._min = value_or_range
            elif value_or_range > self.max: self._max = value_or_range

    # -----------------------------------------------------------------

    def adjust_inwards(self, other_range):

        """
        This function ...
        :param other_range:
        :return:
        """

        if other_range.min > self.min: self._min = other_range.min
        if other_range.max < self.max: self._max = other_range.max

    # -----------------------------------------------------------------

    @property
    def span(self):

        """
        This function ...
        :return:
        """

        return self.max - self.min

    # -----------------------------------------------------------------

    def __contains__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.min <= value <= self.max

# -----------------------------------------------------------------

class IntegerRange(Range):

    """
    This class ...
    """

    def __init__(self, min_value, max_value, inclusive=True, invert=False, rearrange=False):

        """
        This function ...
        :param min_value:
        :param max_value:
        :param inclusive:
        :param invert:
        :param rearrange:
        """

        assert isinstance(min_value, int)
        assert isinstance(max_value, int)

        # Call the constructor of the base class
        super(IntegerRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert, rearrange=rearrange)

    # -----------------------------------------------------------------

    def linear(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(IntegerRange, self).linear(npoints)
        integers = sorted(list(set(map(int, real))))

        if as_list: return integers
        else: return np.array(integers)

    # -----------------------------------------------------------------

    def log(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(IntegerRange, self).log(npoints)
        integers = sorted(list(set(map(int, real))))

        if as_list: return integers
        else: return np.array(integers)

    # -----------------------------------------------------------------

    def sqrt(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(IntegerRange, self).sqrt(npoints)
        integers = sorted(list(set(map(int, real))))

        if as_list: return integers
        else: return np.array(integers)

# -----------------------------------------------------------------

class RealRange(Range):

    """
    This class ...
    """

    def __init__(self, min_value, max_value, inclusive=True, invert=False, rearrange=False):

        """
        The constructor ...
        :param min_value:
        :param max_value:
        :param inclusive:
        :param invert:
        :param rearrange:
        """

        if min_value is None: min_value = float("-inf")
        if max_value is None: max_value = float("+inf")

        if isinstance(min_value, int): min_value = float(min_value)
        if isinstance(max_value, int): max_value = float(max_value)

        assert isinstance(min_value, float)
        assert isinstance(max_value, float)

        # Call the constructor of the base class
        super(RealRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert, rearrange=rearrange)

    # -----------------------------------------------------------------

    @classmethod
    def infinity(cls):

        """
        This function ...
        :return:
        """

        return cls(-float("inf"), float("inf"))

    # -----------------------------------------------------------------

    @classmethod
    def zero(cls):

        """
        This function ...
        :return:
        """

        return cls(-0., 0.)

# -----------------------------------------------------------------

class QuantityRange(Range):

    """
    This function ...
    """

    def __init__(self, min_value, max_value, unit=None, inclusive=True, invert=False, rearrange=False):

        """
        This function ...
        :param min_value:
        :param max_value:
        :param unit:
        :param inclusive:
        :param invert:
        :param rearrange:
        """

        # Convert strings
        if types.is_string_type(min_value): min_value = parse_quantity(min_value)
        if types.is_string_type(max_value): max_value = parse_quantity(max_value)

        # Convert everything so that min_value and max_value are floats in the same unit, and so that 'unit' is the corresponding Unit
        min_is_quantity = hasattr(min_value, "unit")
        max_is_quantity = hasattr(max_value, "unit")

        if min_is_quantity and max_is_quantity:

            unit = min_value.unit
            min_value = min_value.value
            max_value = max_value.to(unit).value

        elif (not min_is_quantity) and (not max_is_quantity):

            if unit is None: raise ValueError("Unit must be specified if min_value and max_value are not quantities")
            elif types.is_string_type(unit): unit = u(unit)

        else: raise ValueError("min_value and max_value must be either both quantities or both floats (with unit specified seperately)")

        # Call the constructor of the base class
        super(QuantityRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert, rearrange=rearrange)

        # Set the unit
        self.unit = unit

    # -----------------------------------------------------------------

    @property
    def min(self):

        """
        This function ...
        :return:
        """

        return self._min * self.unit

    # -----------------------------------------------------------------

    @property
    def max(self):

        """
        This function ...
        :return:
        """

        return self._max * self.unit

    # -----------------------------------------------------------------

    @min.setter
    def min(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self._min = value.to(self.unit).value

    # -----------------------------------------------------------------

    @max.setter
    def max(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self._max = value.to(self.unit).value

    # -----------------------------------------------------------------

    def adjust(self, value_or_range):

        """
        This function ...
        :param value_or_range:
        :return:
        """

        # Range
        if isinstance(value_or_range, Range):

            if value_or_range.min < self.min: self._min = value_or_range.min.to(self.unit).value
            if value_or_range.max > self.max: self._max = value_or_range.max.to(self.unit).value

        # Value
        else:

            # Adjust minimal and maximal speedup
            if value_or_range < self.min: self._min = value_or_range.to(self.unit).value
            elif value_or_range > self.max: self._max = value_or_range.to(self.unit).value

    # -----------------------------------------------------------------

    def adjust_inwards(self, other_range):

        """
        This function ...
        :param other_range:
        :return:
        """

        if other_range.min > self.min: self._min = other_range.min.to(self.unit).value
        if other_range.max < self.max: self._max = other_range.max.to(self.unit).value

    # -----------------------------------------------------------------

    def linear(self, npoints, as_list=False, add_unit=True):

        """
        This function ...
        :param npoints:
        :param as_list:
        :param add_unit:
        :return:
        """

        real = super(QuantityRange, self).linear(npoints)

        if as_list:

            if add_unit:
                result = []
                for num in real: result.append(num * self.unit)
            else: result = list(real)

        else:

            if add_unit: result = real * self.unit
            else: result = real

        # Return the result (list or quantity)
        return result

    # -----------------------------------------------------------------

    def log(self, npoints, as_list=False, add_unit=True):

        """
        This function ...
        :param npoints:
        :param as_list:
        :param add_unit:
        :return:
        """

        real = super(QuantityRange, self).log(npoints)

        if as_list:

            if add_unit:
                result = []
                for num in real: result.append(num * self.unit)
            else: result = list(real)

        else:

            if add_unit: result = real * self.unit
            else: result = real

        # Return the result (list or quantity)
        return result

    # -----------------------------------------------------------------

    def sqrt(self, npoints, as_list=False, add_unit=True):

        """
        This function ...
        :param npoints:
        :param as_list:
        :param add_unit:
        :return:
        """

        real = super(QuantityRange, self).sqrt(npoints)

        if as_list:

            if add_unit:

                result = []
                for num in real: result.append(num * self.unit)

            else: result = list(real)

        else:

            if add_unit: result = real * self.unit
            else: result = real

        # Return the result (list or quantity)
        return result

    # -----------------------------------------------------------------

    def to(self, unit, equivalencies=None):

        """
        This function ...
        :param unit:
        :param equivalencies:
        :return:
        """

        # Rearrange is allowed because it is possible that the old min becomes the new max (eg. frequency -> wavelength)
        return QuantityRange(self.min.to(unit, equivalencies=equivalencies), self.max.to(unit, equivalencies=equivalencies), rearrange=True)

    # -----------------------------------------------------------------

    @property
    def value(self):

        """
        Thisf unction ...
        :return:
        """

        return RealRange(self.min.value, self.max.value)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return str(self.min.value) + " " + ru(self.min.unit) + " > " + str(self.max.value) + " " + ru(self.max.unit)

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return repr(self.min.value) + " " + ru(self.min.unit) + " > " + repr(self.max.value) + " " + ru(self.max.unit)

# -----------------------------------------------------------------

def zip_linear(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    npoints = kwargs.pop("npoints")

    temp = []
    for arg in args:
        if isinstance(arg, QuantityRange): temp.append(arg.linear(npoints, as_list=True))
        else: temp.append(arg.linear(npoints))

    # Zip
    result = zip(*temp)
    return result

# -----------------------------------------------------------------

def zip_log(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    npoints = kwargs.pop("npoints")

    temp = []
    for arg in args:
        if isinstance(arg, QuantityRange): temp.append(arg.log(npoints, as_list=True))
        else: temp.append(arg.log(npoints))

    # Zip
    result = zip(*temp)
    return result

# -----------------------------------------------------------------

def zip_sqrt(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    npoints = kwargs.pop("npoints")

    temp = []
    for arg in args:
        if isinstance(arg, QuantityRange): temp.append(arg.sqrt(npoints, as_list=True))
        else: temp.append(arg.sqrt(npoints))

    # Zip
    result = zip(*temp)
    return result

# -----------------------------------------------------------------

def best_tick(largest, max_nticks):

    """
    This function ...
    :param largest:
    :param max_nticks:
    :return:
    """

    minimum = largest / max_nticks
    magnitude = 10 ** math.floor(math.log(minimum, 10))
    residual = minimum / magnitude
    if residual > 5:
        tick = 10 * magnitude
    elif residual > 2:
        tick = 5 * magnitude
    elif residual > 1:
        tick = 2 * magnitude
    else:
        tick = magnitude
    return tick

# -----------------------------------------------------------------

def best_tick_log(largest, maxnticks):

    """
    This function ...
    :param largest:
    :param maxnticks:
    :return:
    """

    minimum_ticksize = math.pow(largest, 1./maxnticks)
    magnitude = math.floor(minimum_ticksize)
    residual = math.pow(minimum_ticksize, 1./magnitude)
    if residual > 5:
        tick = 10*magnitude
    elif residual > 2:
        tick = 5**magnitude
    elif residual > 1:
        tick = 2**magnitude
    else: tick = magnitude

    # Return
    return tick

# -----------------------------------------------------------------

def best_tick2(largest, max_nticks):

    """
    This function ...
    :param largest:
    :param max_nticks:
    :return:
    """

    minimum = largest / max_nticks
    magnitude = 10 ** math.floor(math.log(minimum, 10))
    residual = minimum / magnitude
    # this table must begin with 1 and end with 10
    table = [1, 1.5, 2, 3, 5, 7, 10]
    tick = table[bisect.bisect_right(table, residual)] if residual < 10 else 10
    return tick * magnitude

# -----------------------------------------------------------------

def round_to_1(x):

    """
    This function ...
    :param x:
    :return:
    """

    return round(x, -int(math.floor(math.log10(abs(x)))))

# -----------------------------------------------------------------
