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

        # Set min and max value
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
    def around_factor(cls, value, factor, inclusive=True, invert=False):

        """
        This function ...
        :param value:
        :param factor:
        :param inclusive:
        :param invert:
        :return:
        """

        return cls.around(value, 1./factor, factor, inclusive=inclusive, invert=invert)

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

    @classmethod
    def within_factor(cls, value, factor, inclusive=True, invert=False):

        """
        This function ...
        :param value:
        :param factor:
        :param inclusive:
        :param invert:
        :return:
        """

        return cls.around_factor(value, math.sqrt(factor), inclusive=inclusive, invert=invert)

    # -----------------------------------------------------------------

    @classmethod
    def within_magnitude(cls, value, magnitude, inclusive=True, invert=False):

        """
        This function ...
        :param value:
        :param magnitude:
        :param inclusive:
        :param invert:
        :return:
        """

        factor = 10**float(magnitude)
        return cls.within_factor(value, factor, inclusive=inclusive, invert=invert)

    # -----------------------------------------------------------------

    @classmethod
    def limits(cls, values, inclusive=True):

        """
        This function ...
        :param values:
        :param inclusive:
        :return:
        """

        min_value = min(values)
        max_value = max(values)
        return cls(min_value, max_value, inclusive=inclusive)

    # -----------------------------------------------------------------

    @property
    def min(self):
        return self._min

    # -----------------------------------------------------------------

    @property
    def max(self):
        return self._max

    # -----------------------------------------------------------------

    @min.setter
    def min(self, value):
        self._min = value

    # -----------------------------------------------------------------

    @max.setter
    def max(self, value):
        self._max = value

    # -----------------------------------------------------------------

    @property
    def log_min(self):
        return np.log10(self._min)

    # -----------------------------------------------------------------

    @property
    def log_max(self):
        return np.log10(self._max)

    # -----------------------------------------------------------------

    @property
    def mean(self):
        return 0.5 * (self.min + self.max)

    # -----------------------------------------------------------------

    @property
    def geometric_mean(self):
        return np.sqrt(self.min * self.max)

    # -----------------------------------------------------------------

    @property
    def center(self):
        return self.mean

    # -----------------------------------------------------------------

    def as_tuple(self):
        return self.min, self.max

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

        values = np.array(sorted(list(set(values))))

        if as_list: return list(values)
        else: return values

    # -----------------------------------------------------------------

    def linear_step(self, step, symmetric=False, center=None):

        """
        This function ...
        :param step:
        :param symmetric:
        :param center:
        :return:
        """

        values = []
        if center is None: center = self.center

        if symmetric:

            below = []
            above = []

            # Add below values
            index = 1
            while True:

                # Calculate the new value
                new_below = center - index * step

                # Stop?
                if new_below < self.min: break

                # Add the new below value
                below.append(new_below)
                index += 1

            # Add above values
            index = 1
            while True:

                # Calcualte the new value
                new_above = center + index * step

                # Stop?
                if self.inclusive:
                    if new_above > self.max: break
                else:
                    if new_above >= self.max: break

                # Add the new above value
                above.append(new_above)
                index += 1

            # Concatenate
            for value in reversed(below): values.append(value)
            values.append(center)
            for value in above: values.append(value)

        else:

            index = 0
            while True:

                # Calculate the new value
                new_value = self.min + index * step

                # Stop?
                if self.inclusive:
                    if new_value > self.max: break
                else:
                    if new_value >= self.max: break

                # Add the new value
                values.append(new_value)
                index += 1

        # Return
        return values

    # -----------------------------------------------------------------

    def log(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        # If only one point, return the geometric mean
        if npoints == 1:
            if as_list: return [np.sqrt(self._min * self._max)]
            else: return np.array([np.sqrt(self._min * self._max)])

        else: values = np.logspace(self.log_min, self.log_max, npoints, endpoint=self.inclusive)
        if self.invert: values = np.flipud(values)

        if as_list: return list(values)
        else: return values

    # -----------------------------------------------------------------

    def sqrt(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

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

    def extend(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        self.min = self.center - factor * self.radius
        self.max = self.center + factor * self.radius

    # -----------------------------------------------------------------

    def extended(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        return self.__class__(self.center - factor * self.radius, self.center + factor * self.radius)

    # -----------------------------------------------------------------

    def compress(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        self.min = self.center - self.radius / factor
        self.max = self.center + self.radius / factor

    # -----------------------------------------------------------------

    def compressed(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        return self.__class__(self.center - self.radius / factor, self.center + self.radius / factor)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        raise RuntimeError("Should be implemented in derived class")

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

        raise RuntimeError("Should be implemented in derived class")

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

    @property
    def radius(self):

        """
        This function ...
        :return:
        """

        return 0.5 * self.span

    # -----------------------------------------------------------------

    # NO: NOT POSSIBLE: __LEN__ ALWAYS RETURNS INTEGERS!
    # def __len__(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return self.span

    # -----------------------------------------------------------------

    def __contains__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if self.inclusive: return self.min <= value <= self.max
        else: return self.min < value < self.max

    # -----------------------------------------------------------------

    def values_in_range(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        in_range = []
        for value in values:
            if value in self: in_range.append(value)
        return in_range

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

        # Check
        if not types.is_integer_type(min_value): raise ValueError("Value must be integer")
        if not types.is_integer_type(max_value): raise ValueError("Value must be integer")

        # Call the constructor of the base class
        super(IntegerRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert, rearrange=rearrange)

    # -----------------------------------------------------------------

    @classmethod
    def limits(cls, values, inclusive=True):

        """
        This function ...
        :param values:
        :param inclusive:
        :return:
        """

        min_value = int(min(values)) # if list of reals
        max_value = int(math.ceil(max(values))) # if list of reals
        return cls(min_value, max_value, inclusive=inclusive)

    # -----------------------------------------------------------------

    @Range.min.setter
    def min(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if not types.is_integer_type(value): raise ValueError("Value must be integer")
        self._min = value

    # -----------------------------------------------------------------

    @Range.max.setter
    def max(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if not types.is_integer_type(value): raise ValueError("Value must be integer")
        self._max = value

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

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from ..tools import types

        if types.is_quantity(value): return QuantityRange(self.min * value, self.max * value)
        elif types.is_real_type(value): return RealRange(self.min * value, self.max * value)
        elif types.is_integer_type(value): return self.__class__(self.min * value, self.max * value)
        else: raise ValueError("Value must be Quantity, real or integer value")

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from ..tools import types

        if types.is_quantity(value): return QuantityRange(self.min / value, self.max / value)
        elif types.is_real_type(value) or types.is_integer_type(value): return RealRange(self.min / float(value), self.max / float(value))
        else: raise ValueError("Value must be Quantity, real or integer value")

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

        if types.is_integer_type(min_value): min_value = float(min_value)
        if types.is_integer_type(max_value): max_value = float(max_value)

        if not types.is_real_type(min_value): raise ValueError("Value must be real")
        if not types.is_real_type(max_value): raise ValueError("Value must be real")

        # Call the constructor of the base class
        super(RealRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert, rearrange=rearrange)

    # -----------------------------------------------------------------

    @classmethod
    def limits(cls, values, inclusive=True):

        """
        This function ...
        :param values:
        :param inclusive:
        :return:
        """

        min_value = float(min(values)) # if list of ints
        max_value = float(max(values)) # if list of ints
        return cls(min_value, max_value, inclusive=inclusive)

    # -----------------------------------------------------------------

    @classmethod
    def infinity(cls):
        return cls(-float("inf"), float("inf"))

    # -----------------------------------------------------------------

    @classmethod
    def zero(cls):
        return cls(-0., 0.)

    # -----------------------------------------------------------------

    @Range.min.setter
    def min(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Check
        if types.is_integer_type(value): value = float(value)
        if not types.is_real_type(value): raise ValueError("Value must be real")

        # Set
        self._min = value

    # -----------------------------------------------------------------

    @Range.max.setter
    def max(self, value):

        """
        Thisn function ...
        :param value:
        :return:
        """

        # Check
        if types.is_integer_type(value): value = float(value)
        if not types.is_real_type(value): raise ValueError("Value must be real")

        # Set
        self._max = value

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This fucntion ...
        :param value:
        :return:
        """

        from ..tools import types

        if types.is_quantity(value): return QuantityRange(self.min * value, self.max * value)
        elif types.is_real_type(value): return self.__class__(self.min * value, self.max * value)
        elif types.is_integer_type(value): return self.__class__(self.min * value, self.max * value)
        else: raise ValueError("Value must be Quantity, real or integer value")

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from ..tools import types

        if types.is_quantity(value): return QuantityRange(self.min / value, self.max / value)
        elif types.is_real_type(value) or types.is_integer_type(value): return self.__class__(self.min / float(value), self.max / float(value))
        else: raise ValueError("Value must be Quantity, real or integer value")

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

        else: raise ValueError("min_value and max_value must be either both quantities or both floats (with unit specified separately)")

        # Call the constructor of the base class
        super(QuantityRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert, rearrange=rearrange)

        # Set the unit
        self.unit = unit

    # -----------------------------------------------------------------

    @property
    def min(self):
        return self._min * self.unit

    # -----------------------------------------------------------------

    @property
    def max(self):
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

    def __mul__(self, value):

        """
        This fucntion ...
        :param value:
        :return:
        """

        from ..tools import types

        if types.is_quantity(value): return QuantityRange(self.min * value, self.max * value)
        elif types.is_real_type(value): return self.__class__(self.min * value, self.max * value)
        elif types.is_integer_type(value): return self.__class__(self.min * value, self.max * value)
        else: raise ValueError("Value must be Quantity, real or integer value")

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from ..tools import types

        if types.is_quantity(value): return QuantityRange(self.min / value, self.max / value)
        elif types.is_real_type(value) or types.is_integer_type(value): return self.__class__(self.min / float(value), self.max / float(value))
        else: raise ValueError("Value must be Quantity, real or integer value")

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

def round_to_1(x):

    """
    This function ...
    :param x:
    :return:
    """

    return round(x, -int(math.floor(math.log10(abs(x)))))

# -----------------------------------------------------------------
