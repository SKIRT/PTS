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
import numpy as np

# Import astronomical modules
from astropy.units import Unit, Quantity

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

    if isinstance(value, int): return IntegerRange.around(value, rel_min, rel_max)
    elif isinstance(value, float): return RealRange.around(value, rel_min, rel_max)
    elif isinstance(value, Quantity): return QuantityRange.around(value, rel_min, rel_max)
    else: raise ValueError("Value has unknown type '" + str(type(value)) + "'")

# -----------------------------------------------------------------

class Range(object):
    
    """
    This class ...
    """
    
    def __init__(self, min_value, max_value, inclusive=True, invert=False):
        
        """
        The constructor ...
        :param min_value:
        :param max_value:
        :param inclusive:
        """

        self.min = min_value
        self.max = max_value
        self.inclusive = inclusive
        self.invert = invert

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

    @property
    def log_min(self):

        """
        This function ...
        :return:
        """

        return np.log10(self.min)

    # -----------------------------------------------------------------

    @property
    def log_max(self):

        """
        This function ...
        :return:
        """

        return np.log10(self.max)

    # -----------------------------------------------------------------

    def linear(self, npoints, as_list=False):

        """
        This function ...
        :param as_list:
        """

        values = np.linspace(self.min, self.max, npoints, endpoint=self.inclusive)
        if self.invert: values = np.flipud(values)

        if as_list: return list(values)
        else: return values

    # -----------------------------------------------------------------

    def log(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        values = np.logspace(self.log_min, self.log_max, npoints, endpoint=self.inclusive)
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

        width = self.max - self.min
        normalized = np.linspace(0.0, 1.0, npoints, endpoint=self.inclusive)
        values = self.min + normalized * width
        if self.invert: values = np.flipud(values)

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

class IntegerRange(Range):

    """
    This class ...
    """

    def __init__(self, min_value, max_value, inclusive=True, invert=False):

        """
        This function ...
        :param min_value:
        :param max_value:
        :param inclusive:
        :param invert:
        """

        assert isinstance(min_value, int)
        assert isinstance(max_value, int)

        # Call the constructor of the base class
        super(IntegerRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert)

    # -----------------------------------------------------------------

    def linear(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(IntegerRange, self).linear(npoints)
        integers = list(set(map(int, real)))

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
        integers = list(set(map(int, real)))

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
        integers = list(set(map(int, real)))

        if as_list: return integers
        else: return np.array(integers)

# -----------------------------------------------------------------

class RealRange(Range):

    """
    This class ...
    """

    def __init__(self, min_value, max_value, inclusive=True, invert=False):

        """
        The constructor ...
        :param min_value:
        :param max_value:
        :param inclusive:
        """

        assert isinstance(min_value, float)
        assert isinstance(max_value, float)

        # Call the constructor of the base class
        super(RealRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert)

# -----------------------------------------------------------------

class QuantityRange(Range):

    """
    This function ...
    """

    def __init__(self, min_value, max_value, unit=None, inclusive=True, invert=False):

        """
        This function ...
        :param min_value:
        :param max_value:
        :param unit:
        :param inclusive:
        :param invert:
        """

        # Convert everything so that min_value and max_value are floats in the same unit, and so that 'unit' is the corresponding Unit
        min_is_quantity = hasattr(min_value, "unit")
        max_is_quantity = hasattr(max_value, "unit")

        if min_is_quantity and max_is_quantity:

            unit = min_value.unit
            min_value = min_value.value
            max_value = max_value.to(unit).value

        elif (not min_is_quantity) and (not max_is_quantity):

            if unit is None: raise ValueError("Unit must be specified if min_value and max_value are not quantities")
            elif isinstance(unit, basestring): unit = Unit(unit)

        else: raise ValueError("min_value and max_value must be either both quantities or both floats (with unit specified seperately)")

        # Call the constructor of the base class
        super(QuantityRange, self).__init__(min_value, max_value, inclusive=inclusive, invert=invert)

        # Set the unit
        self.unit = unit

    # -----------------------------------------------------------------

    def linear(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(QuantityRange, self).linear(npoints)

        if as_list:
            result = []
            for num in real: result.append(num * self.unit)
        else: result = real * self.unit

        # Return the result (list or quantity)
        return result

    # -----------------------------------------------------------------

    def log(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(QuantityRange, self).log(npoints)

        if as_list:
            result = []
            for num in real: result.append(num * self.unit)
        else: result = real * self.unit

        # Return the result (list or quantity)
        return result

    # -----------------------------------------------------------------

    def sqrt(self, npoints, as_list=False):

        """
        This function ...
        :param npoints:
        :param as_list:
        :return:
        """

        real = super(QuantityRange, self).sqrt(npoints)

        if as_list:
            result = []
            for num in real: result.append(num * self.unit)
        else: result = real * self.unit

        # Return the result (list or quantity)
        return result

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
