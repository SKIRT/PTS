#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.errorbar Contains the ErrorBar class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# -----------------------------------------------------------------

def sum_errorbars_quadratically(*args):

    """
    This function ...
    :param args:
    :return:
    """

    lowers = []
    uppers = []

    # Loop over the given error bars
    for arg in args:

        lowers.append(abs(arg.lower))
        uppers.append(arg.upper)

    lowers = np.array(lowers)
    uppers = np.array(uppers)

    lower = - np.sqrt(np.sum(lowers**2))
    upper = np.sqrt(np.sum(uppers**2))

    # Set and return the new error bar
    return ErrorBar(lower, upper)

# -----------------------------------------------------------------

class ErrorBar(object):
    
    """
    This class ...
    """
    
    def __init__(self, lower, upper=None, at=None):
        
        """
        The constructor ...
        """

        # If only a single value is specified for the errorbar, assume both lower and upper are equal
        if upper is None:

            # 'lower' should be positive (first argument in this case is actually the 'upper' ...)
            if at is None:

                if lower < 0: raise ValueError("When passing only one value (symmetric errorbar), this value should be positive")

                self.lower = - lower
                self.upper = lower

            # 'lower' value minus 'at' value should be positive (first argument in this case is actually the 'upper'...)
            else:

                test = lower - at
                if test < 0: raise ValueError("When passing only one value (symmetric errorbar), this value should be positive")

                self.lower = - test
                self.upper = test

        # A lower and upper value is given
        else:

            # The 'at' value is not specified, assume lower and upper are
            if at is None:

                # Check that lower is negative and upper is positive (but zero is allowed for both)
                if lower > 0: raise ValueError("Lower limit should be negative or zero: " + str(lower))
                if upper < 0: raise ValueError("Upper limit should be positive or zero: " + str(upper))

                self.lower = lower
                self.upper = upper

            else:

                # Check that lower is smaller than or equal to 'at' and that upper is greater than or equal to 'at'
                if lower > at: raise ValueError("Lower limit should be smaller than or equal to 'at' value: " + str(lower))
                if upper < at: raise ValueError("Upper limit should be greater than or equal to 'at' value: " + str(upper))

                self.lower = lower - at
                self.upper = upper - at

        # Final sanity check
        assert self.lower <= 0 and self.upper >= 0, "Lower="+str(self.lower)+", upper="+str(self.upper)+" (input: "+str(lower)+","+str(upper)+","+str(at)+")"

    # -----------------------------------------------------------------

    @property
    def unit_lower(self):
        if hasattr(self.lower, "unit"): return self.lower.unit
        else: return None

    # -----------------------------------------------------------------

    @property
    def unit_upper(self):
        if hasattr(self.upper, "unit"): return self.upper.unit
        else: return None

    # -----------------------------------------------------------------

    @property
    def lower_has_unit(self):
        return self.unit_lower is not None

    # -----------------------------------------------------------------

    @property
    def upper_has_unit(self):
        return self.unit_upper is not None

    # -----------------------------------------------------------------

    @property
    def lower_scalar(self):
        if self.lower_has_unit: return self.lower.value
        else: return self.lower

    # -----------------------------------------------------------------

    @property
    def upper_scalar(self):
        if self.upper_has_unit: return self.upper.value
        else: return self.upper

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        # No unit
        if self.unit_lower is None and self.unit_upper is None: return None

        # Check
        if self.unit_lower != self.unit_upper: raise ValueError("Units of lower and upper limit are different")
        return self.unit_lower

    # -----------------------------------------------------------------

    @property
    def has_unit(self):
        return self.unit is not None

    # -----------------------------------------------------------------

    @classmethod
    def zero(cls):

        """
        Thisn function ...
        :return:
        """

        return cls(0.0)

    # -----------------------------------------------------------------

    def to(self, unit):

        """
        This function ...
        :param unit:
        :return:
        """

        return ErrorBar(self.lower.to(unit), self.upper.to(unit))

    # -----------------------------------------------------------------

    def as_tuple(self, scalar=False):

        """
        This function ...
        :param scalar:
        :return:
        """

        if scalar: return (self.lower_scalar, self.upper_scalar,)
        else: return (self.lower, self.upper,)

    # -----------------------------------------------------------------

    @property
    def average(self):

        """
        This function ...
        :return:
        """

        value = 0.5 * (abs(self.lower_scalar) + self.upper_scalar)
        if self.has_unit: return value * self.unit
        else: return value

    # -----------------------------------------------------------------

    @property
    def has_upper(self):

        """
        This function ...
        :return:
        """

        from ..tools import numbers
        return self.upper != numbers.inf

    # -----------------------------------------------------------------

    @property
    def has_lower(self):

        """
        This function ...
        :return:
        """

        from ..tools import numbers
        return self.lower != numbers.min_inf

    # -----------------------------------------------------------------

    @property
    def no_upper(self):

        """
        This function ...
        :return:
        """

        return not self.has_upper

    # -----------------------------------------------------------------

    @property
    def no_lower(self):

        """
        This function ...
        :return:
        """

        return not self.has_lower

    # -----------------------------------------------------------------

    @property
    def only_upper(self):

        """
        This function ...
        :return:
        """

        return self.has_upper and self.no_lower

    # -----------------------------------------------------------------

    @property
    def only_lower(self):

        """
        This function ...
        :return:
        """

        return self.has_lower and self.no_upper

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return ErrorBar(self.lower * value, self.upper * value)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.lower *= value
        self.upper *= value

        # Return a reference to this object
        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return ErrorBar(self.lower / value, self.upper / value)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.lower /= value
        self.upper /= value

        # Return a reference to this object
        return self

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __itruediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__idiv__(value)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        from ..tools import stringify
        return stringify.stringify(self.lower)[1] + ">" + stringify.stringify(self.upper)[1]

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function
        """

        return '<' + self.__class__.__name__ + ' -{0}, +{1}>'.format(abs(self.lower), abs(self.lower))

# -----------------------------------------------------------------
