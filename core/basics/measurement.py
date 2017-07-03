#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.measurement Contains the Measurement class,
#  representing floating point values with a certain uncertainty

# -----------------------------------------------------------------

# Import standard modules
import math

# -----------------------------------------------------------------

class Measurement(object):
    
    """
    This class ...
    """
    
    def __init__(self, value, error=None):
        
        """
        The constructor ...
        """
        
        # Set the attributes
        self.value = value
        self.error = error

    # -----------------------------------------------------------------

    @property
    def relative_error(self):

        """
        This function ...
        :return:
        """

        try: return self.error / self.value
        except ZeroDivisionError: return float("nan")

    # -----------------------------------------------------------------

    def __add__(self, quantity):

        """
        This function ...
        :param quantity:
        :return:
        """

        value = self.value + quantity.value
        error = math.sqrt(math.pow(self.error, 2) + math.pow(quantity.error, 2))
        return Measurement(value, error)

    # -----------------------------------------------------------------

    def __sub__(self, quantity):

        """
        This function ...
        :param quantity:
        :return:
        """

        value = self.value - quantity.value
        error = math.sqrt(math.pow(self.error, 2) + math.pow(quantity.error, 2))
        return Measurement(value, error)

    # -----------------------------------------------------------------

    def __mul__(self, quantity):

        """
        This function ...
        :param quantity:
        :return:
        """

        value = self.value * quantity.value
        error = math.sqrt(math.pow(quantity.value * self.error, 2) + math.pow(self.value * quantity.error, 2))
        return Measurement(value, error)

    # -----------------------------------------------------------------

    def __div__(self, quantity):

        """
        This function ...
        :param quantity:
        :return:
        """

        value = self.value / quantity.value
        error = math.fabs(value) * math.sqrt(math.pow(self.relative_error, 2) + math.pow(quantity.relative_error, 2))
        return Measurement(value, error)

    # -----------------------------------------------------------------

    def __truediv__(self, quantity):

        """
        This function ...
        :param quantity:
        :return:
        """

        value = self.value / quantity.value
        error = math.fabs(value) * math.sqrt(math.pow(self.relative_error, 2) + math.pow(quantity.relative_error, 2))
        return Measurement(value, error)

# -----------------------------------------------------------------
