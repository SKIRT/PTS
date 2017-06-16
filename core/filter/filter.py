#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.filter.filter Contains the Filter class.

# -----------------------------------------------------------------

# Import standard modules
from abc import ABCMeta, abstractproperty

# Import the relevant PTS classes and modules
from ..tools import types

# -----------------------------------------------------------------

def parse_filter(argument, name=None):

    """
    This function ...
    :param argument:
    :param name:
    :return:
    """

    # Import subclasses
    from .narrow import NarrowBandFilter
    from .broad import BroadBandFilter

    # If the argument that is passed is already a Filter instance
    if isinstance(argument, Filter): return argument

    try: fltr = BroadBandFilter(argument, name=name)
    except ValueError:
        try: fltr = NarrowBandFilter(argument, name=name)
        except ValueError: raise ValueError("Could not parse the string '" + argument + "' as a filter")

    # Return the filter
    return fltr

# -----------------------------------------------------------------

def represent_filter(fltr, delimiter=" "):

    """
    This function ...
    :param fltr:
    :param delimiter:
    :return:
    """

    return str(fltr).replace(" ", delimiter)

# -----------------------------------------------------------------

def stringify_filter(fltr, delimiter=" "):

    """
    This function ...
    :param fltr:
    :param delimiter:
    :return:
    """

    # Import subclasses
    from .narrow import NarrowBandFilter
    from .broad import BroadBandFilter

    #print(delimiter)

    # Stringify
    if isinstance(fltr, NarrowBandFilter): return "narrow_band_filter", represent_filter(fltr, delimiter=delimiter)
    elif isinstance(fltr, BroadBandFilter): return "broad_band_filter", represent_filter(fltr, delimiter=delimiter)
    else: raise ValueError("Invalid argument: must be narrow or broad band filter")

# -----------------------------------------------------------------

class Filter(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, filter_id, description):

        """
        This function ...
        :param filter_id:
        :param description:
        """

        # Set attributes
        self._FilterID = filter_id
        self._Description = description

    # -----------------------------------------------------------------

    def filterID(self):

        """
        This fucntion returns a unique identifer for the filter
        :return:
        """

        return self._FilterID

    # -----------------------------------------------------------------

    def description(self):

        """
        This function returns a human-readable description for the filter.
        :return:
        """

        return self._Description

    # -----------------------------------------------------------------

    @abstractproperty
    def __repr__(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def name(self):

        """
        This property ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def observatory(self):

        """
        This property ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def instrument(self):

        """
        This property ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def band(self):

        """
        This property ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def mean(self):

        """
        This property returns the mean wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def effective(self):

        """
        This property returns the effective wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def min(self):

        """
        This function returns the minimum wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def max(self):

        """
        This fucntion returns the maximum wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def center(self):

        """
        This property returns the center wavelength
        :return:
        """

    # -----------------------------------------------------------------

    @abstractproperty
    def pivot(self):

        """
        This property returns the pivot wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def range(self):

        """
        This property returns the wavelength range
        :return:
        """

        from ..basics.range import QuantityRange
        return QuantityRange(self.min, self.max)

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function compares two filter instances
        :param other:
        :return:
        """

        if types.is_string_type(other): other = parse_filter(other)
        return str(other) == str(self)

    # -----------------------------------------------------------------

    def __ne__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        if types.is_string_type(other): other = parse_filter(other)
        return str(other) != str(self)

    # -----------------------------------------------------------------

    def __hash__(self):

        """
        This function ...
        :return: 
        """

        return hash(str(self))

# -----------------------------------------------------------------
