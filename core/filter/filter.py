#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.filter.filter Contains the Filter class.

# -----------------------------------------------------------------

# Import standard modules
from abc import ABCMeta, abstractproperty

# -----------------------------------------------------------------

def parse_filter(argument, name=None):

    """
    This function ...
    :param argument:
    :param name:
    :return:
    """

    # import subclasses
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

class Filter(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, filter_id, description, true_filter):

        """
        This function ...
        :param filter_id:
        :param description:
        :param true_filter:
        """

        # Set attributes
        self._FilterID = filter_id
        self._Description = description
        self.true_filter = true_filter

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

    @property
    def name(self):

        """
        This property ...
        :return:
        """

        if self.true_filter: return self._FilterID.split("/")[1]
        else: return self._FilterID.split("_")[0]

    # -----------------------------------------------------------------

    @property
    def observatory(self):

        """
        This property ...
        :return:
        """

        if self.true_filter: return self._FilterID.split("/")[0]
        else: return None

    # -----------------------------------------------------------------

    @property
    def instrument(self):

        """
        This property ...
        :return:
        """

        if self.true_filter: return self._FilterID.split("/")[1].split(".")[0]
        else: return None

    # -----------------------------------------------------------------

    @property
    def band(self):

        """
        This property ...
        :return:
        """

        if self.true_filter: return self._FilterID.split("/")[1].split(".")[1].replace("_ext", "")
        elif self._FilterID.startswith("Uniform"): return None
        else: return self._FilterID.split("_")[0]

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

        from .range import QuantityRange
        return QuantityRange(self.min, self.max)

# -----------------------------------------------------------------

# -----------------------------------------------------------------
