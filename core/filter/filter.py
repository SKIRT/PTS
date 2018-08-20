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

# Import astronomical modules
from astropy.units import spectral

# Import the relevant PTS classes and modules
from ..tools import types
from ..tools.utils import lazyproperty

# -----------------------------------------------------------------

def parse_filter_from_instrument_and_band(instrument, band):

    """
    Thisn function ...
    :param instrument:
    :param band:
    :return:
    """

    # Import subclasses
    from .narrow import NarrowBandFilter
    from .broad import BroadBandFilter

    # Parse
    try: fltr = BroadBandFilter.from_instrument_and_band(instrument, band)
    except ValueError:
        try: fltr = NarrowBandFilter.from_instrument_and_band(instrument, band)
        except ValueError: raise ValueError("Could not find a filter for instrument '" + instrument + "' and band '" + band + "'")

    # Return the filter
    return fltr

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

    # Parse
    try: fltr = BroadBandFilter(argument, name=name)
    except ValueError:
        try: fltr = NarrowBandFilter(argument, name=name)
        except ValueError: raise ValueError("Could not parse the string '" + argument + "' as a filter")

    # Return the filter
    return fltr

# -----------------------------------------------------------------

def represent_filter(fltr, delimiter=" ", short=False):

    """
    This function ...
    :param fltr:
    :param delimiter:
    :param short:
    :return:
    """

    if short: return fltr.band
    else: return str(fltr).replace(" ", delimiter)

# -----------------------------------------------------------------

def stringify_filter(fltr, **kwargs):

    """
    This function ...
    :param fltr:
    :param kwargs:
    :return:
    """

    # Import subclasses
    from .narrow import NarrowBandFilter
    from .broad import BroadBandFilter

    # Get delimiter
    delimiter = kwargs.pop("delimiter", " ")

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

    @classmethod
    def from_instrument_and_band(cls, instrument, band):

        """
        This function ...
        :param instrument:
        :param band:
        :return:
        """

        return cls(instrument + "." + band)

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

    @lazyproperty
    def mean_frequency(self):
        return self.mean.to("Hz", equivalencies=spectral())

    # -----------------------------------------------------------------

    @abstractproperty
    def effective(self):

        """
        This property returns the effective wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def effective_frequency(self):
        return self.effective.to("Hz", equivalencies=spectral())

    # -----------------------------------------------------------------

    @abstractproperty
    def min(self):

        """
        This function returns the minimum wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def min_frequency(self):
        return self.max.to("Hz", equivalencies=spectral()) # MAX WAVELENGTH

    # -----------------------------------------------------------------

    @abstractproperty
    def max(self):

        """
        This fucntion returns the maximum wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def max_frequency(self):
        return self.min.to("Hz", equivalencies=spectral()) # MIN WAVELENGTH

    # -----------------------------------------------------------------

    @abstractproperty
    def center(self):

        """
        This property returns the center wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def center_frequency(self):
        return self.center.to("Hz", equivalencies=spectral())

    # -----------------------------------------------------------------

    @abstractproperty
    def pivot(self):

        """
        This property returns the pivot wavelength
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def pivot_frequency(self):
        return self.pivot.to("Hz", equivalencies=spectral())

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

    @lazyproperty
    def frequency_range(self):

        """
        This function ...
        :return:
        """

        from ..basics.range import QuantityRange
        return QuantityRange(self.min_frequency, self.max_frequency)

    # -----------------------------------------------------------------

    @lazyproperty
    def frequency(self):
        return self.wavelength.to("Hz", equivalencies=spectral())

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

def is_uv(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_uv(fltr.wavelength)

# -----------------------------------------------------------------

def is_optical(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_optical(fltr.wavelength)

# -----------------------------------------------------------------

def is_ir(fltr):

    """
    Thisj function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_ir(fltr.wavelength)

# -----------------------------------------------------------------

def is_nir(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_nir(fltr.wavelength)

# -----------------------------------------------------------------

def is_mir(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_mir(fltr.wavelength)

# -----------------------------------------------------------------

def is_fir(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_fir(fltr.wavelength)

# -----------------------------------------------------------------

def is_submm(fltr):

    """
    Thisn function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_submm(fltr.wavelength)

# -----------------------------------------------------------------

def is_fir_or_submm(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    from ...magic.tools import wavelengths
    return wavelengths.is_fir_or_submm(fltr.wavelength)

# -----------------------------------------------------------------
